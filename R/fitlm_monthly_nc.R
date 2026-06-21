#' @title Monthly Distribution Fitting on Gridded Data
#'
#' @description Gridded counterpart of \code{\link{fitlm_monthly}} that fits a
#'   list of candidate distributions to every calendar month of a NetCDF raster
#'   or an \code{sxts} object via the L-moments method. The input is first
#'   normalised to an \code{sxts} object (accepting \code{sxts}, \code{Raster*},
#'   or NetCDF filename and variable name), and each calendar month present in
#'   the record is extracted as a sub-grid. Per-month fitting is dispatched to
#'   \code{\link{fitlm_nc}}, which delegates to \code{\link{fitlm_nxts}} for
#'   per-cell L-moment fitting and returns parameter rasters, theoretical and
#'   sample L-moment rasters, GoF rasters, and a density GoF comparison plot.
#'   Parallelism is controlled by \code{parallel_by}: \code{"cells"} parallelises
#'   the cell-level fits within each month,
#'   whereas \code{"months"} parallelises across months but is serial across
#'   cells. The \code{shared_memory} flag determines whether the grid is shared
#'   with workers via a file-backed \code{big.matrix} or
#'   serialised in column chunks (suitable for multi-node plans).
#'   This is reccomended for single machine usage and especially windows for efficiency 
#' and reduced RAM consumption.The returned
#'   list is named by month and contains the full \code{fitlm_nc} output for each.
#'
#' @param data An \code{sxts} object, or a \code{Raster*} object. Leave
#'   \code{NULL} when supplying \code{filename} and \code{varname}.
#' @param filename A NetCDF file name to import if \code{data} is not provided.
#' @param varname The name of the variable to extract from \code{filename}.
#' @param candidates A character vector of distribution names to fit.
#' @param ignore_zeros A logical value, if \code{TRUE} zeros will be ignored.
#'   Default is \code{FALSE}.
#' @param zero_threshold The threshold below which values are considered zero.
#'   Default is 0.01.
#' @param parallel Logical, whether to use parallel processing.
#' @param ncores Number of cores to use for parallel computations. Default is 2.
#' @param shared_memory Logical, when parallel on the cell axis, share the grid
#'   with workers via a file-backed \code{big.matrix} (mmap, single machine)
#'   instead of serialising column chunks. Set \code{FALSE} for multi-node
#'   \code{plan(cluster)} setups. Default \code{TRUE}.
#' @param parallel_by Character, the axis to parallelise over when
#'   \code{parallel = TRUE}: \code{"cells"} (default) parallelises the grid-cell
#'   fits within each month via \code{fitlm_nxts}; \code{"months"} parallelises
#'   the per-month fits, each run serially across cells.
#' @param order Optional named list mapping a candidate name to the vector of
#'   L-moment orders matched by its optimiser, e.g.
#'   \code{list(gengamma = 1:5, expweibull = 1:3)}. Only the numerically-fitted
#'   distributions accept it; passed through to \code{\link{fitlm_nc}}.
#'   Default \code{NULL}.
#' @param ... Additional arguments passed to \code{\link{nc2xts}} when
#'   \code{filename} and \code{varname} are provided.
#'
#' @return A named list with one element per calendar month present in the data
#'   (named after \code{month.name}). Each element is the standard
#'   \code{\link{fitlm_nc}} output: a list with \code{fit_results} (per-candidate
#'   rasters) and \code{gof_plots}.
#'
#' @examples
#' \dontrun{
#' # Simulated 3-cell grid over two years
#' n <- 730
#' dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = n)
#' vals <- cbind(cell1 = rgamma(n, shape = 0.8, scale = 3),
#'               cell2 = rgamma(n, shape = 1.2, scale = 2),
#'               cell3 = rgamma(n, shape = 0.6, scale = 4))
#' coords <- data.frame(lon = c(10, 11, 12), lat = c(45, 46, 47))
#' grid <- sxts(vals, order.by = dates, coords = coords,
#'              projection = "+proj=longlat +datum=WGS84")
#'
#' monthly_fits <- fitlm_monthly_nc(grid, candidates = c("exp", "gamma3"),
#'                                   ignore_zeros = TRUE)
#' names(monthly_fits)
#' }
#'
#' @importFrom lubridate month parse_date_time
#'
#' @export
#'
fitlm_monthly_nc = function(data = NULL, filename = NA, varname = NA,
                            candidates = 'norm', ignore_zeros = FALSE, zero_threshold = 0.01,
                            parallel = FALSE, ncores = 2, shared_memory = TRUE,
                            parallel_by = c("cells", "months"), order = NULL, ...){

  parallel_by <- match.arg(parallel_by)

  # ---- Normalize any input to a single sxts (load once) --------------------
  if (!is.na(filename)) {
    ncdf_sxts <- nc2xts(filename = filename, varname = varname, ...)

  } else if (is.sxts(data)) {
    ncdf_sxts <- data

  } else {
    t         <- raster::rasterToPoints(data)
    tt        <- t(t)
    coords    <- t(tt[c(1, 2), ])
    tt        <- tt[-c(1, 2), ]
    dates     <- rownames(tt)
    temp_dates <- gsub("X", replacement = "", x = dates)
    funs <- c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "ymd H", "dmy H", "mdy H",
              "ydm H", "ymd HM", "dmy HM", "mdy HM", "ydm HM", "ymd HMS", "dmy HMS",
              "mdy HMS", "ydm HMS")
    dates     <- parse_date_time(temp_dates, orders = funs)
    proj_str  <- raster::projection(data)
    ncdf_sxts <- sxts(data = tt, order.by = dates, coords = as.data.frame(coords),
                      projection = proj_str)
    rm(t, tt)
  }

  # ---- Months present (handles missing months) ----------------------------
  months_idx <- lubridate::month(zoo::index(ncdf_sxts))
  i_months   <- sort(unique(months_idx))

  # Pre-subset each month in the master so the month-parallel path serializes
  # only each month's slice (~1/12 of rows) to a worker, not the whole grid 12x.
  monthly_subsets <- lapply(i_months, function(m) ncdf_sxts[which(months_idx == m), ])

  # Worker capturing only small scalars (never the grid) - closure hygiene as in
  # fitlm_nxts. `par_cells` toggles the cell-axis parallelism within fitlm_nc.
  fit_month <- function(sub, par_cells)
    fitlm_nc(data = sub, candidates = candidates,
             ignore_zeros = ignore_zeros, zero_threshold = zero_threshold,
             parallel = par_cells, ncores = ncores, shared_memory = shared_memory,
             order = order)

  # ---- Dispatch over months -----------------------------------------------
  if (parallel && parallel_by == "months") {
    # Respect a user-set plan; otherwise default per-OS, sized by the month count,
    # and restore on exit (same pattern as fitlm_nxts).
    if (inherits(future::plan(), "sequential")) {
      strategy <- if (.Platform$OS.type == "windows") future::multisession else future::multicore
      oplan <- future::plan(strategy, workers = min(ncores, length(i_months)))
      on.exit(future::plan(oplan), add = TRUE)
    }
    oopt <- options(future.globals.maxSize = Inf); on.exit(options(oopt), add = TRUE)
    monthly_results <- future.apply::future_lapply(monthly_subsets, fit_month,
                         par_cells = FALSE, future.seed = TRUE,
                         future.packages = "anyFit")
  } else {
    monthly_results <- lapply(monthly_subsets, fit_month,
                              par_cells = (parallel && parallel_by == "cells"))
  }

  names(monthly_results) <- month.name[i_months]
  return(monthly_results)
}
