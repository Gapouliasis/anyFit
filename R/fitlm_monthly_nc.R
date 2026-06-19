#' @title fitlm_monthly_nc
#'
#' @description This function fits a list of candidate distributions on a monthly
#' basis to gridded data (a NetCDF raster file or an sxts object) using the
#' L-moments method. For each calendar month present in the data it returns the
#' same output as \code{\link{fitlm_nc}}: per-candidate rasters of the fitted
#' distribution parameters, the theoretical and sample L-moments, the
#' goodness-of-fit statistics, and a GoF density plot.
#'
#' Parallelism can be applied on either axis via \code{parallel_by}: across grid
#' cells (the \code{fitlm_nxts} engine, best for large grids) or across the
#' calendar months (a coarse outer loop of at most twelve tasks).
#'
#' @param data An sxts object, or a raster file. Leave NULL when supplying filename/varname.
#' @param filename (optional) A NetCDF file name to import if data is not provided.
#' @param varname (optional) The name of the variable to extract from 'filename'.
#' @param candidates A list of distribution to fit.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param parallel Logical, whether to use parallel processing.
#' @param ncores Number of cores to use in the case of parallel computations.
#' @param shared_memory Logical, when parallel on the cell axis, share the grid
#'   with workers via a filebacked big.matrix (mmap, single machine) instead of
#'   serializing column chunks. Set FALSE for multi-node \code{plan(cluster)}
#'   setups. Default TRUE.
#' @param parallel_by Character, the axis to parallelize over when
#'   \code{parallel = TRUE}: \code{"cells"} (default) parallelizes the grid-cell
#'   fits within each month via \code{fitlm_nxts}; \code{"months"} parallelizes
#'   the per-month fits, each run serially across cells.
#' @param order Optional named list mapping a candidate name to the vector of L-moment
#'   orders matched by its optimiser, e.g. \code{list(gengamma = 1:5, expweibull = 1:3)}.
#'   Only the numerically-fitted distributions accept it; passed through to
#'   \code{\link{fitlm_nc}}. Default NULL.
#' @param ... Additional arguments to pass to 'nc2xts' function (if 'filename' and 'varname' are provided).
#'
#' @return A named list with one element per calendar month present in the data
#' (named after \code{month.name}). Each element is the standard
#' \code{\link{fitlm_nc}} output: a list with \code{fit_results} (per-candidate
#' rasters) and \code{gof_plots}.
#'
#' @examples TO BE FILLED
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
