#' @title monthly_stats_nc
#'
#' @description This function calculates basic statistics on a monthly basis for
#' gridded data (a NetCDF raster file or an sxts object). For each calendar month
#' present in the data it returns the same output as \code{\link{basic_stats_nc}}:
#' a multi-layer raster of per-cell statistics (mean, min, max, SD, L-moments,
#' probability dry, quantiles, and, unless \code{ignore_zeros = TRUE}, wet/dry
#' transition statistics).
#'
#' Computation can optionally be parallelized across the calendar months (a
#' coarse outer loop of at most twelve tasks). The per-cell statistics are
#' already vectorized inside \code{basic_stats_nc}, so no cell-axis parallelism
#' is offered here.
#'
#' @param data An sxts object, or a raster file. Leave NULL when supplying filename/varname.
#' @param filename (optional) A NetCDF file name to import if data is not provided.
#' @param varname (optional) The name of the variable to extract from 'filename'.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param parallel Logical, whether to compute the per-month statistics in
#'   parallel across calendar months.
#' @param ncores Number of cores to use in the case of parallel computations.
#' @param ... Additional arguments to pass to 'nc2xts' function (if 'filename' and 'varname' are provided).
#'
#' @return A named list with one element per calendar month present in the data
#' (named after \code{month.name}). Each element is the \code{\link{basic_stats_nc}}
#' raster of per-cell statistics for that month.
#'
#' @examples TO BE FILLED
#'
#' @importFrom lubridate month parse_date_time
#' @importFrom xts xts
#'
#' @export
#'
monthly_stats_nc = function(data = NULL, filename = NA, varname = NA,
                            ignore_zeros = FALSE, zero_threshold = 0.01,
                            parallel = FALSE, ncores = 2, ...){

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

  # Pre-subset each month in the master so a worker serializes only that month's
  # slice (~1/12 of rows), not the whole grid. Row-subsetting preserves coords.
  monthly_subsets <- lapply(i_months, function(m) ncdf_sxts[which(months_idx == m), ])

  # Worker capturing only small scalars (never the grid) - closure hygiene.
  stats_month <- function(sub)
    basic_stats_nc(data = sub, ignore_zeros = ignore_zeros,
                   zero_threshold = zero_threshold)

  # ---- Dispatch over months -----------------------------------------------
  if (parallel) {
    # Respect a user-set plan; otherwise default per-OS, sized by the month count,
    # and restore on exit (same pattern as fitlm_monthly_nc).
    if (inherits(future::plan(), "sequential")) {
      strategy <- if (.Platform$OS.type == "windows") future::multisession else future::multicore
      oplan <- future::plan(strategy, workers = min(ncores, length(i_months)))
      on.exit(future::plan(oplan), add = TRUE)
    }
    oopt <- options(future.globals.maxSize = Inf); on.exit(options(oopt), add = TRUE)
    monthly_results <- future.apply::future_lapply(monthly_subsets, stats_month,
                         future.seed = TRUE, future.packages = "anyFit")
  } else {
    monthly_results <- lapply(monthly_subsets, stats_month)
  }

  names(monthly_results) <- month.name[i_months]
  return(monthly_results)
}
