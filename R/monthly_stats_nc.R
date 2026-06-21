#' @title monthly_stats_nc
#'
#' @description Computes calendar-month statistics on gridded data — either an
#'   sxts object, a NetCDF file, or a raster object. For each calendar month
#'   present in the dataset, the function calls \code{\link{basic_stats_nc}},
#'   returning a multi-layer raster of per-cell statistics: mean, min, max,
#'   standard deviation, L-moments, probability dry, empirical quantiles, and
#'   (unless \code{ignore_zeros = TRUE}) wet/dry transition statistics. The
#'   outer loop over up to twelve calendar months is separable from the
#'   already-vectorised per-cell computations inside \code{basic_stats_nc}, so
#'   parallelisation is offered at the month level only. NetCDF input is loaded
#'   once via \code{\link{nc2xts}} and then sliced by month index; raster input
#'   is converted to sxts with automatic date parsing.
#'
#' @param data An sxts object, or a Raster* object. Leave \code{NULL} when
#'   supplying \code{filename} and \code{varname}.
#' @param filename Optional; path to a NetCDF file to import when \code{data}
#'   is not provided.
#' @param varname Optional; name of the variable to extract from
#'   \code{filename}.
#' @param ignore_zeros Logical; if \code{TRUE}, zeros are ignored in the
#'   per-cell statistics. Default \code{FALSE}.
#' @param zero_threshold Numeric; threshold below which values are treated as
#'   zero. Default \code{0.01}.
#' @param parallel Logical; whether to compute per-month statistics in parallel.
#'   Default \code{FALSE}.
#' @param ncores Integer; number of cores for parallel computation. Default
#'   \code{2}.
#' @param ... Additional arguments passed to \code{\link{nc2xts}} when
#'   \code{filename} and \code{varname} are supplied.
#'
#' @return A named list with one element per calendar month present in the data,
#'   named after \code{\link[base]{month.name}}. Each element is the
#'   \code{\link{basic_stats_nc}} output raster for that month.
#'
#' @examples
#' # Synthetic sxts with 4 cells and 2 years of daily data
#' set.seed(123)
#' n <- 365 * 2
#' dates <- seq(as.POSIXct("2000-01-01", tz = "UTC"), by = "day", length.out = n)
#' vals <- matrix(pmax(0, rnorm(n * 4, mean = 3, sd = 5)), nrow = n, ncol = 4)
#' coords <- data.frame(x = c(0, 1, 0, 1), y = c(50, 50, 51, 51))
#' ts_sxts <- sxts(data = vals, order.by = dates, coords = coords,
#'                 projection = "+proj=longlat +datum=WGS84")
#'
#' ms_nc <- monthly_stats_nc(data = ts_sxts)
#' names(ms_nc)
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
