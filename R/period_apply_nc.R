#' Temporal aggregation of gridded data
#'
#' @description
#' Aggregates gridded time series data to coarser temporal scales. Accepts an \code{sxts}
#' object, a \code{Raster*} object, or a NetCDF file path; when a filename is
#' supplied, data are imported through \code{\link{nc2xts}} first. The
#' aggregation period is controlled by \code{period} and
#' \code{period_multiplier}. For common
#' aggregation functions (\code{"mean"}, \code{"sum"}, \code{"min"}, \code{"max"},
#' \code{"median"}, \code{"var"}, \code{"sd"}), the function very efficient
#' column-wise function from \code{matrixStats} package, which ensures good computational
#' efficiency. Unsupported functions (e.g. custom functions) fall back to a per-column
#' lapply which is significantly slower. 
#'
#' @param data An \code{sxts} object, or a \code{Raster*} object. Leave
#'   \code{NULL} when supplying \code{filename} and \code{varname}.
#' @param filename Optional NetCDF file path to import if \code{data} is not provided.
#' @param varname Optional variable name to extract from \code{filename}.
#' @param period Period string passed to \code{xts::endpoints} (e.g.\code{"months"}, \code{"years"}).
#' @param period_multiplier Integer multiplier for custom period lengths (default 1).
#' @param FUN Summary function as a string (\code{"mean"}, \code{"sum"}, etc.) or a function.
#' @param ... Additional arguments passed to \code{\link{nc2xts}} when
#'   \code{filename} and \code{varname} are supplied.
#'
#' @return An \code{sxts} object when the input is \code{sxts} or NetCDF; an
#'   \code{xts} when the input is a \code{Raster*}.
#'
#' @examples
#' # Synthetic sxts
#' set.seed(123)
#' dates <- seq(as.POSIXct("2000-01-01"), as.POSIXct("2000-12-31"), by = "day")
#' coords <- data.frame(x = c(1, 2, 3, 4), y = c(1, 1, 1, 1))
#' vals <- matrix(rnorm(length(dates) * 4, 10, 5), nrow = length(dates), ncol = 4)
#' ncdf_sxts <- sxts(vals, order.by = dates, coords = coords,
#'                   projection = "+proj=longlat")
#'
#' # Monthly mean
#' result <- period_apply_nc(ncdf_sxts, period = "months", FUN = "mean")
#'
#' # Quarterly sum
#' result <- period_apply_nc(ncdf_sxts, period = "months",
#'                           period_multiplier = 3, FUN = "sum")
#'
#' @importFrom xts endpoints period.apply
#' @importFrom matrixStats colMeans2 colSums2 colMaxs colMins colMedians colVars colSds
#'
#' @export
period_apply_nc = function(data = NULL, filename = NA, varname = NA, period = "months",
                            period_multiplier = 1, FUN = "mean", ...) {

  if (!is.na(filename)) {
    ncdf_sxts <- nc2xts(filename = filename, varname = varname, ...)

  } else if (is.sxts(data)) {
    ncdf_sxts <- data

  } else {
    ncdf_sxts <- sxtsFromRaster(data)
  }

  spec_period <- endpoints(ncdf_sxts, on = period, k = period_multiplier)

  # For supported summaries use the column-wise matrixStats reducers, which apply
  # to all columns in a single pass and are much faster than looping column by
  # column. Anything else falls back to the per-column period.apply loop.
  reducer <- colwise_reducer(FUN)
  if (!is.null(reducer)) {
    ncdf_stats <- xts::period.apply(ncdf_sxts, spec_period, FUN = reducer)
    # matrixStats reducers return unnamed vectors, so reattach the column names.
    colnames(ncdf_stats) <- colnames(ncdf_sxts)
  } else {
    ncdf_stats <- do.call(cbind, lapply(1:ncol(ncdf_sxts), FUN = function(x) {
      period.apply(ncdf_sxts[, x], spec_period, FUN = FUN)
    }))
  }

  if (is.sxts(ncdf_sxts)) {
    ncdf_stats <- restore_sxts(ncdf_stats, ncdf_sxts)
  }

  return(ncdf_stats)
}

#' Map a summary function to a matrixStats column-wise reducer
#'
#' @description
#' Maps a summary function string (e.g. \code{"mean"}, \code{"sum"}) to the
#' equivalent \code{matrixStats} column-wise reducer, enabling single-pass
#' aggregation across all spatial columns. Returns \code{NULL} when no
#' column-wise equivalent exists, signalling the caller to use the generic
#' per-column fallback.
#'
#' @noRd
# Map a summary function (given as a string or the function itself) to a
# single-pass column-wise matrixStats reducer. Returns NULL when there is no
# column-wise equivalent, signalling the caller to use the generic fallback.
colwise_reducer <- function(FUN) {
  if (identical(FUN, "mean")   || identical(FUN, mean))
    return(function(x) matrixStats::colMeans2(x, na.rm = TRUE))
  if (identical(FUN, "sum")    || identical(FUN, sum))
    return(function(x) matrixStats::colSums2(x, na.rm = TRUE))
  if (identical(FUN, "max")    || identical(FUN, max))
    return(function(x) matrixStats::colMaxs(x, na.rm = TRUE))
  if (identical(FUN, "min")    || identical(FUN, min))
    return(function(x) matrixStats::colMins(x, na.rm = TRUE))
  if (identical(FUN, "median") || identical(FUN, stats::median))
    return(function(x) matrixStats::colMedians(x, na.rm = TRUE))
  if (identical(FUN, "var")    || identical(FUN, stats::var))
    return(function(x) matrixStats::colVars(x, na.rm = TRUE))
  if (identical(FUN, "sd")     || identical(FUN, stats::sd))
    return(function(x) matrixStats::colSds(x, na.rm = TRUE))
  NULL
}
