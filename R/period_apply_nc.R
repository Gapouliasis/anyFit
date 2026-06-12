#' @title period_apply_nc
#'
#' @description This function applies a summary function (e.g., mean, sum) to spatiotemporal data
#' over specified time periods (e.g., months, years). Accepts an sxts object, a raster file,
#' or a NetCDF file path. When the input is an sxts object the output preserves spatial
#' metadata (coordinates, projection) and is returned as an sxts object.
#'
#' @param data An sxts object, or a raster file. Leave NULL when supplying filename/varname.
#' @param filename (optional) A NetCDF file name to import if data is not provided.
#' @param varname (optional) The name of the variable to extract from 'filename'.
#' @param period The time period to apply the function (e.g., "months", "years").
#' @param period_multiplier a multiplier to define arbitrary periods based on the period argument
#' @param FUN The summary function to apply (e.g., "mean", "sum").
#' @param ... Additional arguments to pass to 'nc2xts' function (if 'filename' and 'varname' are provided).
#'
#' @return An sxts object (when input is sxts or NetCDF) or an xts object (when input is raster)
#'   with the summary statistics for the specified time periods.
#'
#' @examples
#' # Example usage 1: Using an sxts object directly
#' result <- period_apply_nc(ncdf_sxts, period = "months", FUN = "mean")
#'
#' # Example usage 2: Using a filename and variable name
#' nc_file <- "data/ncfile.nc"
#' result <- period_apply_nc(filename = nc_file, varname = "tp", period = "months", FUN = "sum")
#' nc_ggplot(result)
#'
#' @importFrom xts endpoints period.apply
#'
#' @export
period_apply_nc = function(data = NULL, filename = NA, varname = NA, period = "months",
                            period_multiplier = 1, FUN = "mean", ...) {

  if (!is.na(filename)) {
    temp      <- nc2xts(filename = filename, varname = varname, ...)
    ncdf_sxts <- temp$ncdf_sxts

  } else if (is.sxts(data)) {
    ncdf_sxts <- data

  } else {
    ncdf_sxts <- sxtsFromRaster(data)
  }

  spec_period <- endpoints(ncdf_sxts, on = period, k = period_multiplier)

  # For mean/sum use the column-wise functions, which apply to all columns in a
  # single pass and are much faster than looping column by column.
  if (identical(FUN, "mean") || identical(FUN, mean)) {
    ncdf_stats <- xts::period.apply(ncdf_sxts, spec_period,
                                    FUN = function(x) colMeans(x, na.rm = TRUE))
  } else if (identical(FUN, "sum") || identical(FUN, sum)) {
    ncdf_stats <- xts::period.apply(ncdf_sxts, spec_period,
                                    FUN = function(x) colSums(x, na.rm = TRUE))
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
