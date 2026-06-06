#' @title period_apply_nc
#'
#' @description This function applies a summary function (e.g., mean, sum) to a NetCDF raster file
#' over specified time periods (e.g., months, years).
#'
#' @param raster_file A raster file
#' @param filename (optional) A NetCDF file name to import if raster_file is not provided.
#' @param varname (optional) The name of the variable to extract from 'filename'.
#' @param period The time period to apply the function (e.g., "months", "years").
#' @param period_multiplier a multiplier to define arbitrary periods based on the period argument
#' @param FUN The summary function to apply (e.g., "mean", "sum").
#' @param ... Additional arguments to pass to 'nc2xts' function (if 'filename' and 'varname' are provided).
#'
#' @return A raster object with the summary statistics for the specified time periods.
#'
#' @examples
#' # Load required libraries
#' library(raster)
#' library(xts)
#'
#' # Example usage 1: Using a raster file directly
#' result <- period_apply_nc(raster, period = "months", FUN = "mean")
#' nc_ggplot(result)
#'
#' # Example usage 2: Using a filename and variable name
#' nc_file <-"data/ncfile.nc"
#' result <- period_apply_nc(filename = nc_file, varname = "tp", period = "months", FUN = "sum")
#' nc_ggplot(result)
#'
#' @import raster
#' @import xts
#' @importFrom xts endpoints
#' @importFrom lubridate ymd ydm mdy myd dmy dym ymd_h dmy_h mdy_h ydm_h ymd_hm dmy_hm mdy_hm ydm_hm ymd_hms dmy_hms mdy_hms ydm_hms
#'

#' @export
period_apply_nc = function (data, filename = NA, varname = NA, period = "months",
                            period_multiplier = 1, FUN = "mean", ...){
  if (!is.na(filename)) {
    temp = nc2xts(filename = filename, varname = varname,
                  ...)
    ncdf_sxts = temp$ncdf_sxts
  }
  else if ("raster" %in% class(data)){
    ncdf_sxts <- rasterFromSxts.sxts(data)
  } else {
    ncdf_sxts <- data
  }

  spec_period <- endpoints(ncdf_sxts, on = period, k = period_multiplier)
  ncdf_stats = lapply(1:ncol(ncdf_sxts), FUN = function(x) {
    period.apply(ncdf_sxts[, x], spec_period, FUN = FUN)
  })
  ncdf_stats = do.call(cbind, ncdf_stats)
  # ncdf_stats = cbind(coords, ncdf_stats)
  # raster_fun = raster::rasterFromXYZ(ncdf_stats)
  # projection(raster_fun) = raster::projection(raster_file)
  return(ncdf_stats)
}


