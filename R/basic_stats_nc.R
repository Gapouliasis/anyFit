#' @title basic_stats_nc
#'
#' @description This function calculates basic statistics (e.g., mean, min, max, SD) for a NetCDF raster file.
#' It returns the following statistics in raster format.
#' \itemize{
#' \item Number of Data - Number of missing data - Percentage of missing Data
#' \item Minimum - Maximum - Average - Variance - Coefficient of Variation - Standard Deviation - Third moment - Skewness
#' - Kurtosis
#' \item L-mean, L-scale , 3rd and 4th L-coefficients
#' \item Probability Dry
#' \item Quantiles 5th , 25th, 50th, 75th, 95th - Interquartile range.
#' \item Mean value and variance from a dry (defined by zero_threshold argument) to a wet state
#' \item Mean value and variance from a wet to a dry (defined by zero_threshold argument) state
#' \item Mean value and variance from a wet to a wet state
#' \item Transition probability from a wet to a wet state
#' \item Transition probability from a dry to a dry state
#' }
#'
#' @param data An sxts object, or a raster file. Leave NULL when supplying filename/varname.
#' @param filename (optional) A NetCDF file name to import if data is not provided.
#' @param varname (optional) The name of the variable to extract from 'filename'.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param parallel Logical, whether to use parallel processing.
#' @param ncores Number of cores to use in the case of parallel computations
#' @param ... Additional arguments to pass to 'nc2xts' function (if 'filename' and 'varname' are provided).
#'
#' @return A raster object containing basic statistics.
#'
#' @examples TO BE FILLED
#'
#' @importFrom lubridate parse_date_time
#' @importFrom parallel mclapply
#' @importFrom xts xts
#'
#' @export
#'
basic_stats_nc = function(data = NULL, filename = NA, varname = NA,
                    ignore_zeros = FALSE, zero_threshold = 0.01,
                    parallel = FALSE, ncores = 2, ...){

  if (!is.na(filename)) {
    temp      <- nc2xts(filename = filename, varname = varname, ...)
    ncdf_sxts <- temp$ncdf_sxts
    coords    <- coords(ncdf_sxts)
    proj_str  <- projection(ncdf_sxts)

  } else if (is.sxts(data)) {
    ncdf_sxts <- data
    coords    <- coords(ncdf_sxts)
    proj_str  <- projection(ncdf_sxts)

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
    ncdf_sxts <- xts(x = tt, order.by = dates)
    proj_str  <- raster::projection(data)
  }

  if (Sys.info()['sysname'] == "Windows" & parallel) {
    ncdf_stats <- parallelsugar::mclapply(1:ncol(ncdf_sxts), FUN = function(x) { basic_stats(ncdf_sxts[, x])$stats_table }, mc.cores = ncores)
  } else if (parallel) {
    ncdf_stats <- parallel::mclapply(1:ncol(ncdf_sxts), FUN = function(x) { basic_stats(ncdf_sxts[, x])$stats_table }, mc.cores = ncores)
  } else {
    ncdf_stats <- lapply(1:ncol(ncdf_sxts), FUN = function(x) { basic_stats(ncdf_sxts[, x])$stats_table })
  }

  ncdf_stats  <- t(do.call(cbind, ncdf_stats))
  ncdf_stats  <- cbind(coords, ncdf_stats)
  raster_stats <- raster::rasterFromXYZ(ncdf_stats)
  raster::projection(raster_stats) <- proj_str

  return(raster_stats)
}
