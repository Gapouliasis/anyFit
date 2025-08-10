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
#' @param raster_file A raster file
#' @param filename (optional) A NetCDF file name to import if raster_file is not provided.
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
#' @import raster
#' @import xts
#' @importFrom lubridate ymd ydm mdy myd dmy dym ymd_h dmy_h mdy_h ydm_h ymd_hm dmy_hm mdy_hm ydm_hm ymd_hms dmy_hms mdy_hms ydm_hms
#' @importFrom parallel mclapply
#'
#' @export
#'
basic_stats_nc = function(raster_file, filename = NA, varname = NA,
                    ignore_zeros = FALSE, zero_threshold = 0.01 ,
                    parallel = FALSE, ncores = 2, ...){

  if (!is.na(filename)){
    temp = nc2xts(filename = filename, varname = varname,...)
    raster_file = temp$raster
    ncdf_xts = temp$ncdf_xts
  }else if (is.list(eobs_data) & all(c("ncdf_xts", "coordinates") %in% names(eobs_data))){
    ncdf_xts = eobs_data$ncdf_xts
    coords = eobs_data$coordinates
    dates = eobs_data$dates
  }else{
    t = raster::rasterToPoints(raster_file)
    tt = t(t)
    coords = t(tt[c(1,2),])
    tt = tt[-c(1,2),]
    dates = rownames(tt)
    temp_dates = gsub("X",replacement = "",x=dates)
    funs = c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "ymd H", "dmy H", "mdy H",
             "ydm H", "ymd HM", "dmy HM", "mdy HM", "ydm HM", "ymd HMS", "dmy HMS",
             "mdy HMS", "ydm HMS")

    dates=parse_date_time(temp_dates, orders = funs)
    ncdf_xts = xts(x = tt,order.by = dates)
  }

  if(Sys.info()['sysname'] == "Windows" & parallel){
    ncdf_stats = parallelsugar::mclapply(1:ncol(ncdf_xts), FUN = function(x){basic_stats(ncdf_xts[,x])$stats_table}, mc.cores = ncores)
  }else if(parallel){
    ncdf_stats = parallel::mclapply(1:ncol(ncdf_xts), FUN = function(x){basic_stats(ncdf_xts[,x])$stats_table}, mc.cores = ncores)
  }else{
    ncdf_stats = lapply(1:ncol(ncdf_xts), FUN = function(x){basic_stats(ncdf_xts[,x])$stats_table})
  }
  ncdf_stats = t(do.call(cbind,ncdf_stats))
  ncdf_stats = cbind(coords,ncdf_stats)
  #rownames(ncdf_stats) = as.character(dates)
  raster_stats = raster::rasterFromXYZ(ncdf_stats)
  projection(raster_stats) = projection(raster_file)

  return(raster_stats)
}
