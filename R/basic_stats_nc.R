#' @title basic_stats_nc
#'
#' @description This function calculates basic statistics (e.g., mean, min, max, SD) for a NetCDF raster file.
#' It returns the statistics in raster format.
#'
#' @param raster_file A raster file
#' @param filename (optional) A NetCDF file name to import if raster_file is not provided.
#' @param varname (optional) The name of the variable to extract from 'filename'.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param parallel Logical, whether to use parallel processing.
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
                    parallel = FALSE,...){

  if (!is.na(filename)){
    temp = nc2xts(filename = filename, varname = varname,...)
    raster_file = temp$raster
    ncdf_xts = temp$ncdf_xts
  }else{
    t = raster::rasterToPoints(raster_file)
    tt = t(t)
    coords = t(tt[c(1,2),])
    tt = tt[-c(1,2),]
    dates = rownames(tt)
    temp_dates = gsub("X",replacement = "",x=dates)
    funs = c(ymd, ydm, mdy, myd, dmy, dym,
             ymd_h, dmy_h, mdy_h, ydm_h,
             ymd_hm, dmy_hm, mdy_hm, ydm_hm,
             ymd_hms, dmy_hms, mdy_hms, ydm_hms)
    for (tfun in funs){
      param_list = list(data = temp_dates)
      param_list$tz = 'UTC'
      dates = tryCatch({do.call(tfun,param_list)},
                       warning = function(w) {})
      if (!is.null(dates)){
        break
      }
    }
    ncdf_xts = xts(x = tt,order.by = dates)
  }

  if(Sys.info()['sysname'] == "Windows" & parallel){
    ncdf_stats = parallelsugar::mclapply(1:ncol(ncdf_xts), FUN = function(x){basic_stats(ncdf_xts[,x])$stats_table})
  }else if(parallel){
    ncdf_stats = parallel::mclapply(1:ncol(ncdf_xts), FUN = function(x){basic_stats(ncdf_xts[,x])$stats_table})
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
