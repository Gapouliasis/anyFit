#' @title fitlm_nc
#'
#' This function fits a list of candidate distributions using the L-moments method NetCDF raster file.
#' It returns fitted distribution parameters, the theoretical L-moments of the fitted distribution,
#' the sample L-Moments and goodness-of-fit statistics in raster format.
#'
#' @param raster_file A raster file
#' @param filename (optional) A NetCDF file name to import if raster_file is not provided.
#' @param varname (optional) The name of the variable to extract from 'filename'.
#' @param candidates A list of distribution to fit.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param parallel Logical, whether to use parallel processing.
#' @param ... Additional arguments to pass to 'nc2xts' function (if 'filename' and 'varname' are provided).
#'
#' @return A list of raster objects containing the fitted distribution parameters, the theoretical L-moments of the fitted distribution,
#' the sample L-Moments and goodness-of-fit statistics in raster format.
#'
#' @examples
#'
#' @import raster
#' @import xts
#' @importFrom lubridate ymd ydm mdy myd dmy dym ymd_h dmy_h mdy_h ydm_h ymd_hm dmy_hm mdy_hm ydm_hm ymd_hms dmy_hms mdy_hms ydm_hms
#' @importFrom parallel mclapply
#'
#' @export
#'
fitlm_nc = function(raster_file, filename = NA, varname = NA,
                    candidates = 'norm',ignore_zeros = FALSE, zero_threshold = 0.01 ,
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
    ncdf_fits = parallelsugar::mclapply(1:ncol(ncdf_xts),
                       FUN = function(x){fitlm_nxts(ncdf_xts[,x],ignore_zeros = ignore_zeros,
                                                    candidates = candidates, zero_threshold = zero_threshold)$params[[1]]})
  }else if(parallel){
    ncdf_fits = parallel::mclapply(1:ncol(ncdf_xts),
                        FUN = function(x){fitlm_nxts(ncdf_xts[,x],ignore_zeros = ignore_zeros,
                                                     candidates = candidates, zero_threshold = zero_threshold)$params[[1]]})
  }else{
    ncdf_fits = lapply(1:ncol(ncdf_xts),
                       FUN = function(x){fitlm_nxts(ncdf_xts[,x],ignore_zeros = ignore_zeros,
                                                    candidates = candidates, zero_threshold = zero_threshold)$params[[1]]})
  }

  ncdf_params = lapply(ncdf_fits, function(x){unlist(x[[1]][[1]]$Param)})
  ncdf_params = t(do.call(cbind,ncdf_params))
  ncdf_params = cbind(coords,ncdf_params)
  #rownames(ncdf_fits) = as.character(dates)
  raster_params = raster::rasterFromXYZ(ncdf_params)
  projection(raster_params) = projection(raster_file)

  ncdf_TheorLMom = lapply(ncdf_fits, function(x){unlist(x[[1]][[1]]$TheorLMom)})
  ncdf_TheorLMom = t(do.call(cbind,ncdf_TheorLMom))
  ncdf_TheorLMom = cbind(coords,ncdf_TheorLMom)
  #rownames(ncdf_fits) = as.character(dates)
  raster_TheorLMom = raster::rasterFromXYZ(ncdf_TheorLMom)
  projection(raster_TheorLMom) = projection(raster_TheorLMom)

  ncdf_DataLMom = lapply(ncdf_fits, function(x){unlist(x[[1]][[1]]$DataLMom)})
  ncdf_DataLMom = t(do.call(cbind,ncdf_DataLMom))
  ncdf_DataLMom = cbind(coords,ncdf_DataLMom)
  #rownames(ncdf_fits) = as.character(dates)
  raster_DataLMom = raster::rasterFromXYZ(ncdf_DataLMom)
  projection(raster_DataLMom) = projection(raster_DataLMom)

  ncdf_GoF = lapply(ncdf_fits, function(x){unlist(x[[1]][[1]]$GoF)})
  ncdf_GoF = t(do.call(cbind,ncdf_GoF))
  ncdf_GoF = cbind(coords,ncdf_GoF)
  #rownames(ncdf_fits) = as.character(dates)
  raster_GoF = raster::rasterFromXYZ(ncdf_GoF)
  projection(raster_GoF) = projection(raster_GoF)

  list_out = list(raster_params = raster_params, raster_TheorLMom = raster_TheorLMom,
                  raster_DataLMom = raster_DataLMom, raster_GoF = raster_GoF)

  return(list_out)
}


