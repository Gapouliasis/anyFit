#' @title fitlm_nc
#'
#' @description This function fits a list of candidate distributions using the L-moments method NetCDF raster file.
#' It returns fitted distribution parameters, the theoretical L-moments of the fitted distribution,
#' the sample L-Moments and goodness-of-fit statistics in raster format.
#'
#' @param raster_file A raster file.
#' @param filename (optional) A NetCDF file name to import if raster_file is not provided.
#' @param varname (optional) The name of the variable to extract from 'filename'.
#' @param candidates A list of distribution to fit.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param parallel Logical, whether to use parallel processing.
#' @param ncores Number of cores to use in the case of parallel computations
#' @param ... Additional arguments to pass to 'nc2xts' function (if 'filename' and 'varname' are provided).
#'
#' @return A list of raster objects containing the fitted distribution parameters, the theoretical L-moments of the fitted distribution,
#' the sample L-Moments and goodness-of-fit statistics in raster format.
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
fitlm_nc = function(raster_file, filename = NA, varname = NA,
                    candidates = 'norm',ignore_zeros = FALSE, zero_threshold = 0.01 ,
                    parallel = FALSE, ncores = 2, ...){

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
    funs = c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "ymd H", "dmy H", "mdy H",
             "ydm H", "ymd HM", "dmy HM", "mdy HM", "ydm HM", "ymd HMS", "dmy HMS",
             "mdy HMS", "ydm HMS")

    dates=parse_date_time(temp_dates, orders = funs)
    ncdf_xts = xts(x = tt,order.by = dates)
  }

  if(Sys.info()['sysname'] == "Windows" & parallel){
    ncdf_fits = parallelsugar::mclapply(1:ncol(ncdf_xts),
                                        FUN = function(x){fitlm_nxts(ncdf_xts[,x],ignore_zeros = ignore_zeros,
                                                                     candidates = candidates, zero_threshold = zero_threshold)$params[[1]]}, mc.cores = ncores)
  }else if(parallel){
    ncdf_fits = parallel::mclapply(1:ncol(ncdf_xts),
                                   FUN = function(x){fitlm_nxts(ncdf_xts[,x],ignore_zeros = ignore_zeros,
                                                                candidates = candidates, zero_threshold = zero_threshold)$params[[1]]}, mc.cores = ncores)
  }else{
    ncdf_fits = lapply(1:ncol(ncdf_xts),
                       FUN = function(x){fitlm_nxts(ncdf_xts[,x],ignore_zeros = ignore_zeros,
                                                    candidates = candidates, zero_threshold = zero_threshold)$params[[1]]})
  }

  result_list <- list()
  gof_list <- list()

  for (candidate in candidates){
    ncdf_params = lapply(ncdf_fits, function(x) {
      unlist(x[[1]][[candidate]]$Param)
    })
    ncdf_params = t(do.call(cbind, ncdf_params))
    ncdf_params = cbind(coords, ncdf_params)
    raster_params = raster::rasterFromXYZ(ncdf_params)
    projection(raster_params) = projection(raster_file)
    ncdf_TheorLMom = lapply(ncdf_fits, function(x) {
      unlist(x[[1]][[candidate]]$TheorLMom)
    })
    ncdf_TheorLMom = t(do.call(cbind, ncdf_TheorLMom))
    ncdf_TheorLMom = cbind(coords, ncdf_TheorLMom)
    raster_TheorLMom = raster::rasterFromXYZ(ncdf_TheorLMom)
    projection(raster_TheorLMom) = projection(raster_file)
    ncdf_DataLMom = lapply(ncdf_fits, function(x) {
      unlist(x[[1]][[candidate]]$DataLMom)
    })
    ncdf_DataLMom = t(do.call(cbind, ncdf_DataLMom))
    ncdf_DataLMom = cbind(coords, ncdf_DataLMom)
    raster_DataLMom = raster::rasterFromXYZ(ncdf_DataLMom)
    projection(raster_DataLMom) = projection(raster_file)
    ncdf_GoF = lapply(ncdf_fits, function(x) {
      unlist(x[[1]][[candidate]]$GoF)
    })
    ncdf_GoF = t(do.call(cbind, ncdf_GoF))
    ncdf_GoF = cbind(coords, ncdf_GoF)
    raster_GoF = raster::rasterFromXYZ(ncdf_GoF)
    projection(raster_GoF) = projection(raster_file)
    ncdf_GoF <- as.data.frame(ncdf_GoF)
    ncdf_GoF$Distribution = as.factor(candidate)
    gof_list[[candidate]] = ncdf_GoF
    result_list[[candidate]] = list(raster_params = raster_params, raster_TheorLMom = raster_TheorLMom,
                                    raster_DataLMom = raster_DataLMom, raster_GoF = raster_GoF)
  }

  gof_df = do.call(rbind, gof_list)
  gof_df = gof_df[, c("MLE", "CM" ,"KS" ,"MSEquant", "DiffOfMax", "MeanDiffOf10Max", "Distribution")]
  gof_long = reshape2::melt(gof_df, id.vars = "Distribution")
  gof_plot = ggplot(data = gof_long, aes(x=value)) + geom_density(alpha=0.3, aes(fill = Distribution)) +
    theme_bw() + facet_wrap(vars(variable),  scales = "free")

  list_out <- list(fit_results = result_list, gof_plots = gof_plot)
  return(list_out)
}


