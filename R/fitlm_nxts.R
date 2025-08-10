#' @title fitlm_nxts
#'
#' @description This function fits a list of candidate distributions using the L-moments method
#' to an arbitrary number of timeseries in xts format. Additionally to the list of the fitted parameters,
#'  goodness-of-fit metric, PP and QQ plots.
#'
#' @param ts A xts object containing the time series data.
#' @param candidates A list of distribution to fit.
#' @param nrow Number of rows for plotting.
#' @param ncol Number of columns for plotting.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param parallel Logical, whether to use parallel processing.
#' @param ncores Number of cores to use in the case of parallel computations
#'
#' @return A list with the estimated parameters, diagnostic plots, QQ plots and PP plots.
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' candidates <- list('exp','expweibull', 'gamma3')
#'
#' fits_all <- fitlm_nxts(data, candidates, nrow = 5, ncol = 4, ignore_zeros = TRUE)
#' fits_all$QQ_panels[[1]]
#'
#' @export

fitlm_nxts <- function(ts, candidates, nrow = 5, ncol = 4, ignore_zeros = FALSE,
                       zero_threshold = 0.01, parallel = FALSE, ncores = 2){
  variables <- colnames(ts)
  ts_list <- list()
  for (i in 1:ncol(ts)){
    ts_list <- c(ts_list, list(na.omit(ts[,i])))
  }

  if(Sys.info()['sysname'] == "Windows" & parallel){
    multi_fits <- parallelsugar::mclapply(ts_list, FUN = fitlm_multi,candidates = candidates, ignore_zeros = ignore_zeros,
                         zero_threshold = zero_threshold, mc.cores = ncores)
  }else if(parallel){
    multi_fits <- parallel::mclapply(ts_list, FUN = fitlm_multi,candidates = candidates, ignore_zeros = ignore_zeros,
                                          zero_threshold = zero_threshold, mc.cores = ncores)
  }else{
    multi_fits <- lapply(ts_list, FUN = fitlm_multi,candidates = candidates, ignore_zeros = ignore_zeros,
                         zero_threshold = zero_threshold)
  }

  params <- list()
  diagnostic_plots <- list()
  QQ_plots <- list()
  PP_plots <- list()
  for (i in 1:ncol(ts)){
    params <- c(params,list(list(params = multi_fits[[i]]$parameter_list,GoF_summary = multi_fits[[i]]$GoF_summary)))
    diagnostic_plots <- c(diagnostic_plots,list(multi_fits[[i]]$diagnostics))
    QQ_plots <- c(QQ_plots,l =  list(multi_fits[[i]]$QQplot + ggtitle(variables[i])))
    PP_plots <- c(PP_plots,l =  list(multi_fits[[i]]$PPplot + ggtitle(variables[i])))
  }

  names(params) <- variables
  names(diagnostic_plots) <- variables
  names(QQ_plots) <- variables
  names(PP_plots) <- variables

  n_plots <- nrow*ncol
  panels <- pracma::ceil(ncol(ts)/n_plots)
  QQ_panels <- list()
  PP_panels <- list()
  for (i in 1:panels){
    if ((i*(n_plots)) < ncol(ts)){
      temp_QQ <- QQ_plots[(1 + (i-1)*n_plots):(i*(n_plots))]
      temp_PP <- PP_plots[(1 + (i-1)*n_plots):(i*(n_plots))]
    }else{
      temp_QQ <- QQ_plots[(1 + (i-1)*n_plots):ncol(ts)]
      temp_PP <- PP_plots[(1 + (i-1)*n_plots):ncol(ts)]
    }
    QQ_panels <- c(QQ_panels, list(patchwork::wrap_plots(plotlist = temp_QQ, nrow = nrow, ncol = ncol)+
                                     plot_layout(guides = "collect") & theme(legend.position = 'bottom')))

    PP_panels <- c(PP_panels, list(patchwork::wrap_plots(plotlist = temp_PP, nrow = nrow, ncol = ncol)+
                                     plot_layout(guides = "collect") & theme(legend.position = 'bottom')))
  }

  return(list(params = params, diagnostic_plots = diagnostic_plots, QQ_plots = QQ_plots,
              PP_plots = PP_plots, QQ_panels = QQ_panels, PP_panels = PP_panels))
}

