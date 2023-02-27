#' @title fitlm_monthly
#' 
#' @description This function fits a list of candidate distributions on a monthly basis using the L-moments method 
#' to a timeseries in xts format. Additionally to the list of the fitted parameters,
#'  goodness-of-fit metric, PP and QQ plots.
#' 
#' @param ts A xts object containing the time series data. 
#' @param candidates A list of distribution to fit.
#' @param nrow Number of rows for plotting.
#' @param ncol Number of columns for plotting.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' 
#' @return A list containing the monthly fitted parameters, Goodness-of-Fit Summary, and diagnostic plots.
#' 
#' @examples
#'file <- "KNMI_Daily.csv"
#'file_path <- file.path(cwd,file)
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path, 
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'                  
#'candidates <- list('exp','expweibull', 'gamma3')
#'monthly_fits <- fitlm_monthly(data[,4],candidates = candidates, ignore_zeros = TRUE)            
#' 
#' monthly_fits$monthly_QQplot
#' 
#' @export
#' 

#Function to fit multiple distributions to raw scale at a monthly level
#With considerations for missing months
fitlm_monthly <- function(ts,candidates,ignore_zeros = FALSE, zero_threshold = 0.01,
                          nrow = 4, ncol = 3){
  i_months <- unique(month(ts))
  month_name <- rep(0,length(i_months))
  for (i in i_months){
    month_name[i] <- month.name[i]
    I <- which(month(ts) == i)
    monthly_ts <- ts[I]
    monthly_fit <- fitlm_multi(monthly_ts, candidates = candidates,
                               ignore_zeros = ignore_zeros, zero_threshold = zero_threshold)
    
    monthly_fit$QQplot <- monthly_fit$QQplot + ggtitle(month.name[i]) + 
      theme(legend.position = 'bottom')
    assign(sprintf('plotQQ_%i',i),monthly_fit$QQplot)
    
    monthly_fit$PPplot <- monthly_fit$PPplot + ggtitle(month.name[i]) + 
      theme(legend.position = 'bottom')
    assign(sprintf('plotPP_%i',i),monthly_fit$PPplot)
    
    if (i == 1){
      params_monthly = list()
      GoF_monthly = list()
      for (i in 1:length(candidates)){
        params_monthly <- c(params_monthly, list(as.data.frame(unlist(monthly_fit$parameter_list[[i]]$Param))))
        GoF_monthly <- c(GoF_monthly, list(as.data.frame(monthly_fit$parameter_list[[i]]$GoF)))
      }
    }else{
      for (i in 1:length(candidates)){
        temp <- params_monthly[[i]]
        temp <- cbind(temp, as.data.frame(unlist(monthly_fit$parameter_list[[i]]$Param)))
        params_monthly[[i]] <- temp
        temp <- GoF_monthly[[i]]
        temp <- cbind(temp, as.data.frame(monthly_fit$parameter_list[[i]]$GoF))
        GoF_monthly[[i]] <- temp
      }
    }
  }
  
  names(params_monthly) <- candidates
  names(GoF_monthly) <- candidates
  for (i in 1:length(candidates)){
    temp <- params_monthly[[i]]
    colnames(temp) <- month.name
    params_monthly[[i]] <- temp
    temp <- GoF_monthly[[i]]
    colnames(temp) <- month.name
    GoF_monthly[[i]] <- temp
  }
  
  if (length(i_months)<12){
    i_complete <- seq(1,12)
    Iexist <- which(i_complete == i_months)
    i_empty <- i_complete[-Iexist,]
    for (i in i_empty){
      assign(sprintf('plotQQ_%i',i),ggplot())
      assign(sprintf('plotPP_%i',i),ggplot())
    }
  }
  
  monthly_QQplot <- ggpubr::ggarrange(plotQQ_1, plotQQ_2, plotQQ_3, plotQQ_4,
                                        plotQQ_5, plotQQ_6, plotQQ_7, plotQQ_8,
                                        plotQQ_9, plotQQ_10, plotQQ_11, plotQQ_12,
                                      nrow = nrow, ncol = ncol,
                                      common.legend = TRUE, legend = 'bottom')
  
  monthly_PPplot <- ggpubr::ggarrange(plotPP_1, plotPP_2, plotPP_3, plotPP_4,
                                      plotPP_5, plotPP_6, plotPP_7, plotPP_8,
                                      plotPP_9, plotPP_10, plotPP_11, plotPP_12,
                                      nrow = nrow, ncol = ncol,
                                      common.legend = TRUE, legend = 'bottom')
  
  out <- list('params_monthly' = params_monthly,'GoF_monthly' = GoF_monthly,
              'monthly_QQplot' = monthly_QQplot, 'monthly_PPplot' = monthly_PPplot)
}
