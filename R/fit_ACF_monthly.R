#' @title fit_ACF_monthly
#'
#' @description Function to fit theoretical Autocorrelation Functions (ACFs) to monthly values of timeseries data.
#' Three functions are available, the Cauchy-type autocorrelation structure,
#' the Hurst - Kolmogorov and the Short Range dependence structure.
#'
#' @param ts A xts object containing the time series data.
#' @param lag_max Maximum lag to use in fitting.
#' @param type A list of character strings, containing the type of ACF to fit.
#'   Options include: \code{'CAS'}, \code{'HK'}, and \code{'SRD'}.
#' @param nrow Number of rows for plotting.
#' @param ncol Number of columns for plotting.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @return A list containing a vector of monthly fitted parameters, a data frame of
#'   fitted ACF values, and a plot of the fitted ACFs.
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' monthly_acfs <- fit_ACF_monthly(data[,4], lag = 10)
#' monthly_acfs$ACF_monthly_plot
#'
#' @export
#'

fit_ACF_monthly <- function(ts,lag,type = list('CAS','HK','SRD'), nrow = 4, ncol = 3){
  i_months <- unique(month(ts))
  month_name <- rep(0,length(i_months))
  for (i in i_months){
    month_name[i] <- month.name[i]
    I <- which(month(ts) == i)
    monthly_ts <- ts[I]
    monthly_acf <- fit_ACF(ts = monthly_ts, lag = lag, type = type)
    monthly_acf$ACF_plot <- monthly_acf$ACF_plot + ggtitle(month.name[i]) +
      theme(legend.position = "bottom")
    assign(sprintf('plot_%i',i),monthly_acf$ACF_plot)
    if (i == 1){
      ACF_params_monthly <- unlist( monthly_acf$ACF_params)
    }else{
      ACF_params_monthly <- rbind(ACF_params_monthly, unlist(monthly_acf$ACF_params))
    }
  }

  rownames(ACF_params_monthly) <- month_name

  if (length(i_months)<12){
    i_complete <- seq(1,12)
    Iexist <- which(i_complete == i_months)
    i_empty <- i_complete[-Iexist,]
    for (i in i_empty){
      assign(sprintf('plot_%i',i),ggplot())
    }
  }

  ACF_monthly_plot <- ggpubr::ggarrange(plot_1, plot_2, plot_3, plot_4,
                                        plot_5, plot_6, plot_7, plot_8,
                                        plot_9, plot_10, plot_11, plot_12, nrow = nrow, ncol = ncol,
                                        common.legend = TRUE, legend = "bottom")

  out <- list("ACF_params_monthly" = ACF_params_monthly,
              "ACF_monthly_plot" = ACF_monthly_plot)
  return(out)
}
