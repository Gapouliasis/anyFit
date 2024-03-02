#' @title monthly_fdiagnostics
#'
#' @description Performs a goodness-of-fit test to diagnose fitted distributions by month.
#'
#' @param ts A xts object containing the time series data.
#' @param distr The distribution function to be tested.
#' @param params A list object containing the parameters of the fitted distribution.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param nrow Number of rows of the plots.
#' @param ncol Number of columns of the plots.
#'
#' @return A list containing a Q-Q plot and a P-P plot for each month and a goodness-of-fit table.
#' The GoF metric calculated are the Kramer von-Mises and Kolmogorov-Smirnov. NEEDS EXPANSION
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#'monthly_fits <- fitlm_monthly(smonthly_ts,candidates = 'gamma3', ignore_zeros = TRUE)
#'
#'monthly_fcheck <- monthly_fdiagnostics(simXts, distr = 'gamma3', params = monthly_fits$params_monthly, ignore_zeros = TRUE)
#'
#' @export
#'
#'

monthly_fdiagnostics <- function(ts, distr, params, ignore_zeros = FALSE, zero_threshold = 0.01,
                                  nrow = 3, ncol = 4){
  pfunction = paste0('p',distr)
  qfunction = paste0('q',distr)
  i_months <- unique(month(ts))
  month_name <- rep(0,length(i_months))
  for (i in i_months){
    month_name[i] <- month.name[i]
    I <- which(month(ts) == i)
    monthly_ts <- ts[I]
    params_temp <- list(param = eval(parse(text = paste0('params$',month.name[i]))))
    params_temp <- c(params_temp, list(ts = monthly_ts, dist = distr,
                                       ignore_zeros = ignore_zeros, zero_threshold = zero_threshold))
    monthly_fit <- do.call(fit_diagnostics, params_temp)

    monthly_fit$QQplot <- monthly_fit$QQplot + ggtitle(month.name[i]) +
      theme(legend.position = 'bottom')
    assign(sprintf('plotQQ_%i',i),monthly_fit$QQplot)

    monthly_fit$PPplot <- monthly_fit$PPplot + ggtitle(month.name[i]) +
      theme(legend.position = 'bottom')
    assign(sprintf('plotPP_%i',i),monthly_fit$PPplot)

    if (i == 1){
      GoF_monthly <- list(monthly_fit$GoF)
    }else{
      GoF_monthly <- c(GoF_monthly, list(monthly_fit$GoF))
    }
  }

  names(GoF_monthly) <- month_name

  if (length(i_months)<12){
    i_complete <- seq(1,12)
    Iexist <- which(i_complete == i_months)
    i_empty <- i_complete[-Iexist,]
    for (i in i_empty){
      assign(sprintf('plotQQ_%i',i),ggplot())
      assign(sprintf('plotPP_%i',i),ggplot())
    }
  }

  monthly_QQplot <- patchwork::wrap_plots(plotQQ_1, plotQQ_2, plotQQ_3, plotQQ_4,
                                      plotQQ_5, plotQQ_6, plotQQ_7, plotQQ_8,
                                      plotQQ_9, plotQQ_10, plotQQ_11, plotQQ_12,
                                      nrow = nrow, ncol = ncol)

  monthly_PPplot <- patchwork::wrap_plots(plotPP_1, plotPP_2, plotPP_3, plotPP_4,
                                      plotPP_5, plotPP_6, plotPP_7, plotPP_8,
                                      plotPP_9, plotPP_10, plotPP_11, plotPP_12,
                                      nrow = nrow, ncol = ncol)

  out <- list('monthly_QQplot' = monthly_QQplot, 'monthly_PPplot' = monthly_PPplot, 'GoF_monthly' = GoF_monthly)

}
