#' Monthly Distribution Diagnostics
#'
#' @description Computes Q-Q and P-P diagnostic plots and Cramer-von Mises /
#' Kolmogorov-Smirnov goodness-of-fit statistics for a fitted distribution,
#' applied separately to each calendar month and assembled into 12-panel grids.
#'
#' @param ts An xts object containing the time series data.
#' @param distr Character string naming the distribution (e.g. \code{'gamma3'}).
#' @param params A named list of monthly parameter vectors, as returned by
#'   \code{\link{fitlm_monthly}}$params_monthly.
#' @param ignore_zeros Logical. If \code{TRUE}, zeros are excluded. Default
#'   \code{FALSE}.
#' @param zero_threshold Numeric. Values below this threshold are treated as
#'   zero. Default 0.01.
#' @param nrow Number of rows in the plot grid. Default 3.
#' @param ncol Number of columns in the plot grid. Default 4.
#'
#' @return A list with elements \code{monthly_QQplot} (12-panel Q-Q grid),
#'   \code{monthly_PPplot} (12-panel P-P grid), and \code{GoF_monthly} (list of
#'   per-month CM/KS statistics).
#'
#' @examples
#' # Synthetic daily data
#' set.seed(42)
#' ts <- xts::xts(rgamma(3650, shape = 2, scale = 5),
#'           order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = 3650))
#' monthly_fits <- fitlm_monthly(ts, candidates = 'gamma3', ignore_zeros = TRUE)
#' monthly_fcheck <- monthly_fdiagnostics(ts, distr = 'gamma3',
#'                                        params = monthly_fits$params_monthly,
#'                                        ignore_zeros = TRUE)
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
    param_col <- params[[distr]][[month.name[i]]]
    names(param_col) <- rownames(params[[distr]])
    params_temp <- list(params = as.list(param_col))
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
