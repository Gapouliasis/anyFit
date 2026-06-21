#' @title Monthly Distribution Fitting
#'
#' @description Fits a list of candidate distributions to each calendar month of
#'   a time series via the L-moments method. The series is partitioned by month,
#'   and \code{\link{fitlm_multi}} is applied to each monthly subset independently.
#'   Parameter estimates for each candidate are aggregated into monthly tables,
#'   and the six goodness-of-fit metrics (MLE, CM, KS, MSEquant, DiffOfMax,
#'   MeanDiffOf10Max) are collected into per-candidate matrices with months as
#'   rows. Twelve-panel Q-Q and P-P diagnostic grids are assembled with
#'   \code{patchwork::wrap_plots}, using the \code{nrow} and \code{ncol} arguments
#'   to control the layout. Months absent from the data are filled with empty
#'   \code{ggplot} panels so that the grid always has twelve cells. The returned
#'   list provides the monthly parameter tables, GoF matrices, and the combined
#'   Q-Q and P-P panel plots.
#'
#' @param ts An xts object containing the time series data.
#' @param candidates A character vector of distribution names to fit.
#' @param nrow Number of rows for the diagnostic plot grid. Default is 4.
#' @param ncol Number of columns for the diagnostic plot grid. Default is 3.
#' @param ignore_zeros A logical value, if \code{TRUE} zeros will be ignored.
#'   Default is \code{FALSE}.
#' @param zero_threshold The threshold below which values are considered zero.
#'   Default is 0.01.
#' @param order Optional named list mapping a candidate name to the vector of
#'   L-moment orders matched by its optimiser, e.g.
#'   \code{list(gengamma = 1:5, expweibull = 1:3)}. Only the numerically-fitted
#'   distributions accept it; passed through to \code{\link{fitlm_multi}}.
#'   Default \code{NULL}.
#'
#' @return A list with components \code{params_monthly} (a list of data frames,
#'   one per candidate, with columns for each month), \code{GoF_monthly} (a list
#'   of transposed matrices, one per candidate, with months as rows and GoF
#'   metrics as columns), \code{monthly_QQplot}, and \code{monthly_PPplot}
#'   (twelve-panel diagnostic grids).
#'
#' @examples
#' # Two years of daily precipitation-like data
#' x <- xts::xts(rgamma(730, shape = 0.8, scale = 3),
#'          order.by = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 730))
#' x[sample(1:730, 200)] <- 0
#'
#' candidates <- c('exp', 'gamma3')
#' monthly_fits <- fitlm_monthly(x, candidates = candidates, ignore_zeros = TRUE)
#' monthly_fits$monthly_QQplot
#'
#' @export
#'

#Function to fit multiple distributions to raw scale at a monthly level
#With considerations for missing months
fitlm_monthly <- function(ts,candidates,ignore_zeros = FALSE, zero_threshold = 0.01,
                          nrow = 4, ncol = 3, order = NULL){
  i_months <- unique(month(ts))
  month_name <- rep(0,length(i_months))
  for (i in i_months){
    month_name[i] <- month.name[i]
    I <- which(month(ts) == i)
    monthly_ts <- ts[I]
    monthly_fit <- fitlm_multi(monthly_ts, candidates = candidates,
                               ignore_zeros = ignore_zeros, zero_threshold = zero_threshold,
                               order = order)

    monthly_fit$QQplot <- monthly_fit$QQplot + ggtitle(month.name[i]) +
      theme(legend.position = 'bottom')
    assign(sprintf('plotQQ_%i',i),monthly_fit$QQplot)

    monthly_fit$PPplot <- monthly_fit$PPplot + ggtitle(month.name[i]) +
      theme(legend.position = 'bottom')
    assign(sprintf('plotPP_%i',i),monthly_fit$PPplot)

    if (i == 1){
      params_monthly = list()
      GoF_monthly = list()
      for (j in 1:length(candidates)){
        params_monthly <- c(params_monthly, list(as.data.frame(unlist(monthly_fit$parameter_list[[j]]$Param))))
        GoF_monthly <- c(GoF_monthly, list(as.data.frame(monthly_fit$parameter_list[[j]]$GoF)))
      }
    }else{
      for (j in 1:length(candidates)){
        temp <- params_monthly[[j]]
        temp <- cbind(temp, as.data.frame(unlist(monthly_fit$parameter_list[[j]]$Param)))
        params_monthly[[j]] <- temp
        temp <- GoF_monthly[[j]]
        temp <- rbind(temp, as.data.frame(monthly_fit$parameter_list[[j]]$GoF))
        GoF_monthly[[j]] <- temp
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
    rownames(temp) <- month.name
    GoF_monthly[[i]] <- t(temp)
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

  monthly_QQplot <-  patchwork::wrap_plots(plotQQ_1, plotQQ_2, plotQQ_3, plotQQ_4,
                                        plotQQ_5, plotQQ_6, plotQQ_7, plotQQ_8,
                                        plotQQ_9, plotQQ_10, plotQQ_11, plotQQ_12,
                                      nrow = nrow, ncol = ncol) +
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')

  monthly_PPplot <-  patchwork::wrap_plots(plotPP_1, plotPP_2, plotPP_3, plotPP_4,
                                      plotPP_5, plotPP_6, plotPP_7, plotPP_8,
                                      plotPP_9, plotPP_10, plotPP_11, plotPP_12,
                                      nrow = nrow, ncol = ncol) +
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')

  out <- list('params_monthly' = params_monthly,'GoF_monthly' = GoF_monthly,
              'monthly_QQplot' = monthly_QQplot, 'monthly_PPplot' = monthly_PPplot)
}
