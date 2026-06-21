#' @title fit_ACF_monthly
#'
#' @description Wraps \code{\link{fit_ACF}} to fit ACF models separately
#'   to each calendar-month sub-series of an xts object. The input series
#'   is split by calendar month, \code{fit_ACF} is applied to each
#'   resulting sub-series, and the fitted parameters are assembled into a
#'   parameter matrix with months as rows. A 12-panel patchwork plot of
#'   fitted ACF curves is produced, with empty panels for months absent
#'   from the data. The grid layout is controlled by \code{nrow} and
#'   \code{ncol}. The three model types (CAS, HK, SRD) supported by
#'   \code{fit_ACF} are available and may be specified singly or in
#'   combination.
#'
#' @param ts An xts object containing the time series data.
#' @param lag Integer. Maximum lag passed to \code{fit_ACF}.
#' @param type Character vector of ACF model types. Options are
#'   \code{"CAS"}, \code{"HK"}, and \code{"SRD"}. Default all three.
#' @param nrow Integer. Number of rows in the panel grid. Default \code{4}.
#' @param ncol Integer. Number of columns in the panel grid. Default
#'   \code{3}.
#'
#' @return A list with elements \code{ACF_params_monthly} (matrix of
#'   fitted parameters with months as rows) and
#'   \code{ACF_monthly_plot} (12-panel patchwork of fitted ACF
#'   \code{ggplot} objects with a shared legend).
#'
#' @examples
#' ts <- xts::xts(rnorm(365), order.by = as.Date("2020-01-01") + 0:364)
#' monthly_acfs <- fit_ACF_monthly(ts, lag = 10)
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
    monthly_acf <- fit_ACF(ts = monthly_ts, lag_max = lag, type = type)
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

  ACF_monthly_plot <-  patchwork::wrap_plots(plot_1, plot_2, plot_3, plot_4,
                                        plot_5, plot_6, plot_7, plot_8,
                                        plot_9, plot_10, plot_11, plot_12, nrow = nrow, ncol = ncol) +
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')

  out <- list("ACF_params_monthly" = ACF_params_monthly,
              "ACF_monthly_plot" = ACF_monthly_plot)
  return(out)
}
