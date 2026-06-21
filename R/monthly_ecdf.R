#' Monthly Empirical CDF Grid
#'
#' @description Generates a 12-panel grid of empirical cumulative distribution
#' function plots, one per calendar month. Supports multi-column xts input with
#' series distinguished by colour.
#'
#' @param ts An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, zeros are excluded. Default
#'   \code{FALSE}.
#' @param zero_threshold Numeric. Values below this threshold are treated as
#'   zero. Default 0.01.
#'
#' @return A \pkg{patchwork} grid of ggplot ECDF panels (4 rows by 3 columns).
#'
#' @examples
#' # Synthetic daily data
#' set.seed(42)
#' ts <- xts::xts(rgamma(3650, shape = 2, scale = 5),
#'           order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = 3650))
#' monthly_ecdf(ts)
#'
#' @importFrom ggplot2 ggplot aes geom_point ggtitle scale_color_brewer labs theme element_text
#' @importFrom lubridate month
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom reshape2 melt
#' @importFrom zoo coredata
#' @importFrom matrixStats colRanks
#'
#' @export
monthly_ecdf = function(ts, ignore_zeros = FALSE, zero_threshold = 0.01){
  if (ignore_zeros == TRUE){
    temp <- ts[ts>zero_threshold,]
  }else{
    temp <- ts
  }

  i_months <- unique(month(ts))
  month_name <- rep(0,length(i_months))
  monthly_plots = list()
  for (i in i_months){
    month_name[i] <- month.name[i]
    I <- which(month(ts) == i)
    monthly_ts <- ts[I]
    monthly_ts_long = reshape2::melt(coredata(monthly_ts))

    monthly_ecdf = matrixStats::colRanks(coredata(monthly_ts), ties.method = "average",
                                         preserveShape = TRUE)/nrow(monthly_ts)
    colnames(monthly_ecdf) = colnames(coredata(monthly_ts))
    monthly_ecdf_long = reshape2::melt(monthly_ecdf)
    monthly_ecdf_long$q = monthly_ts_long$value

    monthly_ecdf_plot = ggplot(monthly_ecdf_long, aes(x = q, y = value, color = factor(Var2))) +
      geom_point(shape = 1, size = 1.5, stroke = 1.5) + ggtitle(month.name[i]) +
      scale_color_brewer(palette='Set1') +
      labs(x = 'Empirical Quantile', y = 'Empirical Probability') +
    theme(legend.position = 'bottom')
    monthly_plots = c(monthly_plots, list(monthly_ecdf_plot))
  }

  plot_ecdfs = patchwork::wrap_plots(plotlist = monthly_plots,
                    nrow = 4, ncol = 3) +
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  return(plot_ecdfs)
}
