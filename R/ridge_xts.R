#' Ridge density plots for time series
#'
#' @description
#' Generates ggridges density ridge plots for all series combined and
#' per-series per-calendar-month variants. Supports optional zero-value
#' exclusion via \code{ignore_zeros} and \code{zero_threshold}.
#'
#' @param ts An xts object containing the time series data.
#' @param palette Colour palette for the ridge plot (default \code{"Set3"}).
#' @param ignore_zeros Logical; if \code{TRUE}, values below \code{zero_threshold} are excluded.
#' @param zero_threshold Numeric threshold below which values are treated as zero (default 0.01).
#'
#' @return A list with two elements: \code{plot_all} (full ridge plot for all
#'   series) and \code{plot_monthly} (named list of per-month ridge plots per
#'   series).
#'
#' @examples
#' # Synthetic xts
#' set.seed(123)
#' dates <- seq(as.POSIXct("2000-01-01"), as.POSIXct("2000-12-31"), by = "day")
#' vals <- matrix(rnorm(length(dates) * 4, 10, 5), nrow = length(dates), ncol = 4)
#' colnames(vals) <- c("S1", "S2", "S3", "S4")
#' ts <- xts::xts(vals, order.by = dates)
#'
#' ridges <- ridge_plots(ts, ignore_zeros = TRUE)
#' ridges$plot_all
#'
#' @export

ridge_plots <- function(ts,palette='Set3', ignore_zeros = FALSE,
                        zero_threshold = 0.01){
  ts = coredata(data.frame(Y=as.matrix(ts), date=time(ts)))
  ts = reshape2::melt(ts, id.vars= 'date')
  if (ignore_zeros == TRUE){
    ts <- ts[ts$value > zero_threshold,]
  }

  # Plot all the stations
  ridge_plot_all <- ggplot(ts, aes(x = `value`, y = `variable`, fill = ..x..)) +
    ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)+
    viridis::scale_fill_viridis()+ggplot2::ylab('')+xlab('')+theme(legend.position = 'none')
  #scale_fill_brewer(palette=palette)
  #viridis::scale_fill_viridis()
  #Plot monthly ridge plots for every station
  ridge_plot_month = list()
  for (station in unique(ts$variable)){
    dt = ts[ts$variable==station,]

    dt$month = month(dt$date)
    dt$month = month.abb[dt$month]

    ridge_plot_month[[station]]<- ggplot(dt, aes(x = `value`, y = `month`, fill = ..x..)) +
      ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)+
      viridis::scale_fill_viridis()+ggplot2::ylab('Month')+xlab('')+ggtitle(station) + theme(legend.position = 'none')

  }

  return(list(plot_all=ridge_plot_all, plot_monthly = ridge_plot_month))

}
