#' @title ridge_xts
#'
#' @description Function that generates ridge plots from timeseries.
#' It will also generate ridge plots on the monthly scale.
#'
#' @param ts A xts object containing the time series data.
#' @param palette Color palette to use for ridge plot. Default is 'Set3'
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @return A list with two elements: plot_all and plot_monthly. plot_all is the full ridge plot for the data, plot_monthly is a list of ridge plots for each month of the data.
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#'ridges <- ridge_plots(data[,seq(1,10)], ignore_zeros = TRUE)
#'ridges$plot_all
#'ridges$plot_monthly$Y.Xtemp.2
#'
#' @export

ridge_plots <- function(ts,palette='Set3', ignore_zeros = FALSE,
                        zero_threshold = 0.01){
  ts = coredata(data.frame(Y=as.matrix(data), date=time(data)))
  ts = reshape2::melt(ts, id.vars= 'date')
  if (ignore_zeros == TRUE){
    ts <- ts[ts$value > zero_threshold,]
  }

  # Plot all the stations
  ridge_plot_all <- ggplot(ts, aes(x = `value`, y = `variable`, fill = ..x..)) +
    ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01)+
    viridis::scale_fill_viridis()+ylab('')+xlab('')+theme(legend.position = 'none')
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
      viridis::scale_fill_viridis()+ylab('Month')+xlab('')+ggtitle(station) + theme(legend.position = 'none')

  }

  return(list(plot_all=ridge_plot_all, plot_monthly = ridge_plot_month))

}
