#' @title ridge_xts
#' 
#' @param data A vector or matrix of data
#' @param palette Color palette to use for ridge plot. Default is 'Set3'
#' @param ignore_zeros Boolean specifying whether to ignore values below zero_threshold
#' @param zero_threshold Threshold value for ignoring data. Default is 0.01
#' 
#' @return A list with two elements: plot_all and plot_monthly. plot_all is the full ridge plot for the data, plot_monthly is a list of ridge plots for each month of the data.
#' 
#' @examples
#' data(airquality)
#' ridge_plots(airquality$Ozone)
#' 
#' @export

ridge_plots <- function(data,palette='Set3', ignore_zeros = FALSE,
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
