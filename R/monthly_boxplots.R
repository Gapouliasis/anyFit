#' @title monthly_boxplots
#' 
#' @description Plot monthly boxplots for a given time series.
#'
#' @param ts A xts object containing the time series data.
#' @param palette A color palette. Default is 'Set1'.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @return A ggplot object with monthly boxplots
#' 
#' @example 
#' 
#'file <- "KNMI_Daily.csv"
#'file_path <- file.path(cwd,file)
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path, 
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'                  
#'boxes <- monthly_boxplots(data[,4], palette = 'Set3', ignore_zeros = TRUE)                  
#'
#' @export

monthly_boxplots <- function(ts, palette = 'Set1', ignore_zeros = FALSE,
                             zero_threshold = 0.01){

  month_ID <- month(ts)
  varnames <- colnames(ts)
  ts <- coredata(ts)
  ts <- as.data.frame(ts)
  colnames(ts) <- varnames
  ts$month <- factor(month.name[month_ID], levels = month.name)
  ts_melt <- reshape2::melt(ts, id.vars = (length(varnames) + 1))
  if (ignore_zeros == TRUE){
    ts_melt <- ts_melt[ts_melt$value > zero_threshold,]
  }
  
  if (length(varnames) > 1){
    boxes <- ggplot(data = ts_melt, aes(x=month, y=value, fill=variable)) + 
      geom_boxplot(width=0.8) + scale_fill_brewer(palette= palette) +
      xlab('Month') + ylab('Data') + labs(fill = 'Legend') +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  }else{
    boxes <- ggplot(data = ts_melt, aes(x=month, y=value, fill=month)) + 
      geom_boxplot(width=0.8) + scale_fill_brewer(palette= palette) +
      xlab('Month') + ylab('Data') + 
      theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  }
  
  return(boxes)
}
