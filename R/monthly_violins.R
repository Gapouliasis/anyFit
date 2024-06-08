#' @title monthly_violins
#'
#' @description  Plots a violin plot of monthly values
#'
#' @param ts A xts object containing the time series data
#' @param palette The color palette to use. Default is 'Set3'.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @return A ggplot2 object.
#'
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' violins <- monthly_violins(data[,4], palette = 'Set3', ignore_zeros = TRUE)
#'
#' @export
#'

monthly_violins <- function(ts,palette='Set3',ignore_zeros = FALSE,
                            zero_threshold = 0.01){
  month_ID <- month(ts)
  varnames <- colnames(ts)
  ts <- coredata(ts)
  ts <- as.data.frame(ts)
  colnames(ts) <- varnames
  if (is.null(varnames)){
    varnames = rep('V', ncol(ts))
    varnames = paste0(varnames,seq(1,ncol(ts)))
  }
  ts$month <- factor(month.name[month_ID], levels = month.name)
  ts_melt <- reshape2::melt(ts, id.vars = (length(varnames) + 1))
  if (ignore_zeros == TRUE){
    ts_melt <- ts_melt[ts_melt$value > zero_threshold,]
  }

  if (length(varnames) > 1){
    violin <- ggplot(data = ts_melt, aes(x=month, y=value, fill=variable)) +
      geom_violin(trim = FALSE) + scale_fill_brewer(palette=palette) +
      xlab('Month') + ylab('Data') + labs(fill = 'Legend')
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  }else{
    violin <- ggplot(data = ts_melt, aes(x=month, y=value, fill=month)) +
      geom_violin(trim = FALSE) +
      geom_boxplot(width=0.1) + scale_fill_brewer(palette=palette) +
      xlab('Month') + ylab('Data') +
      theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  }

  return(violin)
}
