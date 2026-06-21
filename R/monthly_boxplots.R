#' Monthly Boxplots
#'
#' @description Produces a ggplot boxplot panel grouped by calendar month.
#' Multi-column xts input is faceted by variable; single-column input is
#' coloured by month.
#'
#' @param ts An xts object containing the time series data.
#' @param palette An RColorBrewer palette name. Default \code{'Set1'}.
#' @param ignore_zeros Logical. If \code{TRUE}, zeros are excluded. Default
#'   \code{FALSE}.
#' @param zero_threshold Numeric. Values below this threshold are treated as
#'   zero. Default 0.01.
#'
#' @return A ggplot object with monthly boxplots.
#'
#' @examples
#' # Synthetic daily data
#' set.seed(42)
#' ts <- xts::xts(rgamma(3650, shape = 2, scale = 5),
#'           order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = 3650))
#' boxes <- monthly_boxplots(ts, palette = 'Set3', ignore_zeros = TRUE)
#'
#' @export

monthly_boxplots <- function(ts, palette = 'Set1', ignore_zeros = FALSE,
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
    boxes <- ggplot(data = ts_melt, aes(x=month, y=value, fill=variable)) +
      ggplot2::geom_boxplot(width=0.8) + ggplot2::scale_fill_brewer(palette= palette) +
      xlab('Month') + ggplot2::ylab('Data') + labs(fill = 'Legend') +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  }else{
    boxes <- ggplot(data = ts_melt, aes(x=month, y=value, fill=month)) +
      ggplot2::geom_boxplot(width=0.8) + ggplot2::scale_fill_brewer(palette= palette) +
      xlab('Month') + ggplot2::ylab('Data') +
      theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  }

  return(boxes)
}
