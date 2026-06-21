#' @title monthly_violins
#'
#' @description Produces a ggplot2 violin plot of monthly values from an xts
#'   time series. Single-column series overlay a boxplot inside each violin;
#'   multi-column series facet by variable using a brewer fill palette.
#'
#' @param ts An xts object containing the time series data. Multi-column xts
#'   are supported; each column is treated as a separate variable.
#' @param palette Character; the RColorBrewer palette name for the fill scale.
#'   Default \code{"Set3"}.
#' @param ignore_zeros Logical; if \code{TRUE}, values at or below
#'   \code{zero_threshold} are excluded. Default \code{FALSE}.
#' @param zero_threshold Numeric; threshold below which values are treated as
#'   zero. Default \code{0.01}.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' # Synthetic daily precipitation
#' set.seed(123)
#' n <- 365 * 2
#' dates <- seq(as.POSIXct("2000-01-01", tz = "UTC"), by = "day", length.out = n)
#' precip <- pmax(0, rnorm(n, mean = 3, sd = 5))
#' ts <- xts::xts(precip, order.by = dates)
#'
#' monthly_violins(ts, palette = "Set3")
#'
#' # Exclude zeros
#' monthly_violins(ts, ignore_zeros = TRUE, zero_threshold = 0.1)
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
      ggplot2::geom_violin(trim = FALSE) + ggplot2::scale_fill_brewer(palette=palette) +
      xlab('Month') + ggplot2::ylab('Data') + labs(fill = 'Legend')
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  }else{
    violin <- ggplot(data = ts_melt, aes(x=month, y=value, fill=month)) +
      ggplot2::geom_violin(trim = FALSE) +
      ggplot2::geom_boxplot(width=0.1) + ggplot2::scale_fill_brewer(palette=palette) +
      xlab('Month') + ggplot2::ylab('Data') +
      theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  }

  return(violin)
}
