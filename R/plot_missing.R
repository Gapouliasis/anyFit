#' Plot missing value positions
#'
#' @description
#' Replaces NAs with the series mean to make gaps visually identifiable as a
#' horizontal line at the mean level. Intended for quick diagnostic checks.
#'
#' @param ts An xts time series object.
#'
#' @return A \code{ggplot} object marking measured and missing values.
#'
#' @examples
#' # Synthetic xts with artificial gaps
#' set.seed(123)
#' dates <- seq(as.POSIXct("2000-01-01"), as.POSIXct("2000-12-31"), by = "day")
#' vals <- matrix(rnorm(length(dates) * 2, 10, 5), nrow = length(dates), ncol = 2)
#' colnames(vals) <- c("A", "B")
#' vals[sample(length(vals), 20)] <- NA
#' ts <- xts::xts(vals, order.by = dates)
#'
#' plot_missing(ts)
#'
#' @export
#'


plot_missing <- function(ts){
  ylabel <- colnames(ts)
  INA <- which(is.na(ts))
  ts_df <- zoo::fortify.zoo(ts)
  names(ts_df) <- c('Index','Values')
  data_mean <- mean(ts_df$Values, na.rm = TRUE)
  ts_df$Legend <- 'Measured'
  ts_df[INA,3] <- 'Missing'
  ts_df[INA,2] <- data_mean
  f <- ggplot(data = ts_df, aes(x = Index, y = Values, colour = Legend)) + geom_point() +
    ggplot2::scale_color_manual(values = c('black','red')) + ggplot2::ylab(ylabel) +
    xlab('Date') + ggtitle('Missing Values check')
  return(f)
}
