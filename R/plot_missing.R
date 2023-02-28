#' @title plot_missing
#'
#' @description Function to plot missing values position.
#' Assigns the mean value of timeseries to NA for easy identification of empties.
#' Meant for quick and dirty (i.e. Engineering) visual checks
#'
#' @param ts An xts time series object.
#'
#' @return A `ggplot` object with the missing values.
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#'
#' plot_missing(data[,seq(1,4)])
#'
#' @export
#'


plot_missing <- function(ts){
  ylabel <- colnames(ts)
  INA <- which(is.na(ts))
  ts_df <- fortify.zoo(ts)
  names(ts_df) <- c('Index','Values')
  data_mean <- mean(ts_df$Values, na.rm = TRUE)
  ts_df$Legend <- 'Measured'
  ts_df[INA,3] <- 'Missing'
  ts_df[INA,2] <- data_mean
  f <- ggplot(data = ts_df, aes(x = Index, y = Values, colour = Legend)) + geom_point() +
    scale_color_manual(values = c('black','red')) + ylab(ylabel) +
    xlab('Date') + ggtitle('Missing Values check')
  return(f)
}

