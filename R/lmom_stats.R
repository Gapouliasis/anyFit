#' @title lmom_stats
#'
#' @description Calculate the L-Mean, L-Scale, L-Skew, L-Kurtosis and L-CV of a time series.
#'
#' @param ts A xts object containing the time series data.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @return A dataframe with the L-Mean, L-Scale, L-Skew, L-Kurtosis and L-CV of the time series
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#' ts_lmoms <- lmom_stats(data, ignore_zeros = TRUE)
#'
#'
#' @export
#'

lmom_stats <- function(ts, ignore_zeros = FALSE, zero_threshold = 0.01){
  variables <- colnames(ts)
  ts <- coredata(ts)
  ts_list <- list()
  for (i in 1:ncol(ts)){
    if (ignore_zeros == TRUE){
      I <- which(ts[,i] > zero_threshold)
      ts_list <- c(ts_list, list(na.omit(ts[I,i])))
    }else{
      ts_list <- c(ts_list, list(na.omit(ts[,i])))
    }

  }

  lstats <- function(x){
    lmoms <- round(lmom::samlmu(x,nmom = 4, ratios = FALSE),2)
    lCV <- round(lmoms[2]/lmoms[1],4)
    lskew <- round(lmoms[3]/lmoms[2],4)
    lkurtosis <- round(lmoms[4]/lmoms[2],4)
    lmom_basic <- data.frame(Value=c(lmoms[1],lmoms[2],lskew,lkurtosis,lCV))
    row.names(lmom_basic) <- c('L-Mean','L-Scale',
                               'L-Skew','L-Kurtosis','L-CV')
    return(lmom_basic)
  }

  lmom_basic <- lapply(ts_list, lstats)
  lmom_basic <- as.data.frame(lmom_basic)
  colnames(lmom_basic) <- variables

  return(lmom_basic)

}
