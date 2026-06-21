#' L-Moment Statistics
#'
#' @description Computes the L-mean, L-scale, L-skewness, L-kurtosis, and L-CV
#' for each column of an xts object using \code{lmom::samlmu}.
#'
#' @param ts An xts object containing the time series data.
#' @param ignore_zeros Logical. If \code{TRUE}, zeros are excluded before
#'   computing L-moments. Default \code{FALSE}.
#' @param zero_threshold Numeric. Values below this threshold are treated as
#'   zero. Default 0.01.
#'
#' @return A data frame with rows for L-mean, L-scale, L-skewness, L-kurtosis,
#'   and L-CV, and one column per input series.
#'
#' @examples
#' # Synthetic daily data
#' set.seed(42)
#' ts <- xts::xts(rgamma(3650, shape = 2, scale = 5),
#'           order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = 3650))
#' ts_lmoms <- lmom_stats(ts, ignore_zeros = TRUE)
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
      ts_list <- c(ts_list, list(stats::na.omit(ts[I,i])))
    }else{
      ts_list <- c(ts_list, list(stats::na.omit(ts[,i])))
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
