#' @title normalise_xts
#'
#' @description Empirical normal-score transform of an xts time series: values
#'   are mapped to standard normal quantiles via their empirical cumulative
#'   distribution function (ranking with average ties). Optionally applied
#'   per calendar month to preserve seasonal structure in cyclo-stationary
#'   processes.
#'
#' @param ts An xts object containing the time series data. Multi-column xts
#'   are supported; each column is normalised independently.
#' @param dist_period Character; the distribution period for normalisation.
#'   Use \code{"monthly"} for per-calendar-month normalisation or \code{NA}
#'   for the entire series. Default \code{"monthly"}.
#' @param ignore_zeros Logical; if \code{TRUE}, zero values (at or below
#'   \code{zero_threshold}) are excluded before normalisation. Default
#'   \code{FALSE}.
#' @param zero_threshold Numeric; threshold below which values are treated as
#'   zero. Default \code{0.01}.
#'
#' @return A named list of normalised xts objects, one per column of the input.
#'
#' @examples
#' # Synthetic daily precipitation
#' set.seed(42)
#' n <- 365 * 3
#' dates <- seq(as.POSIXct("2000-01-01", tz = "UTC"), by = "day", length.out = n)
#' precip <- pmax(0, rnorm(n, mean = 3, sd = 5))
#' ts <- xts::xts(precip, order.by = dates)
#'
#' # Monthly normal-score transform
#' ns <- normalise_xts(ts, dist_period = "monthly")
#' hist(zoo::coredata(ns[[1]]), main = "Normalised values")
#'
#' # Full-series transform
#' ns_full <- normalise_xts(ts, dist_period = NA)
#'
#' @importFrom lubridate month
#' @importFrom xts xts rbind.xts
#' @importFrom zoo coredata index
#'
#' @export

normalise_xts = function(ts, dist_period = 'monthly',
                         ignore_zeros = FALSE, zero_threshold = 0.01){

  out_list = list()
  ts_orig = ts

  for (c in 1:ncol(ts_orig)){
    ts = ts_orig[,c]
    if (ignore_zeros == TRUE){
      Izeros = which(ts<=zero_threshold)
      ts <- ts[-Izeros,]
    }

    if (isTRUE(dist_period == 'monthly')){
      temp_list = list()
      for (i in 1:12){
        I <- which(month(ts) == i)
        temp_data = coredata(ts[I])
        Ina = which(is.na(temp_data))
        u_emp = rank(temp_data, na.last = NA, ties.method = "average")/(length(temp_data)+1)
        z_emp = stats::qnorm(p = u_emp, mean = 0, sd = 1)
        temp_dates = index(ts)[I]
        if(length(Ina)>0){temp_dates = temp_dates[-Ina]}
        temp_xts = xts(x = z_emp, order.by = temp_dates)
        temp_list = c(temp_list, list(temp_xts))
      }

      # if (ignore_zeros == TRUE){
      #   temp_xts = xts(x = z_emp, order.by = temp_dates)
      #   temp_list = c(temp_list, list(temp_xts))
      # }

      normal_xts = do.call(rbind.xts, temp_list)

    }else{
      temp_data = coredata(ts)
      u_emp = rank(temp_data, na.last = NA, ties.method = "average")/(length(temp_data)+1)
      z_emp = stats::qnorm(p = u_emp, mean = 0, sd = 1)
      normal_xts = xts(x = z_emp, order.by = index(ts))
    }

    out_list = c(out_list, list(normal_xts))
  }

  names(out_list) = colnames(ts_orig)
  return(out_list)
 }
