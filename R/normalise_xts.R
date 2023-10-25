#' @title normalise_xts
#'
#' @description This function normalizes an xts time series by transforming the data into standard normal
#' distribution, either on a monthly basis or for the entire time series.
#' Tranforming the data in the standard normal space can assist in visually identifying outliers.
#' Using the monthly distribution to normalise the data can be useful in cyclo-stationary processes.
#' In other cases, 'dist_period' should be set to NA.
#'
#'
#' @param ts A xts object containing the time series data.
#' @param dist_period The distribution period for normalization ('monthly' or NA).
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored when computing the statistics. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#'
#' @return A list of normalized xts time series.
#'
#' @examples
#'
#' @import xts
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

    if (dist_period == 'monthly'){
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
      u_emp = do.call(cbind, u_emp)
      z_emp = stats::qnorm(p = u_emp, mean = 0, sd = 1)
      normal_xts = xts(x = z_emp, order.by = index(ts))
    }

    out_list = c(out_list, list(normal_xts))
  }

  names(out_list) = colnames(ts_orig)
  return(out_list)
 }
