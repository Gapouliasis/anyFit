normalise_xts = function(ts, dist_period = 'monthly',
                         ignore_zeros = FALSE, zero_threshold = 0.01){

  col_list = list()
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

    col_list = c(col_list, list(normal_xts))
  }

  normal_xts = do.call(cbind.xts, col_list)
}
