get_xts = function(variable, filename){
  nc_brick = brick(filename, varname = variable)
  ts = t(raster::extract(nc_brick, lonlat[Imin,]))
  colnames(ts) = variable
  dates = as.POSIXct(substring(rownames(ts),2), format = "%Y.%m.%d.%H.%M.%S", tz = "UTC")
  temp_xts = xts(x = ts, order.by = dates)
  return(temp_xts)
}

for (filename in filenames_int){
  temp_month = lapply(variables, FUN = get_xts, filename = filename)
  temp_month = do.call(cbind.xts, temp_month)
  temp_list = c(temp_list, list(temp_month))
}

ERA5_xts = do.call(rbind.xts, temp_list)

saveRDS(ERA5_xts, file.path('data_out','ERA5_Integrated.rds'))
