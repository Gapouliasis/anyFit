#' @title nc_xts
#'
#' This function extracts one or more variables from NetCDF files from one or more locations at specified coordinates
#' and creates an xts time series. This function performs no interpolation/extrapolation of the data.
#' It return the values of the cell which contains the point coordinates. If the specified points fall outside of the dataset
#' it will return NANs.
#'
#' @param filename The path to the NetCDF file. It can be a character vector of multiple NetCDF files.
#' This is useful when data are downloaded in time chunks. The function will automatically order by date.
#' @param varname A character vector of variable names to extract.
#' @param coords A matrix or data.frame of coordinates (longitude and latitude) for extraction.
#'
#' @return An xts time series object containing the extracted variables.
#'
#' @examples
#'
#' rnd_xts = nc2xts_nn(filename = filename, varname = "rr", coords = coords)
#'
#' @import xts
#' @import raster
#'
#' @export

nc2xts_nn = function(filename, varname, coords){
  get_xts = function(variable, filename){
    nc_brick = raster::brick(filename, varname = variable)
    ts = t(raster::extract(nc_brick, coords, method='simple'))
    #colnames(ts) = variable
    dates = rownames(ts)
    temp_dates = gsub("X",replacement = "",x=dates)
    funs = c(ymd, ydm, mdy, myd, dmy, dym,
             ymd_h, dmy_h, mdy_h, ydm_h,
             ymd_hm, dmy_hm, mdy_hm, ydm_hm,
             ymd_hms, dmy_hms, mdy_hms, ydm_hms)
    for (tfun in funs){
      param_list = list(data = temp_dates)
      param_list$tz = 'UTC'
      dates = tryCatch({do.call(tfun,param_list)},
                       warning = function(w) {})
      if (!is.null(dates)){
        break
      }
    }
    temp_xts = xts(x = ts, order.by = dates)
    return(temp_xts)
  }

  temp_list = list()
  for (file in filename){
    temp_xts2 = lapply(varname, FUN = get_xts, filename = filename)
    temp_xts2 = do.call(cbind.xts, temp_xts2)
    if (nrow(coords) > 1){
      temp_names = lapply(paste0('Var',seq(1,length(varname))),function(x){paste(paste0('Point', seq(nrow(coords))),x,sep = '.')})
      colnames(temp_xts2) = unlist(temp_names)
    }else{
      colnames(temp_xts2) = varname
    }

    temp_list = c(temp_list, list(temp_xts2))
  }
  out_xts = do.call(rbind.xts, temp_list)
  return(out_xts)
}



