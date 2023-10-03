#' @title nc2xts
#'
#' This function reads a NetCDF raster file, optionally masks it by country, continent, or user defined shapefile,
#' and converts it into an xts object for further analysis.
#'
#' @param filename The path to the NetCDF file.
#' @param varname The name of the variable to extract from the NetCDF file.
#' @param shapefile (optional) The path to a shapefile for masking.
#' @param country (optional) The name of a country for masking.
#' @param continent (optional) The name of a continent for masking.
#'
#' @return A list containing the raster, xts object, dates, and coordinates.
#'
#' @examples
#' # Load required libraries
#' library(ncdf4)
#' library(raster)
#' library(xts)
#'
#' # Example usage 1: Convert NetCDF raster to xts object
#' nc_file <-"data/ncfile.nc"
#' result <- nc2xts(filename = nc_file, varname = "tp")
#' nc_ggplot(result$raster)
#'
#' # Example usage 2: Spatial subsetting by country
#' nc_file <-"data/ncfile.nc"
#' result <- nc2xts(filename = nc_file, varname = "tp", country = "France")
#' nc_ggplot(result$raster)
#'
#' @import ncdf4
#' @import raster
#' @import xts
#' @importFrom lubridate ymd ydm mdy myd dmy dym ymd_h dmy_h mdy_h ydm_h ymd_hm dmy_hm mdy_hm ydm_hm ymd_hms dmy_hms mdy_hms ydm_hms
#' @importFrom rgdal readOGR
#'
#' @export

#filename = "C:/Users/Admin/Downloads/rr_ens_mean_0.25deg_reg_2011-2022_v27.0e.nc"
#filename = "C:/Users/Admin/Documents/ERA5/ERA5_2022_TotalPrecipitation_Europe.nc"

nc2xts = function(filename, varname, shapefile = NA, country = NA, continent = NA){
  nc_data = ncdf4::nc_open(filename = filename)

  #names(nc_data$var)

  nc_brick = raster::brick(filename, varname = varname, level = 1)

  world_data = readRDS("data/world_data.rds")
  if (!is.na(country)){
    countries = world_data$name
    if (country %in% countries){
      mask = subset(world_data, name==country)
      r2 <- raster::crop(nc_brick, raster::extent(mask))
      r3 <- raster::mask(r2, mask = mask)
    }else{
      stop("Country name is incorrect")
    }
  }else if(!is.na(continent)){
    continents = world_data$continent
    if (continent %in% continents){
      mask = subset(world_data, continent==continent)
      r2 <- raster::crop(nc_brick, raster::extent(mask))
      r3 <- raster::mask(r2, mask = mask)
    }else{
      stop("Continent name is incorrect")
    }

  }else if(!is.na(continent)){
    mask = rgdal::readOGR(shapefile)
    r2 <- raster::crop(nc_brick, raster::extent(mask))
    r3 <- raster::mask(r2, mask = mask)
  }else{
    r3 = nc_brick
  }

  t = raster::rasterToPoints(r3)
  tt = t(t)
  coords = t(tt[c(1,2),])
  tt = tt[-c(1,2),]

  dates = rownames(tt)
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

  ncdf_xts = xts(x = tt,order.by = dates)

  list_out = list(raster = r3, ncdf_xts = ncdf_xts, dates = dates, coordinates = coords)

  return(list_out)
}




