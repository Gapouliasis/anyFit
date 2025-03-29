#' @title nc2xts
#'
#' @description This function reads a NetCDF raster file, optionally masks it by country, continent, or user defined shapefile,
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
#'
#' @export

nc2xts = function (filename, varname, shapefile = NA, country = NA, continent = NA) {
  nc_data = ncdf4::nc_open(filename = filename)
  nc_brick = raster::brick(filename, varname = varname, level = 1)
  if (!is.na(country)) {
    countries = world_data$name
    if (country %in% countries) {
      mask = subset(world_data, name == country)
      r2 <- raster::crop(nc_brick, raster::extent(mask))
      r3 <- raster::mask(r2, mask = mask)
    }
    else {
      stop("Country name is incorrect")
    }
  }
  else if (!is.na(continent)) {
    continents = world_data$continent
    if (continent %in% continents) {
      mask = subset(world_data, continent == continent)
      r2 <- raster::crop(nc_brick, raster::extent(mask))
      r3 <- raster::mask(r2, mask = mask)
    }
    else {
      stop("Continent name is incorrect")
    }
  }
  else if (!is.na(shapefile)) {
    mask = raster::shapefile(shapefile)
    r2 <- raster::crop(nc_brick, raster::extent(mask))
    r3 <- raster::mask(r2, mask = mask)
  }
  else {
    r3 = nc_brick
  }
  t = raster::rasterToPoints(r3)
  tt = t(t)
  coords = t(tt[c(1, 2), ])
  tt = tt[-c(1, 2), ]

  var_names = names(nc_data$dim)

  time_var = NULL
  for (var in var_names) {
    attr_units = ncdf4::ncatt_get(nc_data, var, "units")$value
    if (!is.null(attr_units) && grepl("since", attr_units)) {
      time_var = var
      base_date = as.POSIXct(sub(".*since ", "", attr_units), tz = "UTC")
      base_time = sub(" since.*", "", attr_units)
      break
    }
  }

  # Display the identified time variable
  # print(paste0("Time variable name is:", time_var))
  raw_dates = ncdf4::ncvar_get(nc_data, time_var)
  dates = base_date + do.call(base_time, list(raw_dates))

  r4 = raster::setZ(r3, as.character(dates))
  names(r4) = as.character(dates)
  ncdf_sxts <- sxts(data = tt, order.by = dates, coords = coords, projection = raster::projection(r4))
  list_out = list(raster = r4, ncdf_sxts = ncdf_sxts)
  return(list_out)
}



