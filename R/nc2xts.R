#' @title nc2xts
#'
#' @description Reads a NetCDF file and converts it into an sxts object.
#' Optionally masks the result by country, continent, or a user-supplied shapefile.
#' Reads data directly with ncdf4 — no raster intermediate — making it
#' substantially faster and more memory-efficient for large files.
#'
#' The CRS is detected automatically from the NetCDF metadata (grid_mapping
#' variable attributes or global crs_wkt/spatial_ref attributes). If no CRS
#' is found in the file the \code{projection} argument must be supplied.
#' When both sources are present they are compared and an error is raised if
#' they disagree.
#'
#' @param filename Path to the NetCDF file.
#' @param varname Name of the variable to extract.
#' @param shapefile (optional) Path to a shapefile for masking.
#' @param country (optional) Country name for masking (matches world_data$name).
#' @param continent (optional) Continent name for masking (matches world_data$continent).
#' @param projection (optional) CRS string (PROJ4 or WKT) used as a fallback when
#'   the NetCDF file contains no projection metadata.  If the file does contain CRS
#'   metadata both values are compared and the function errors on disagreement.
#'   Defaults to \code{NA}.
#'
#' @return A list with element \code{ncdf_sxts}: an sxts object with rows = time
#'   steps and columns = spatial locations.
#'
#' @examples
#' result <- nc2xts(filename = "data/ncfile.nc", varname = "tp",
#'                  projection = "+proj=longlat +datum=WGS84")
#' result <- nc2xts(filename = "data/ncfile.nc", varname = "tp", country = "France",
#'                  projection = "+proj=longlat +datum=WGS84")
#'
#' @importFrom ncdf4 nc_open nc_close ncatt_get ncvar_get
#' @importFrom sf st_crs st_read
#'
#' @export

nc2xts <- function(filename, varname, shapefile = NA, country = NA, continent = NA,
                   projection = "+proj=longlat +datum=WGS84") {
  nc_data <- ncdf4::nc_open(filename)
  on.exit(ncdf4::nc_close(nc_data))

  # --- 1. Identify lon / lat / time dimensions from the variable's own dim list ---
  var_dims <- nc_data$var[[varname]]$dim
  if (is.null(var_dims)){
    stop("Variable '", varname, "' not found in ", filename)
  }

  lon_idx <- lat_idx <- time_idx <- NULL
  lon_var <- lat_var <- time_var <- NULL
  base_date <- base_time <- NULL

  for (i in seq_along(var_dims)) {
    dim_name <- var_dims[[i]]$name
    units    <- ncdf4::ncatt_get(nc_data, dim_name, "units")$value
    axis     <- ncdf4::ncatt_get(nc_data, dim_name, "axis")$value
    sname    <- ncdf4::ncatt_get(nc_data, dim_name, "standard_name")$value
    dname    <- tolower(dim_name)

    is_time <- (!is.null(units)  && grepl("since", units, fixed = TRUE)) ||
               (!is.null(axis)   && axis  == "T") ||
               (!is.null(sname)  && sname == "time") ||
               dname %in% c("time", "t")

    is_lon  <- (!is.null(axis)   && axis  == "X") ||
               (!is.null(sname)  && grepl("longitude", sname)) ||
               dname %in% c("lon", "longitude", "x", "rlon")

    is_lat  <- (!is.null(axis)   && axis  == "Y") ||
               (!is.null(sname)  && grepl("latitude", sname)) ||
               dname %in% c("lat", "latitude", "y", "rlat")

    if (is_time) {
      time_idx  <- i
      time_var  <- dim_name
      base_date <- as.POSIXct(sub(".*since ", "", units), tz = "UTC")
      base_time <- sub(" since.*", "", units)
    }
    if (is_lon) { lon_idx <- i; lon_var <- dim_name }
    if (is_lat) { lat_idx <- i; lat_var <- dim_name }
  }

  if (is.null(time_idx) || is.null(lon_idx) || is.null(lat_idx)){
    stop("Could not identify lon/lat/time dimensions for variable '", varname,
         "'. Identified: lon=", lon_var, " lat=", lat_var, " time=", time_var)
  }

  # --- 2. Read coordinates and dates ---
  lons      <- ncdf4::ncvar_get(nc_data, lon_var)
  lats      <- ncdf4::ncvar_get(nc_data, lat_var)
  raw_times <- ncdf4::ncvar_get(nc_data, time_var)
  dates     <- base_date + get(base_time, envir = asNamespace("lubridate"))(raw_times)

  # --- 3. Read variable array; ncdf4 automatically applies scale/offset and _FillValue->NA ---
  data_array <- ncdf4::ncvar_get(nc_data, varname)

  # --- 4. Permute to [time, lon, lat] then reshape to [n_time x n_locations] ---
  # aperm indices reference the position of each dim in var_dims (i.e. in data_array)
  data_array  <- aperm(data_array, c(time_idx, lon_idx, lat_idx))
  n_time      <- length(dates)
  data_matrix <- matrix(data_array, nrow = n_time, ncol = length(lons) * length(lats))

  # expand.grid(lons, lats) produces rows in lon-varies-fastest order,
  # matching the column order of data_matrix after the aperm above
  coords <- expand.grid(x = lons, y = lats)

  # --- 5. Detect CRS from NetCDF metadata ---
  nc_crs_string <- NULL

  # Check grid_mapping attribute on the data variable (CF conventions)
  gm_att <- ncdf4::ncatt_get(nc_data, varname, "grid_mapping")
  if (gm_att$hasatt) {
    gm_name <- gm_att$value
    for (crs_attr in c("crs_wkt", "spatial_ref", "proj4", "proj4_params")) {
      att <- ncdf4::ncatt_get(nc_data, gm_name, crs_attr)
      if (isTRUE(att$hasatt) && nchar(trimws(att$value)) > 0) {
        nc_crs_string <- att$value
        break
      }
    }
  }

  # Fall back to global NC attributes (GDAL / other conventions)
  if (is.null(nc_crs_string)) {
    for (crs_attr in c("crs_wkt", "spatial_ref")) {
      att <- ncdf4::ncatt_get(nc_data, 0, crs_attr)
      if (isTRUE(att$hasatt) && nchar(trimws(att$value)) > 0) {
        nc_crs_string <- att$value
        break
      }
    }
  }

  nc_crs <- if (!is.null(nc_crs_string)) {
    tryCatch(sf::st_crs(nc_crs_string), error = function(e) NULL)
  } else {
    NULL
  }

  # --- 6. Cross-check and resolve CRS ---
  is_manual <- !is.na(projection)

  if (is.null(nc_crs) && !is_manual) {
    stop(
      "Could not detect CRS from '", filename, "'. ",
      "Supply the projection manually via the 'projection' argument ",
      "(e.g. projection = \"+proj=longlat +datum=WGS84\")."
    )
  } else if (!is.null(nc_crs) && is_manual) {
    manual_crs <- tryCatch(sf::st_crs(projection), error = function(e) NULL)
    if (!is.null(manual_crs) && nc_crs != manual_crs) {
      stop(
        "CRS mismatch: NetCDF reports '", nc_crs_string,
        "' but the 'projection' argument is '", projection, "'."
      )
    }
    resolved_crs <- projection
  } else if (!is.null(nc_crs)) {
    resolved_crs <- if (!is.na(nc_crs$wkt)) nc_crs$wkt else nc_crs_string
  } else {
    resolved_crs <- projection
  }

  # --- 7. Build sxts ---
  ncdf_sxts <- sxts(
    data       = data_matrix,
    order.by   = dates,
    coords     = coords,
    projection = resolved_crs
  )

  # --- 8. Apply masking ---
  if (!is.na(country)) {
    if (!country %in% world_data$name){
      stop("Country name '", country, "' not found. Check world_data$name for valid names.")
    }
    mask      <- world_data[world_data$name == country, ]
    ncdf_sxts <- mask.sxts(ncdf_sxts, mask = mask)

  } else if (!is.na(continent)) {
    if (!continent %in% world_data$continent){
      stop("Continent name '", continent, "' not found. Check world_data$continent for valid names.")
    }
    mask      <- world_data[world_data$continent == continent, ]
    ncdf_sxts <- mask.sxts(ncdf_sxts, mask = mask)

  } else if (!is.na(shapefile)) {
    mask      <- sf::st_read(shapefile, quiet = TRUE)
    ncdf_sxts <- mask.sxts(ncdf_sxts, mask = mask)
  }

  list(ncdf_sxts = ncdf_sxts)
}
