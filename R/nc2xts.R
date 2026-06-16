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
#' @param xlim (optional) Numeric length-2 vector giving the longitude/x range to
#'   read, expressed in the grid's own coordinate units. Used as an alternative to a
#'   polygon mask to restrict the read to a bounding box. Requires \code{ylim}.
#' @param ylim (optional) Numeric length-2 vector giving the latitude/y range to
#'   read, expressed in the grid's own coordinate units. Requires \code{xlim}.
#' @param projection (optional) CRS string (PROJ4 or WKT) used as a fallback when
#'   the NetCDF file contains no projection metadata.  If the file does contain CRS
#'   metadata both values are compared and the function errors on disagreement.
#'   Defaults to \code{NA}.
#'
#' @details For large files, supplying a \code{country}, \code{continent},
#'   \code{shapefile}, or \code{xlim}/\code{ylim} restricts the NetCDF read to the
#'   bounding box of the requested region, so only that hyperslab is loaded from disk
#'   rather than the whole grid. For polygon masks the bounding box is computed in the
#'   grid's own CRS, so the prefilter is correct for both geographic and projected
#'   grids. The prefilter is skipped (full read, with a warning) for rotated-pole grids.
#'
#' @return An sxts object with rows = time steps and columns = spatial locations.
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
                   xlim = NA, ylim = NA,
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
  rotated   <- FALSE

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

    # Rotated-pole grids: 1D axes are in rotated coordinates, so a geographic/projected
    # bounding box cannot be mapped to a contiguous index window (see prefilter below).
    if (dname %in% c("rlon", "rlat") ||
        (!is.null(sname) && sname %in% c("grid_longitude", "grid_latitude"))) {
      rotated <- TRUE
    }

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

  # --- 2. Read coordinates and dates (cheap: 1D axes only) ---
  lons      <- ncdf4::ncvar_get(nc_data, lon_var)
  lats      <- ncdf4::ncvar_get(nc_data, lat_var)
  raw_times <- ncdf4::ncvar_get(nc_data, time_var)
  dates     <- base_date + get(base_time, envir = asNamespace("lubridate"))(raw_times)

  # --- 3. Detect CRS from NetCDF metadata (needed before the read for the bbox transform) ---
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

  # --- 4. Resolve masking geometry (validate names; obtain an sf mask) ---
  mask_sf <- NULL
  if (!is.na(country)) {
    if (!country %in% world_data$name){
      stop("Country name '", country, "' not found. Check world_data$name for valid names.")
    }
    mask_sf <- world_data[world_data$name == country, ]
  } else if (!is.na(continent)) {
    if (!continent %in% world_data$continent){
      stop("Continent name '", continent, "' not found. Check world_data$continent for valid names.")
    }
    mask_sf <- world_data[world_data$continent == continent, ]
  } else if (!is.na(shapefile)) {
    mask_sf <- sf::st_read(shapefile, quiet = TRUE)
  }

  # --- 5. Determine the lon/lat index window to read ---
  # Default: full extent. Narrow it from the mask bbox (transformed into the grid's CRS)
  # or from an explicit xlim/ylim. The prefilter is skipped for rotated-pole grids.
  lon_keep <- seq_along(lons)
  lat_keep <- seq_along(lats)

  xr <- yr <- NULL
  if (!is.null(mask_sf)) {
    if (rotated) {
      warning("Rotated-pole grid detected; bounding-box prefilter skipped and the full ",
              "grid will be read. Masking results for rotated grids may be unreliable.")
    } else {
      bb_mask <- tryCatch(
        sf::st_transform(sf::st_as_sf(mask_sf), sf::st_crs(resolved_crs)),
        error = function(e) sf::st_as_sf(mask_sf)
      )
      bb <- sf::st_bbox(bb_mask)
      xr <- c(bb[["xmin"]], bb[["xmax"]])
      yr <- c(bb[["ymin"]], bb[["ymax"]])
    }
  } else if (!all(is.na(xlim)) && !all(is.na(ylim))) {
    xr <- range(xlim)
    yr <- range(ylim)
  }

  if (!is.null(xr)) {
    lon_keep <- which(lons >= xr[1] & lons <= xr[2])
    lat_keep <- which(lats >= yr[1] & lats <= yr[2])
    if (length(lon_keep) == 0L || length(lat_keep) == 0L) {
      stop("The requested region does not overlap the grid extent.")
    }
    # 1-cell pad each side for floating-point safety (coords are monotonic -> contiguous run)
    lon_keep <- seq.int(max(1L, min(lon_keep) - 1L), min(length(lons), max(lon_keep) + 1L))
    lat_keep <- seq.int(max(1L, min(lat_keep) - 1L), min(length(lats), max(lat_keep) + 1L))
  }

  lons <- lons[lon_keep]
  lats <- lats[lat_keep]

  # NetCDF axes are commonly stored as floats (e.g. NClimGrid's 1/24-deg grid), so the
  # cell spacing wobbles at ~1e-5 and downstream rasterFromXYZ rejects the grid as
  # irregular. Snap each axis back onto its regular lattice (cheap: 1-D axes only).
  # Matches the clean coordinates the previous raster::brick-based reader produced.
  regularize_axis <- function(v) if (length(v) < 2) v else v[1] + (seq_along(v) - 1) * mean(diff(v))
  lons <- regularize_axis(lons)
  lats <- regularize_axis(lats)

  # --- 6. Read only the requested hyperslab; ncdf4 applies scale/offset and _FillValue->NA ---
  start <- rep(1L, length(var_dims))
  count <- vapply(var_dims, function(d) d$len, integer(1))
  start[lon_idx] <- min(lon_keep); count[lon_idx] <- length(lon_keep)
  start[lat_idx] <- min(lat_keep); count[lat_idx] <- length(lat_keep)
  data_array <- ncdf4::ncvar_get(nc_data, varname, start = start, count = count,
                                 collapse_degen = FALSE)

  # --- 7. Permute to [time, lon, lat] then reshape in place to [n_time x n_locations] ---
  # aperm indices reference the position of each dim in var_dims (i.e. in data_array).
  # After the aperm the buffer is contiguous in column-major [time, lon, lat] order, so a
  # plain dim<- reshapes to the matrix without a second full-array copy.
  data_array <- aperm(data_array, c(time_idx, lon_idx, lat_idx))
  n_time     <- length(dates)
  dim(data_array) <- c(n_time, length(lons) * length(lats))

  # expand.grid(lons, lats) produces rows in lon-varies-fastest order,
  # matching the column order of data_array after the aperm above
  coords <- expand.grid(x = lons, y = lats)

  # --- 8. Build sxts ---
  ncdf_sxts <- sxts(
    data       = data_array,
    order.by   = dates,
    coords     = coords,
    projection = resolved_crs
  )

  # --- 9. Trim to the exact requested region (drops the read's safety-pad cells).
  #         Cheap: runs over the slab's points only, not the whole grid. ---
  if (!is.null(mask_sf)) {
    ncdf_sxts <- mask.sxts(ncdf_sxts, mask = mask_sf)
  } else if (!all(is.na(xlim)) && !all(is.na(ylim))) {
    ncdf_sxts <- mask.sxts(ncdf_sxts, xlim = xlim, ylim = ylim)
  }

  return(ncdf_sxts)
}
