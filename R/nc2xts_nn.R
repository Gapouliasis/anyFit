#' @title nc2xts_nn
#'
#' @description Extracts a point (nearest-neighbour) time series for one variable from one or
#' more NetCDF files at the supplied coordinates and returns a plain \code{xts}. No
#' interpolation/extrapolation is performed: the value of the grid cell nearest to each point is
#' returned. Points falling outside the grid extent yield an all-\code{NA} column.
#'
#' Data are read directly with \code{ncdf4} (no \code{raster} intermediate): only the requested
#' cells are read from disk via an indexed hyperslab, and dates are parsed from the CF time units.
#'
#' @param filename Path to the NetCDF file. May be a character vector of multiple files (e.g. data
#' downloaded in time chunks); the results are row-bound and ordered by date.
#' @param varname Name of the variable to extract (a single variable name).
#' @param coords A matrix or data.frame of coordinates; the first two columns are taken as
#' x (longitude) and y (latitude).
#'
#' @return An xts time series object with one column per supplied coordinate.
#'
#' @examples
#'
#' rnd_xts = nc2xts_nn(filename = filename, varname = "rr", coords = coords)
#'
#' @importFrom xts xts rbind.xts
#' @importFrom ncdf4 nc_open nc_close ncatt_get ncvar_get
#'
#' @export

nc2xts_nn <- function(filename, varname, coords) {
  coords <- as.matrix(coords)[, 1:2, drop = FALSE]

  # Read the point series from a single NetCDF file.
  get_xts <- function(file) {
    nc_data <- ncdf4::nc_open(file)
    on.exit(ncdf4::nc_close(nc_data))

    # --- Identify lon / lat / time dims from the variable's own dim list ---
    var_dims <- nc_data$var[[varname]]$dim
    if (is.null(var_dims)) {
      stop("Variable '", varname, "' not found in ", file)
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

    if (is.null(time_idx) || is.null(lon_idx) || is.null(lat_idx)) {
      stop("Could not identify lon/lat/time dimensions for variable '", varname,
           "'. Identified: lon=", lon_var, " lat=", lat_var, " time=", time_var)
    }

    # --- Read coordinates and dates (cheap: 1D axes only) ---
    lons      <- ncdf4::ncvar_get(nc_data, lon_var)
    lats      <- ncdf4::ncvar_get(nc_data, lat_var)
    raw_times <- ncdf4::ncvar_get(nc_data, time_var)
    dates     <- base_date + get(base_time, envir = asNamespace("lubridate"))(raw_times)
    n_time    <- length(dates)

    # Cell extent (half a cell beyond the outermost centres) for the in-grid test
    rx <- stats::median(abs(diff(lons))); ry <- stats::median(abs(diff(lats)))
    x_lo <- min(lons) - rx / 2; x_hi <- max(lons) + rx / 2
    y_lo <- min(lats) - ry / 2; y_hi <- max(lats) + ry / 2

    # --- Read each point's nearest-cell series via an indexed hyperslab ---
    n_pts <- nrow(coords)
    mat   <- matrix(NA_real_, nrow = n_time, ncol = n_pts)
    for (p in seq_len(n_pts)) {
      px <- coords[p, 1]; py <- coords[p, 2]
      if (px < x_lo || px > x_hi || py < y_lo || py > y_hi) next  # outside grid -> NA

      ix <- which.min(abs(lons - px))
      iy <- which.min(abs(lats - py))

      start <- rep(1L, length(var_dims))
      count <- vapply(var_dims, function(d) d$len, integer(1))
      start[lon_idx] <- ix; count[lon_idx] <- 1L
      start[lat_idx] <- iy; count[lat_idx] <- 1L
      cell <- ncdf4::ncvar_get(nc_data, varname, start = start, count = count,
                               collapse_degen = FALSE)
      mat[, p] <- as.numeric(cell)
    }

    xts::xts(mat, order.by = dates)
  }

  out_list <- lapply(filename, get_xts)
  out_xts  <- if (length(out_list) == 1L) out_list[[1]] else do.call(xts::rbind.xts, out_list)

  if (ncol(out_xts) == 1L) {
    colnames(out_xts) <- varname
  } else {
    colnames(out_xts) <- paste0("Point", seq_len(ncol(out_xts)))
  }
  out_xts
}
