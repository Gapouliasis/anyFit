#' @title sxts Class - Spatial eXtensible Time Series
#'
#' @description A spatial extension of xts class using attributes to store
#' spatial information. This approach is more elegant and memory-efficient.
#' Create a spatial xts (sxts) object
#'
#' @param data Matrix or data.frame of time series data
#' @param order.by POSIXct vector for time index
#' @param coords A data.frame or matrix with coordinates (x, y columns)
#' @param projection Character string specifying the coordinate reference system
#' @param ... Additional arguments passed to xts()
#' @return An sxts object
#' @export
sxts <- function(data, order.by, coords = NULL, projection = NULL) {
  obj <- xts::xts(data, order.by)  # Create an xts obj
  class(obj) <- c("sxts", class(obj))  # Set the class to inherit from xts
  attr(obj, "description") <- "This is a spatial xts (sxts) object"
  attr(obj, "elements") <- ncol(obj)  # Add an 'elements' attribute
  attr(obj, "coords") <- coords  # Add a 'coords' attribute
  # Add spatial attributes
  if (!is.null(coords)) {
    if (ncol(coords) == 2) {
      names(coords)[1:2] <- c("x", "y")
    } else {
      stop("coords must have at least 2 columns for x and y coordinates")
    }

    # Validate dimension compatibility
    if (nrow(coords) != ncol(obj)) {
      stop("Number of coordinate rows must match number of data columns")
    }

    attr(obj, "coords") <- coords
  }
  attr(obj, "projection") <- projection  # Add a 'projection' attribute
  return(obj)
}

#' @rdname sxts
#' @export
str.sxts <- function(object, ...) {
  cat("'sxts' Spatial xts object\n")
  cat(" Description: ", attr(object, "description"), "\n")
  cat(" Projection: ", attr(object, "projection"), "\n")
  cat(" Elements: ", attr(object, "elements"), " spatial locations\n")
  cat("Spatial extent is:", min(attr(object, "coords")[,1]), ",", max(attr(object, "coords")[,1]), ",",
      min(attr(object, "coords")[,2]), ",", max(attr(object, "coords")[,2]), " (xmin, xmax, ymin, ymax)", "\n")
  cat(" Time series data:\n")
  cat("Range of dates is from:", as.character(min(index(object))), "to:",  as.character(max(index(object))), "\n")
  # NextMethod("str")
}


# Custom summary method for 'sxts'
#' @rdname sxts
#' @export
summary.sxts <- function(obj, ...) {
  cat("Summary of sxts object:\n")
  cat("Number of elements:", ncol(obj), "\n")
  cat("Number of dates:", nrow(obj), "\n")
  cat("Range of dates is from:", as.character(min(index(obj))), "to:",  as.character(max(index(obj))), "\n")
  time_diff <- diff(index(obj))
  if (length(unique(time_diff)) == 1) {
    cat("Time step is:",  "Strict", "\n")
    cat("Time step is:",   as.numeric(time_diff[1], units = "hours") , "hours",  "\n")
  } else{
    cat("Time step is:",  "Variable", "\n")
    cat("Min time step is:",  as.numeric(min(time_diff[1]), units = "hours") , "hours", "\n")
    cat("Max time step is:",  as.numeric(max(time_diff[1]), units = "hours") , "hours", "\n")
  }
  cat("Number of dates:", nrow(obj), "\n")
  cat("Spatial projection:", attr(obj, "projection"), "\n")
  cat("Spatial extent is:", min(attr(obj, "coords")[,1]), ",", max(attr(obj, "coords")[,1]), ",",
      min(attr(obj, "coords")[,2]), ",", max(attr(obj, "coords")[,2]), " (xmin, xmax, ymin, ymax)", "\n")
  cat("Min value is:", min(obj),  "\n")
  cat("Max value is:", max(obj),  "\n")
  # NextMethod("summary")  # Call the default summary method for xts
}


# Print method
#' @rdname sxts
#' @export
print.sxts <- function(x, ...) {
  cat("Spatial xts (sxts) object\n")
  cat("Description:", attr(x, "description"), "\n")
  cat("Projection:", attr(x, "projection"), "\n")
  cat("Number of locations:", attr(x, "elements"), "\n")
  cat("Dimensions:", paste(dim(x), collapse = " x "), "\n")

  if (!is.null(coords(x))) {
    cat("\nSpatial extent:\n")
    coord_summary <- coords(x)[, c("x", "y")]
    cat("  X range:", range(coord_summary$x, na.rm = TRUE), "\n")
    cat("  Y range:", range(coord_summary$y, na.rm = TRUE), "\n")
  }

  cat("\nTime series data:\n")
  # Use the default xts print method for the time series part
  # NextMethod("print")
}

# Null-default operator
`%||%` <- function(x, y) if (is.null(x)) y else x


# Method to retrieve the additional attributes
#' @rdname sxts
#' @export
attributes.sxts <- function(obj) {
  list(
    elements = attr(obj, "elements"),
    coords = attr(obj, "coords"),
    projection = attr(obj, "projection")
  )
}


# Accessor functions for attributes
#' @rdname sxts
#' @export
coords <- function(x) {
  UseMethod("coords")
}

#' @rdname sxts
#' @export
coords.sxts <- function(x) {
  attr(x, "coords")
}

#' @rdname sxts
#' @export
projection <- function(x) {
  UseMethod("projection")
}

#' @rdname sxts
#' @export
projection.sxts <- function(x) {
  attr(x, "projection")
}

#' @rdname sxts
#' @export
elements <- function(x) {
  UseMethod("elements")
}

#' @rdname sxts
#' @export
elements.sxts <- function(x) {
  attr(x, "elements")
}

# Subset method that preserves sxts attributes
#' @rdname sxts
#' @export
`[.sxts` <- function(x, i, j, drop = FALSE, ...) {
  # Get the subset using xts method
  result <- NextMethod("[")

  # Preserve sxts class and attributes if still valid
  if (xts::is.xts(result)) {
    class(result) <- c("sxts", class(result))
    attr(result, "description") <- attr(x, "description")
    attr(result, "projection") <- attr(x, "projection")

    # Handle spatial subsetting
    if (!missing(j) && !is.null(coords(x))) {
      # Subset coordinates accordingly
      coords_subset <- coords(x)[j, , drop = FALSE]
      attr(result, "coords") <- coords_subset
      attr(result, "elements") <- nrow(coords_subset)
    } else {
      # No spatial subsetting, keep original coords
      attr(result, "coords") <- coords(x)
      attr(result, "elements") <- attr(x, "elements")
    }
  }

  return(result)
}


# Validation function
#' @rdname sxts
#' @export
is.sxts <- function(x) {
  inherits(x, "sxts")
}

restore_sxts <- function(result, original) {
  if (xts::is.xts(result)) {
    class(result) <- c("sxts", class(result))
    attr(result, "description") <- attr(original, "description")
    attr(result, "projection")  <- attr(original, "projection")
    attr(result, "coords")      <- attr(original, "coords")
    attr(result, "elements")    <- attr(original, "elements")
  }
  result
}

# S3 Methods (these work with NextMethod)
#' @rdname sxts
#' @export
lag.sxts <- function(x, k = 1, ...) {
  result <- NextMethod("lag")
  restore_sxts(result, x)
}

#' @rdname sxts
#' @export
diff.sxts <- function(x, lag = 1, differences = 1, ...) {
  result <- NextMethod("diff")
  restore_sxts(result, x)
}

#' @rdname sxts
#' @export
Ops.sxts <- function(e1, e2) {
  result <- NextMethod("Ops")
  if (inherits(e1, "sxts")) {
    restore_sxts(result, e1)
  } else if (inherits(e2, "sxts")) {
    restore_sxts(result, e2)
  } else {
    result
  }
}

#' Simple mask for sxts
#' @rdname sxts
#' @importFrom sf st_crs st_transform st_as_sf st_within st_union st_coordinates NA_crs_
#' @export
#'
mask.sxts <- function(obj, xlim = NA, ylim = NA, shapefile_name = NA, mask = NA) {
  obj.coords <- coords.sxts(obj)

  if (!all(is.na(xlim)) && !all(is.na(ylim))) {
    # Bounding-box masking — no CRS alignment needed
    Ix <- obj.coords[, 'x'] >= min(xlim) & obj.coords[, 'x'] <= max(xlim)
    Iy <- obj.coords[, 'y'] >= min(ylim) & obj.coords[, 'y'] <= max(ylim)
    masked_obj <- obj[, which(Ix & Iy)]

  } else {
    # Spatial mask: load from file path or accept a directly supplied sf/Spatial* object
    if (!is.na(shapefile_name)) {
      mask <- sf::st_read(shapefile_name, quiet = TRUE)
    } else if (identical(mask, NA)) {
      stop("Provide xlim/ylim for bounding-box masking, a shapefile_name, or a mask object.")
    }

    # Normalise mask to sf
    mask_sf <- if (inherits(mask, "sf")) mask else sf::st_as_sf(mask)

    # CRS alignment
    sxts_proj <- attr(obj, "projection")
    sxts_crs  <- tryCatch(
      if (!is.null(sxts_proj) && !is.na(sxts_proj)) sf::st_crs(sxts_proj) else sf::NA_crs_,
      error = function(e) sf::NA_crs_
    )
    mask_crs <- sf::st_crs(mask_sf)

    if (!is.na(sxts_crs) && !is.na(mask_crs)) {
      if (sxts_crs != mask_crs) {
        mask_sf <- sf::st_transform(mask_sf, sxts_crs)
      }
      # Point-in-polygon via sf spatial predicate (handles multipolygons correctly)
      pts_sf <- sf::st_as_sf(as.data.frame(obj.coords), coords = c("x", "y"), crs = sxts_crs)
      within <- sf::st_within(pts_sf, sf::st_union(mask_sf), sparse = FALSE)
      Imask  <- which(within)
    }

    masked_obj <- obj[, Imask]
  }

  return(masked_obj)
}



# S3 generic for sxts -> raster conversion
#' @export
rasterFromSxts <- function(x, ...) UseMethod("rasterFromSxts")

# Add a sxts to raster method
#' @rdname sxts
#' @export
rasterFromSxts.sxts <- function(obj, ...){
  coords <- attr(obj, "coords")
  proj_str <- attr(obj, "projection")
  dates <- index(obj)
  temp_matrix <- cbind(coords, t(coredata(obj)))
  new_raster <- raster::rasterFromXYZ(temp_matrix)
  raster::projection(new_raster) <- proj_str
  new_raster <- raster::setZ(new_raster, as.character(dates))
  names(new_raster) <- as.character(dates)
  return(new_raster)
}

# Convert a Raster object to an sxts object
#' @rdname sxts
#' @export
sxtsFromRaster <- function(raster, ...){
  t = raster::rasterToPoints(raster)
  tt = t(t)
  coords = t(tt[c(1, 2), ])
  tt = tt[-c(1, 2), ]
  dates = rownames(tt)
  temp_dates = gsub("X", replacement = "", x = dates)
  funs = c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "ymd H", "dmy H", "mdy H",
           "ydm H", "ymd HM", "dmy HM", "mdy HM", "ydm HM", "ymd HMS", "dmy HMS",
           "mdy HMS", "ydm HMS")

  for (tfun in funs) {
    dates = tryCatch({
      dates=parse_date_time(temp_dates, orders = tfun, truncated = 3)
    }, warning = function(w) {
    })
    if (!is.null(dates)) {
      break
    }
  }

  projection <- raster::projection(raster)

  new_sxts <- sxts(data = tt, order.by = dates, coords = coords, projection = projection)
}


#' Aggregate sxts over shapefile polygons (zonal statistics)
#'
#' @description For each polygon, finds all sxts spatial points that fall
#' within it and aggregates them row-wise (per time step) using \code{FUN}.
#' Returns a plain \code{xts} with polygon names as column names. Points that
#' do not fall within any polygon are silently dropped.
#'
#' The polygon source can be supplied in three ways (mutually exclusive):
#' \itemize{
#'   \item \code{shapefile} — an \code{sf}/\code{Spatial*} object or a file path;
#'         \code{name_col} must also be given.
#'   \item \code{country} — a country name matching \code{world_data$name};
#'         yields one output column named after the country.
#'   \item \code{continent} — a continent name matching \code{world_data$continent};
#'         yields one output column per country within the continent.
#' }
#'
#' @param x An \code{sxts} object
#' @param shapefile An \code{sf} object, \code{Spatial*} object, or file path to a shapefile
#' @param name_col Character; column in \code{shapefile} whose values become output column names
#' @param country Character; country name (matches \code{world_data$name})
#' @param continent Character; continent name (matches \code{world_data$continent})
#' @param FUN Aggregation function applied row-wise within each polygon (default: \code{mean})
#' @param ... Additional arguments forwarded to \code{FUN} (e.g. \code{na.rm = TRUE})
#' @return A plain \code{xts} with one column per polygon that contains at least one sxts point
#' @importFrom sf st_read st_as_sf st_crs st_transform st_within NA_crs_
#' @export
zonal_stats <- function(x, ...) {
  UseMethod("zonal_stats")
}

#' @rdname zonal_stats
#' @export
zonal_stats.sxts <- function(x, shapefile = NULL, name_col = NULL,
                              country = NA, continent = NA, FUN = mean, ...) {
  # Resolve polygon source
  if (!is.na(country)) {
    if (!country %in% world_data$name) {
      stop("Country '", country, "' not found. Check world_data$name for valid names.")
    }
    mask_sf  <- world_data[world_data$name == country, ]
    mask_sf  <- if (inherits(mask_sf, "sf")) mask_sf else sf::st_as_sf(mask_sf)
    name_col <- "name"
  } else if (!is.na(continent)) {
    if (!continent %in% world_data$continent) {
      stop("Continent '", continent, "' not found. Check world_data$continent for valid names.")
    }
    mask_sf  <- world_data[world_data$continent == continent, ]
    mask_sf  <- if (inherits(mask_sf, "sf")) mask_sf else sf::st_as_sf(mask_sf)
    name_col <- "name"
  } else {
    if (is.null(shapefile)) {
      stop("Supply one of: shapefile, country, or continent.")
    }
    if (is.character(shapefile)) {
      shapefile <- sf::st_read(shapefile, quiet = TRUE)
    }
    mask_sf <- if (inherits(shapefile, "sf")) shapefile else sf::st_as_sf(shapefile)
    if (is.null(name_col)) {
      stop("name_col must be specified when using a shapefile.")
    }
  }

  if (!name_col %in% names(mask_sf)) {
    avail <- paste(setdiff(names(mask_sf), attr(mask_sf, "sf_column")), collapse = ", ")
    stop(sprintf("Column '%s' not found in shapefile. Available columns: %s", name_col, avail))
  }

  obj_coords <- coords.sxts(x)
  sxts_proj  <- attr(x, "projection")
  sxts_crs   <- tryCatch(
    if (!is.null(sxts_proj) && !is.na(sxts_proj)) sf::st_crs(sxts_proj) else sf::NA_crs_,
    error = function(e) sf::NA_crs_
  )
  mask_crs <- sf::st_crs(mask_sf)
  if (!is.na(sxts_crs) && !is.na(mask_crs) && sxts_crs != mask_crs) {
    mask_sf <- sf::st_transform(mask_sf, sxts_crs)
  }

  pts_sf <- if (!is.na(sxts_crs)) {
    sf::st_as_sf(as.data.frame(obj_coords), coords = c("x", "y"), crs = sxts_crs)
  } else {
    sf::st_as_sf(as.data.frame(obj_coords), coords = c("x", "y"))
  }

  n_polys     <- nrow(mask_sf)
  result_list <- vector("list", n_polys)
  valid_polys <- logical(n_polys)

  for (i in seq_len(n_polys)) {
    within_idx <- which(sf::st_within(pts_sf, mask_sf[i, ], sparse = FALSE)[, 1])
    if (length(within_idx) == 0L) next
    valid_polys[i]   <- TRUE
    masked           <- x[, within_idx]
    result_list[[i]] <- apply(masked, 1, FUN, ...)
  }

  valid_idx <- which(valid_polys)
  if (length(valid_idx) == 0L) {
    stop("No sxts points fall within any polygon.")
  }

  result_matrix           <- do.call(cbind, result_list[valid_idx])
  colnames(result_matrix) <- as.character(mask_sf[[name_col]][valid_idx])
  return(xts::xts(result_matrix, order.by = zoo::index(x)))
}