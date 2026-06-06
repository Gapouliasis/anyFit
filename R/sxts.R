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
  cat("Spatial extent is:", min(attr(obj, "coords")[,1]), ",", max(attr(obj, "coords")[,1]), ",",
      min(attr(obj, "coords")[,2]), ",", max(attr(obj, "coords")[,2]), " (xmin, xmax, ymin, ymax)", "\n")
  cat(" Time series data:\n")
  cat("Range of dates is from:", as.character(min(index(obj))), "to:",  as.character(max(index(obj))), "\n")
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
  if (is.xts(result)) {
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
#' @export
#'
mask.sxts <- function(obj, xlim = NA, ylim = NA, shapefile_name = NA, mask = NA){
  obj.coords <- coords.sxts(obj)

  if (is.na(shapefile)){
    Ix <- obj.coords[,'x'] >= min(xlim) & obj.coords[,'x'] <= max(xlim)
    Iy <- obj.coords[,'y'] >= min(ylim) & obj.coords[,'y'] <= max(ylim)
    Imask <- which(Ix & Iy)
    masked_obj <- obj[,Imask]
  }else{
    if (!is.na(shapefile_name)){
      mask = raster::shapefile(shapefile_name)
    }
    mask.xy <- raster::geom(mask)
    Imask <- sp::point.in.polygon(point.x = obj.coords[,'x'], point.y = obj.coords[,'y'],
                                  pol.x = mask.xy[, 'x'], pol.y = mask.xy[, 'y'])
    masked_obj <- obj[,which(Imask > 0)]
  }
  return(masked_obj)
}



# Add a sxts to raster method
#' @rdname sxts
#' @export
rasterFromSxts.sxts <- function(obj, ...){
  coords <- attr(obj, "coords")
  projection <- attr(obj, "projection")
  dates <- index(obj)
  temp_matrix <- cbind(coords, t(coredata(obj)))
  new_raster <- raster::rasterFromXYZ(temp_matrix)
  projection(new_raster) <- projection
  new_raster <- raster::setZ(new_raster, as.character(dates))
  names(new_raster) <- as.character(dates)
  return(new_raster)
}

# Add a raster to sxts method
#' @rdname sxts
#' @export
sxtsFromRaster.sxts <- function(raster, ...){
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


library(xts)
file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
time_zone <- "UTC"

data <- delim2xts(file_path = file_path,
                  time_zone = "UTC", delim = " ", strict_step = TRUE)

a = fitlm_dagum(data[,1], ignore_zeros = TRUE)

candidates <- list('exp','expweibull', 'gamma3', 'dagum','burr')
fits <- fitlm_multi(data[,1],candidates = candidates, ignore_zeros = TRUE)
fits$diagnostics

european_countries <- c("albania", "andorra", "austria", "belarus", "belgium", "bosnia and herzegovina", "bulgaria",
                        "croatia", "cyprus", "czechia", "denmark", "estonia", "finland", "france", "germany", "greece", "hungary",
                        "ireland", "italy", "kosovo", "latvia", "liechtenstein",  "lithuania", "luxembourg", "malta",
                        "moldova", "monaco", "montenegro", "netherlands", "north macedonia", "norway", "poland", "portugal", "romania",
                        "san marino", "serbia", "slovakia", "slovenia", "spain", "sweden", "switzerland", "ukraine", "UK", "vatican city")

library(lubridate)
library(raster)
library(xts)
library(ncdf4)
library(anyFit)
data_path <- "C:/Users/gapou/Documents/ERA5_tp"
varname <- "tp"
filename <- file.path(data_path, "ERA5_tp_1940_01.nc")
data <- nc2xts(filename = filename, varname = varname)

data_sxts <- data$ncdf_sxts

test_sxts <- data_sxts[,3:6]

nc_ggplot(data$raster[[1]]) &
  borders("world", regions = european_countries, xlim = c(-15, 40), ylim = c(32, 72))

a = mask.sxts(data_sxts, xlim = c(19,28), ylim = c(35,45))
b = rasterFromSxts.sxts(a)

shapefile_name = "C:/Users/gapou/Downloads/grc_adm_2025_ab_shp/grc_admbnda_adm0_AB.shp"
a = mask.sxts(data_sxts, shapefile = shapefile_name)
b = rasterFromSxts.sxts(a)

nc_ggplot(b[[2]]) &
  borders("world", regions = european_countries, xlim = c(-15, 40), ylim = c(32, 72))

spec_period <- endpoints(test_sxts, on = "days", k = 3)
ncdf_stats = period.apply.sxts(test_sxts, spec_period, FUN = "mean")

start_time <- Sys.time()

data_agg <- period_apply_nc(data_sxts, period = "days", period_multiplier = 3)

time_end = Sys.time()

print(time_end - start_time)



data_path <- "C:/Users/gapou/Documents/knmi-eobs"
filename <- file.path(data_path, "rr_ens_mean_0.25deg_reg_v30.0e.nc")

start_time <- Sys.time()

eobs_data = nc2xts(filename = filename, varname = "rr", country = "Italy")

time_end = Sys.time()

print(time_end - start_time)

