# sxts constructor
sxts <- function(data, order.by, coords = NULL, projection = NULL) {
  obj <- xts(data, order.by)  # Create an xts obj
  class(obj) <- c("sxts", class(obj))  # Set the class to inherit from xts
  attr(obj, "description") <- "This is a sxts object"
  attr(obj, "name") <- "poutses"
  attr(obj, "elements") <- ncol(obj)  # Add an 'elements' attribute
  attr(obj, "coords") <- coords  # Add a 'coords' attribute
  attr(obj, "projection") <- projection  # Add a 'projection' attribute
  return(obj)
}

# Method to retrieve the additional attributes
attributes.sxts <- function(obj) {
  list(
    elements = attr(obj, "elements"),
    coords = attr(obj, "coords"),
    projection = attr(obj, "projection")
  )
}

.makeISO8601 <- function(x) {
  paste(start(x), end(x), sep = "/")
}

str.sxts <- function(object){
  is.data.empty <- is.null(dim(object)) || sum(dim(object)) == 0
  is.zero.index <- (length(.index(object)) == 0)

  nr <- NROW(object)
  nc <- ifelse(is.data.empty, 0, NCOL(object))

  # "zero-length" xts
  #    * index length == 0, but tclass and tzone are set
  #    * NROW == 0
  #    * NCOL >  0 and may have column names
  # examples:
  #   str(.xts(1, 1)["1900"])
  #   str(.xts(cbind(a = 1, b = 2), 1)["1900"])
  is.zero.length <- (is.zero.index && nr == 0 && !is.data.empty)

  # "zero-width" xts
  #    * index length > 0
  #    * NROW == 0
  #    * NCOL == 0
  # example:
  #   str(.xts(, 1:5))
  is.zero.width <- (!is.zero.index && is.data.empty)

  # "empty" xts
  #    * index length == 0, but tclass and tzone are set
  #    * NROW == 0
  #    * NCOL == 0
  # example:
  #   str(.xts(, numeric(0)))
  #   str(.xts(matrix()[0,0], numeric(0)))
  is.empty <- (is.zero.index && is.data.empty)

  if (is.empty) {
    header <- "An empty sxts object"
  } else if (is.zero.length) {
    header <- "A zero-length sxts object"
  } else {
    # zero-width and regular xts objects
    if (is.zero.width) {
      header <- "A zero-width sxts object on"
    } else {
      header <- "An sxts object on"
    }
    time.range <- sub("/", " / ", .makeISO8601(object), fixed = TRUE)
    header <- paste(header, time.range)
  }

  cat(header, "\n")
}





# Custom summary method for 'sxts'
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

# # Update the print method to display additional attributes
# print.sxts <- function(obj, ...) {
#   cat("This is a spatial xts (xts) obj:\n")
#   cat("Coordinates:", attr(obj, "coordinates"), "\n")
#   cat("Projection:", attr(obj, "projection"), "\n")
#   NextMethod("print")  # Call the default print method for xts
# }

# Add a sxts to raster method
rasterFromSxts.sxts <- function(obj, ...){
  coords <- attr(obj, "coords")
  projection <- attr(obj, "projection")
  temp_matrix <- cbind(coords, t(coredata(obj)))
  new_raster <- raster::rasterFromXYZ(temp_matrix)
  projection(new_raster) <- projection
  new_raster <- raster::setZ(new_raster, as.character(dates))
  names(new_raster) <- as.character(dates)
  return(new_raster)
}

# Add a raster to sxts method
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


