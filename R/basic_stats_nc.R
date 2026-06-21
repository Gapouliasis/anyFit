#' Compute basic statistics for gridded (NetCDF or raster) data
#'
#' @description
#' Computes a comprehensive set of over 30 summary statistics for every grid
#' cell in a spatial time series. The function accepts an \code{sxts} object,
#' a \code{Raster*} object, or a NetCDF file path with variable name. When
#' \code{filename} and \code{varname} are supplied, the data are first imported
#' via \code{\link{nc2xts}} and additional arguments are forwarded to that
#' function. For raster inputs the layer names are parsed for dates, 
#' and the raster values are converted to an xts matrix. The
#' returned raster has one layer per statistic; layer names correspond to the
#' rows of the \code{basic_stats} output.
#'
#' @param data An \code{sxts} object, or a \code{Raster*} object. Leave as
#'   \code{NULL} when supplying \code{filename} and \code{varname} instead.
#' @param filename A NetCDF file name to import. Ignored when \code{data} is
#'   supplied directly.
#' @param varname The name of the variable to extract from \code{filename}.
#' @param ignore_zeros Logical. If \code{TRUE}, values at or below
#'   \code{zero_threshold} are excluded from distributional statistics.
#'   Default is \code{FALSE}.
#' @param zero_threshold Numeric threshold below which values are treated as
#'   zero. Default is \code{0.01}.
#' @param ... Additional arguments passed to \code{\link{nc2xts}} when
#'   \code{filename} and \code{varname} are provided.
#'
#' @return A \code{RasterBrick} with one layer per statistic. The layer names
#'   correspond to the statistic labels (e.g. \code{"Mean"}, \code{"StDev"},
#'   \code{"Skewness"}, \code{"Pdr"}, etc.). The spatial reference (projection
#'   and extent) is inherited from the input.
#'
#' @examples
#' \dontrun{
#' # From an sxts object (e.g. output of nc2xts)
#' s <- nc2xts(filename = "precip.nc", varname = "pr")
#' r <- basic_stats_nc(data = s, ignore_zeros = TRUE)
#' raster::plot(r[["Mean"]])
#'
#' # Direct from NetCDF (shortcut — nc2xts is called internally)
#' r <- basic_stats_nc(filename = "precip.nc", varname = "pr",
#'                      ignore_zeros = TRUE)
#' }
#'
#' @importFrom lubridate parse_date_time
#' @importFrom xts xts
#'
#' @export
#'
basic_stats_nc = function(data = NULL, filename = NA, varname = NA,
                    ignore_zeros = FALSE, zero_threshold = 0.01, ...){

  if (!is.na(filename)) {
    ncdf_sxts <- nc2xts(filename = filename, varname = varname, ...)
    coords    <- coords(ncdf_sxts)
    proj_str  <- projection(ncdf_sxts)

  } else if (is.sxts(data)) {
    ncdf_sxts <- data
    coords    <- coords(ncdf_sxts)
    proj_str  <- projection(ncdf_sxts)

  } else {
    t         <- raster::rasterToPoints(data)
    tt        <- t(t)
    coords    <- t(tt[c(1, 2), ])
    tt        <- tt[-c(1, 2), ]
    dates     <- rownames(tt)
    temp_dates <- gsub("X", replacement = "", x = dates)
    funs <- c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "ymd H", "dmy H", "mdy H",
              "ydm H", "ymd HM", "dmy HM", "mdy HM", "ydm HM", "ymd HMS", "dmy HMS",
              "mdy HMS", "ydm HMS")
    dates     <- parse_date_time(temp_dates, orders = funs)
    ncdf_sxts <- xts(x = tt, order.by = dates)
    proj_str  <- raster::projection(data)
  }

  ncdf_stats <- t(as.matrix(basic_stats(ncdf_sxts, ignore_zeros = ignore_zeros, zero_threshold = zero_threshold, plot = FALSE)$stats_table))

  ncdf_stats  <- cbind(coords, ncdf_stats)
  raster_stats <- raster::rasterFromXYZ(ncdf_stats)
  raster::projection(raster_stats) <- proj_str

  return(raster_stats)
}
