#'@title nc_ggplot
#'
#' @description This function takes a NetCDF raster file and a variable name as inputs, and generates
#' a ggplot2 raster plot.
#'
#' @param raster_file A NetCDF raster file.
#'
#' @return A ggplot2 object representing the raster plot.
#'
#' @examples
#'
#' # Example usage:
#' # Create a ggplot2 raster plot from a NetCDF raster file
#' nc_file <- system.file("extdata/ncfile.nc", package = "raster")
#' plot <- nc_ggplot(nc_file, "temperature")
#' plot
#'
#' @import ggplot2
#' @importFrom raster as.data.frame
#' @importFrom raster raster
#'
#'
#' @export

nc_ggplot = function(raster_file){
  raster_df <- raster::as.data.frame(raster_file, xy = TRUE)
  rownames(raster_df) <- NULL
  colnames(raster_df)[3] = "variable"

  ncdf_plot = ggplot(data = raster_df) +
    geom_raster(aes_string(x = "x", y = "y", fill = "variable")) +
    scale_fill_viridis_c(na.value = "white") +
    theme_void() +
    theme(
      legend.position = "bottom"
    )
  return(ncdf_plot)
}

