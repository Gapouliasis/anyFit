#'@title nc_ggplot
#'
#' @description This function takes a NetCDF raster file and generates
#' a ggplot2 raster plot. Supports RasterBricks e.g. Timeseries of rasters
#'
#' @param raster_file A raster file.
#' @param title Logical, whether to add variable names as plot titles.
#' @param legend.title The title for the color legend.
#' @param common.legend Add a common legend to all subplots. This will scale fill according to the min and max of the entire dataset.
#' @param viridis.option The viridis color palette option (e.g., "viridis", "magma", "plasma").
#' @param ... Additional arguments to pass to 'wrap_plots' function from the 'patchwork' package.
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
#' @import viridis
#'
#' @export

nc_ggplot = function (raster_file, title = FALSE, legend.title = NA, common_legend = FALSE, viridis.option = "viridis",
                       ...) {
  raster_df <- raster::as.data.frame(raster_file, xy = TRUE)
  rownames(raster_df) <- NULL
  low_lim = min(raster_df[, -c(1,2)], na.rm = TRUE)
  upp_lim = max(raster_df[, -c(1,2)], na.rm = TRUE)
  raw_dates = colnames(raster_df)[-c(1,2)]
  dates = gsub("X", replacement = "", x = raw_dates)
  funs = c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "ymd H", "dmy H", "mdy H",
           "ydm H", "ymd HM", "dmy HM", "mdy HM", "ydm HM", "ymd HMS", "dmy HMS",
           "mdy HMS", "ydm HMS")
  flag = TRUE
  dates = tryCatch({
    dates =parse_date_time(dates, orders = funs)}, warning = function(w) {
      if ( conditionMessage(w) == "All formats failed to parse. No formats found."){
        return(FALSE)
      }})

  if (dates[1] != FALSE){
    colnames(raster_df)[3:ncol(raster_df)] = as.character(dates)
  }

  temp_list = list()
  for (i in 3:ncol(raster_df)) {
    temp_df = cbind(raster_df[, c(1, 2)], raster_df[, i])
    colnames(temp_df)[3] = colnames(raster_df)[i]
    if (common_legend == TRUE){
      temp_plot = ggplot(data = temp_df) + geom_raster(aes(x = x,
                                                           y = y, fill = !!sym(colnames(temp_df)[3]))) + coord_equal() +
        theme_void() + viridis::scale_fill_viridis(na.value = "white", option = viridis.option, limits = c(low_lim, upp_lim)) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    }else{
      temp_plot = ggplot(data = temp_df) + geom_raster(aes(x = x,
                                                           y = y, fill = !!sym(colnames(temp_df)[3]))) + coord_equal() +
        theme_void() + viridis::scale_fill_viridis(na.value = "white", option = viridis.option) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    }

    if (title == TRUE) {
      if (startsWith(colnames(temp_df)[3], "X")) {
        plot_title = gsub("X", replacement = "", x = colnames(temp_df)[3])
      }
      else {
        plot_title = colnames(temp_df)[3]
      }
      temp_plot = temp_plot + ggtitle(plot_title)
    }
    if (!is.na(legend.title)) {
      temp_plot = temp_plot + labs(fill = legend.title)
    }
    temp_list = c(temp_list, list(temp_plot))
  }
  ncdf_plot = patchwork::wrap_plots(plotlist = temp_list, ...)
  if(common_legend == TRUE){
    ncdf_plot = ncdf_plot + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  }
  return(ncdf_plot)
}

