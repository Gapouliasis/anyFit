#'@title nc_ggplot
#'
#' @description This function takes a NetCDF raster file or an sxts object and generates
#' a ggplot2 raster plot. Supports RasterBricks e.g. Timeseries of rasters
#'
#' @param raster_file A raster file or an sxts object.
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
#' @importFrom ggplot2 ggplot aes geom_raster coord_equal theme_void theme element_text ggtitle labs
#' @importFrom lubridate parse_date_time
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom raster as.data.frame
#'
#' @export

nc_ggplot = function (raster_file, title = FALSE, legend.title = NA, common_legend = FALSE, viridis.option = "viridis",
                       ...) {

  if (is.sxts(raster_file)) {
    crd       <- coords(raster_file)
    vals      <- t(zoo::coredata(raster_file))
    raster_df <- as.data.frame(cbind(crd[, c("x", "y")], vals))
    dates     <- zoo::index(raster_file)
    colnames(raster_df)[3:ncol(raster_df)] <- as.character(dates)

  } else {
    raster_df  <- raster::as.data.frame(raster_file, xy = TRUE)
    raw_dates  <- colnames(raster_df)[-c(1, 2)]
    date_strs  <- gsub("X", replacement = "", x = raw_dates)
    funs <- c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "ymd H", "dmy H", "mdy H",
              "ydm H", "ymd HM", "dmy HM", "mdy HM", "ydm HM", "ymd HMS", "dmy HMS",
              "mdy HMS", "ydm HMS")
    dates <- tryCatch({
      parse_date_time(date_strs, orders = funs)
    }, warning = function(w) {
      if (conditionMessage(w) == "All formats failed to parse. No formats found.") {
        return(FALSE)
      }
    })
    if (!isFALSE(dates)) {
      colnames(raster_df)[3:ncol(raster_df)] <- as.character(dates)
    }
  }

  rownames(raster_df) <- NULL
  low_lim <- min(raster_df[, -c(1, 2)], na.rm = TRUE)
  upp_lim <- max(raster_df[, -c(1, 2)], na.rm = TRUE)

  temp_list <- list()
  for (i in 3:ncol(raster_df)) {
    temp_df <- cbind(raster_df[, c(1, 2)], raster_df[, i])
    colnames(temp_df)[3] <- colnames(raster_df)[i]
    if (common_legend == TRUE) {
      temp_plot <- ggplot(data = temp_df) +
        geom_raster(aes(x = x, y = y, fill = !!sym(colnames(temp_df)[3]))) +
        coord_equal() + theme_void() +
        viridis::scale_fill_viridis(na.value = "white", option = viridis.option, limits = c(low_lim, upp_lim)) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    } else {
      temp_plot <- ggplot(data = temp_df) +
        geom_raster(aes(x = x, y = y, fill = !!sym(colnames(temp_df)[3]))) +
        coord_equal() + theme_void() +
        viridis::scale_fill_viridis(na.value = "white", option = viridis.option) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    }
    if (title == TRUE) {
      plot_title <- colnames(temp_df)[3]
      if (startsWith(plot_title, "X")) {
        plot_title <- gsub("X", replacement = "", x = plot_title)
      }
      temp_plot <- temp_plot + ggtitle(plot_title)
    }
    if (!is.na(legend.title)) {
      temp_plot <- temp_plot + labs(fill = legend.title)
    }
    temp_list <- c(temp_list, list(temp_plot))
  }

  ncdf_plot <- patchwork::wrap_plots(plotlist = temp_list, ...)
  if (common_legend == TRUE) {
    ncdf_plot <- ncdf_plot + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  }
  return(ncdf_plot)
}
