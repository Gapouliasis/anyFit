#' @noRd
#' @description Strips a leading \code{"X"} prefix from each name and attempts
#'   to parse the result as a date via \code{\link[lubridate]{parse_date_time}}.
#'   If all names parse to dates, the date strings are returned as layer titles;
#'   otherwise the original names are kept verbatim. This heuristic distinguishes
#'   between timeseries layers (where R prepends \code{"X"} to date column
#'   names) and statistical/parameter layers (e.g. \code{"mean"},
#'   \code{"shape1"}) that should be displayed as-is.
#'
#' @param names A character vector of layer names.
#'
#' @return A character vector of display labels.
#'
#' @importFrom lubridate parse_date_time
# Relabel value-column names as dates when they parse as such; otherwise keep
# them as-is (e.g. statistic/parameter names from basic_stats_nc / fitlm_nc).
.layer_titles <- function(names) {
  date_strs <- gsub("X", replacement = "", x = names)
  funs <- c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "ymd H", "dmy H", "mdy H",
            "ydm H", "ymd HM", "dmy HM", "mdy HM", "ydm HM", "ymd HMS", "dmy HMS",
            "mdy HMS", "ydm HMS")
  dates <- suppressWarnings(parse_date_time(date_strs, orders = funs))
  if (all(is.na(dates))) names else as.character(dates)
}

#' @title nc_ggplot
#'
#' @description Produces ggplot2 raster maps from an sxts object or a Raster*
#'   object (RasterLayer, RasterStack, or RasterBrick). The function handles
#'   both timeseries layers — where individual time steps are mapped as
#'   separate panels — and statistical layers such as the parameter and
#'   goodness-of-fit rasters returned by \code{\link{basic_stats_nc}} or
#'   \code{\link{fitlm_nc}}. Layer names that parse as dates are reformatted
#'   as titles; statistic names are kept verbatim. All rasters are rendered
#'   with viridis colour scales via \pkg{ggplot2} and composited into a
#'   multi-panel layout by \pkg{patchwork}. A common legend across all panels
#'   can be enforced by setting \code{common_legend = TRUE}, which scales the
#'   fill range to the global minimum and maximum of the entire dataset.
#'
#' @param data An sxts object or a raster Raster* object (RasterLayer,
#'   RasterStack, RasterBrick).
#' @param title Logical; if \code{TRUE}, layer names are added as plot titles.
#'   Default \code{FALSE}.
#' @param legend.title Character; title for the colour legend. Use \code{NA}
#'   to omit. Default \code{NA}.
#' @param common_legend Logical; if \code{TRUE}, a single common legend is
#'   collected across all panels with fill scaled to the global data range.
#'   Default \code{FALSE}.
#' @param viridis.option Character; the viridis colour palette option
#'   (\code{"viridis"}, \code{"magma"}, \code{"plasma"}, etc.). Default
#'   \code{"viridis"}.
#' @param ... Additional arguments passed to
#'   \code{\link[patchwork]{wrap_plots}}.
#'
#' @return A ggplot2 object representing the raster map.
#'
#' @examples
#' # Synthetic sxts with 4 cells and 5 daily steps
#' set.seed(42)
#' n <- 5
#' dates <- seq(as.POSIXct("2000-01-01", tz = "UTC"), by = "day", length.out = n)
#' vals <- matrix(rnorm(n * 4, mean = 10, sd = 2), nrow = n, ncol = 4)
#' coords <- data.frame(x = c(0, 1, 0, 1), y = c(0, 0, 1, 1))
#' ts_sxts <- sxts(data = vals, order.by = dates, coords = coords,
#'                 projection = "+proj=longlat +datum=WGS84")
#'
#' nc_ggplot(ts_sxts, title = TRUE)
#'
#' # With a common legend
#' nc_ggplot(ts_sxts, title = TRUE, common_legend = TRUE,
#'           legend.title = "Value")
#'
#' @importFrom ggplot2 ggplot aes geom_raster coord_equal theme_void theme element_text ggtitle labs
#' @importFrom lubridate parse_date_time
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom raster as.data.frame
#'
#' @export

nc_ggplot = function (data, title = FALSE, legend.title = NA, common_legend = FALSE, viridis.option = "viridis",
                       ...) {

  if (is.sxts(data)) {
    crd       <- coords(data)
    vals      <- t(zoo::coredata(data))
    raster_df <- as.data.frame(cbind(crd[, c("x", "y")], vals))
    dates     <- zoo::index(data)
    colnames(raster_df)[3:ncol(raster_df)] <- as.character(dates)

  } else if (inherits(data, "Raster")) {
    raster_df <- raster::as.data.frame(data, xy = TRUE)
    vcols <- seq_len(ncol(raster_df))[-c(1, 2)]
    colnames(raster_df)[vcols] <- .layer_titles(colnames(raster_df)[vcols])

  } else {
    stop("`data` must be an sxts object or a raster Raster* object ",
         "(RasterLayer/RasterStack/RasterBrick), not ", class(data)[1], ".")
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
        geom_raster(aes(x = x, y = y, fill = !!rlang::sym(colnames(temp_df)[3]))) +
        coord_equal() + theme_void() +
        viridis::scale_fill_viridis(na.value = "white", option = viridis.option, limits = c(low_lim, upp_lim)) +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    } else {
      temp_plot <- ggplot(data = temp_df) +
        geom_raster(aes(x = x, y = y, fill = !!rlang::sym(colnames(temp_df)[3]))) +
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
