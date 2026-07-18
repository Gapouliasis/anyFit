#' @details
#' anyFit streamlines exploratory data analysis of point and gridded
#' hydroclimatic time series: loading (NetCDF and delimited files), summary
#' statistics, temporal aggregation, and distribution fitting by the method of
#' L-moments, all built on the \pkg{xts} time-series format and \pkg{ggplot2}
#' visualisation. Gridded data are carried in the [sxts] (spatial xts) class,
#' which stores coordinates and projection alongside the series.
#'
#' The package ships a real example dataset used throughout the vignettes and
#' examples: an E-OBS ensemble-mean daily rainfall grid over Europe
#' (0.25 degrees, 2011-2022), installed at
#' \code{system.file("extdata", "rr_ens_mean_0.25deg_reg_2011-2022_v27.0e.nc", package = "anyFit")}.
#' Load it with [nc2xts()], typically clipped to a region, e.g.
#' \code{nc2xts(f, "rr", country = "U.K. of Great Britain and Northern Ireland")}.
#' A daily station series in
#' \code{system.file("extdata", "KNMI_Daily.csv", package = "anyFit")}
#' illustrates [delim2xts()].
#'
#' @keywords internal
"_PACKAGE"
