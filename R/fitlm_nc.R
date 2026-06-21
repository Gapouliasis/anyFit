#' @title Distribution Fitting on Gridded Data
#'
#' @description Fits a list of candidate distributions to every grid cell of a
#'   NetCDF raster or an \code{sxts} object via the L-moments method. The input
#'   is normalised to an \code{sxts} object, accepting \code{sxts},
#'   \code{Raster*}, or NetCDF filename and variable name. The results are
#'   returned as per-candidate \code{RasterBrick} objects for the parameters, theoretical L-moments, data
#'   L-moments, and GoF metrics. A density comparison plot of the six GoF metrics
#'   across all grid cells is assembled with \code{ggplot2}, faceted by metric
#'   and coloured by candidate distribution. When \code{parallel = TRUE}, the work is distributed across available cores
#' via the \pkg{future} framework, with either file-backed shared memory (bigmemory
#' memory-mapped matrices for single-machine parallelism) or serialised column chunks
#' (for multi-node clusters). File-backed shared memory can be enabled by \code{shared_memory = TRUE}. 
#' This is reccomended for single machine usage and especially windows for efficiency and reduced RAM consumption.
#'
#' @param data An \code{sxts} object, or a \code{Raster*} object. Leave
#'   \code{NULL} when supplying \code{filename} and \code{varname}.
#' @param filename A NetCDF file name to import if \code{data} is not provided.
#' @param varname The name of the variable to extract from \code{filename}.
#' @param candidates A character vector of distribution names to fit.
#' @param ignore_zeros A logical value, if \code{TRUE} zeros will be ignored.
#'   Default is \code{FALSE}.
#' @param zero_threshold The threshold below which values are considered zero.
#'   Default is 0.01.
#' @param parallel Logical, whether to use parallel processing.
#' @param ncores Number of cores to use for parallel computations. Default is 2.
#' @param shared_memory Logical, when parallel, share the grid with workers via a
#'   file-backed \code{big.matrix} (mmap, single machine) instead of serialising
#'   column chunks. Set \code{FALSE} for multi-node \code{plan(cluster)} setups.
#'   Default \code{TRUE}.
#' @param order Optional named list mapping a candidate name to the vector of
#'   L-moment orders matched by its optimiser, e.g.
#'   \code{list(gengamma = 1:5, expweibull = 1:3)}. Only the numerically-fitted
#'   distributions accept it; passed through to \code{\link{fitlm_nxts}}.
#'   Default \code{NULL}.
#' @param ... Additional arguments passed to \code{\link{nc2xts}} when
#'   \code{filename} and \code{varname} are provided.
#'
#' @return A list with components \code{fit_results} (a named list, one element
#'   per candidate, each containing \code{raster_params}, \code{raster_TheorLMom},
#'   \code{raster_DataLMom}, and \code{raster_GoF} \code{RasterBrick} objects)
#'   and \code{gof_plots} (a \code{ggplot} density comparison of the six GoF
#'   metrics across all grid cells, faceted by metric and coloured by candidate).
#'
#' @examples
#' \dontrun{
#' # Simulated 3-cell grid
#' n <- 365
#' dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = n)
#' vals <- cbind(cell1 = rgamma(n, shape = 0.8, scale = 3),
#'               cell2 = rgamma(n, shape = 1.2, scale = 2),
#'               cell3 = rgamma(n, shape = 0.6, scale = 4))
#' coords <- data.frame(lon = c(10, 11, 12), lat = c(45, 46, 47))
#' grid <- sxts(vals, order.by = dates, coords = coords,
#'              projection = "+proj=longlat +datum=WGS84")
#'
#' fits <- fitlm_nc(grid, candidates = c("exp", "gamma3"), ignore_zeros = TRUE)
#' names(fits$fit_results)
#' fits$gof_plots
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_density theme_bw facet_wrap vars
#' @importFrom lubridate parse_date_time
#'
#' @export
#'
fitlm_nc = function(data = NULL, filename = NA, varname = NA,
                    candidates = 'norm', ignore_zeros = FALSE, zero_threshold = 0.01,
                    parallel = FALSE, ncores = 2, shared_memory = TRUE, order = NULL, ...){

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
    ncdf_sxts <- xts::xts(x = tt, order.by = dates)
    proj_str  <- raster::projection(data)
    rm(t, tt)
  }

  keep     <- rep(TRUE, ncol(ncdf_sxts))
  fit_sxts <- ncdf_sxts
  if (ignore_zeros) {
    keep     <- colSums(coredata(ncdf_sxts) > zero_threshold, na.rm = TRUE) >= 1
    fit_sxts <- ncdf_sxts[, keep]          # smaller object -> RAM win in parallel workers
  }

  # Scatter the per-cell fits back onto the full grid (NA for dropped cells) so
  # rasterFromXYZ always receives the complete, regular lattice.
  scatter <- function(mat) {
    full <- matrix(NA_real_, nrow = length(keep), ncol = ncol(mat))
    full[keep, ] <- mat
    colnames(full) <- colnames(mat)
    full
  }

  ncdf_fits <- fitlm_nxts(fit_sxts, ignore_zeros = ignore_zeros,
                           candidates = candidates, zero_threshold = zero_threshold,
                           parallel = parallel, ncores = ncores,
                           shared_memory = shared_memory,
                           diagnostic_plots = FALSE, order = order)$params

  result_list <- list()
  gof_list    <- list()

  for (candidate in candidates) {
    ncdf_params <- lapply(ncdf_fits, function(x) { unlist(x$params[[candidate]]$Param) })
    ncdf_params <- t(do.call(cbind, ncdf_params))
    ncdf_params <- cbind(coords, scatter(ncdf_params))
    raster_params <- raster::rasterFromXYZ(ncdf_params)
    raster::projection(raster_params) <- proj_str

    ncdf_TheorLMom <- lapply(ncdf_fits, function(x) { unlist(x$params[[candidate]]$TheorLMom) })
    ncdf_TheorLMom <- t(do.call(cbind, ncdf_TheorLMom))
    ncdf_TheorLMom <- cbind(coords, scatter(ncdf_TheorLMom))
    raster_TheorLMom <- raster::rasterFromXYZ(ncdf_TheorLMom)
    raster::projection(raster_TheorLMom) <- proj_str

    ncdf_DataLMom <- lapply(ncdf_fits, function(x) { unlist(x$params[[candidate]]$DataLMom) })
    ncdf_DataLMom <- do.call(rbind, ncdf_DataLMom)
    ncdf_DataLMom <- cbind(coords, scatter(ncdf_DataLMom))
    raster_DataLMom <- raster::rasterFromXYZ(ncdf_DataLMom)
    raster::projection(raster_DataLMom) <- proj_str

    ncdf_GoF <- lapply(ncdf_fits, function(x) { unlist(x$params[[candidate]]$GoF) })
    ncdf_GoF <- t(do.call(cbind, ncdf_GoF))
    ncdf_GoF <- cbind(coords, scatter(ncdf_GoF))
    raster_GoF <- raster::rasterFromXYZ(ncdf_GoF)
    raster::projection(raster_GoF) <- proj_str

    ncdf_GoF <- as.data.frame(ncdf_GoF)
    ncdf_GoF$Distribution <- as.factor(candidate)
    gof_list[[candidate]] <- ncdf_GoF
    result_list[[candidate]] <- list(raster_params = raster_params, raster_TheorLMom = raster_TheorLMom,
                                     raster_DataLMom = raster_DataLMom, raster_GoF = raster_GoF)
  }

  gof_df   <- do.call(rbind, gof_list)
  gof_df   <- gof_df[, c("MLE", "CM", "KS", "MSEquant", "DiffOfMax", "MeanDiffOf10Max", "Distribution")]
  gof_long <- reshape2::melt(gof_df, id.vars = "Distribution")
  gof_plot <- ggplot(data = gof_long, aes(x = value)) +
    geom_density(alpha = 0.3, aes(fill = Distribution)) +
    theme_bw() + facet_wrap(vars(variable), scales = "free")

  list_out <- list(fit_results = result_list, gof_plots = gof_plot)
  return(list_out)
}
