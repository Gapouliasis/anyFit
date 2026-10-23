#' World country and continent boundary polygons
#'
#' A global map of country and territory boundaries whose attributes back
#' anyFit's spatial helpers. The `name` and `continent` columns are the ones
#' matched when a `country` or `continent` is supplied to [mask.sxts()] and
#' [zonal_stats()], so the exact spelling in these columns is what those
#' arguments expect (e.g. the United Kingdom is stored as
#' `"U.K. of Great Britain and Northern Ireland"`).
#'
#' @format A [sp::SpatialPolygonsDataFrame] with 256 features and 8 attribute
#'   columns:
#' \describe{
#'   \item{iso3}{Three-letter ISO 3166-1 alpha-3 country code.}
#'   \item{status}{Sovereignty or status of the territory.}
#'   \item{color_code}{Colour group used for cartographic display.}
#'   \item{name}{Country or territory name (matched by `country =` arguments).}
#'   \item{continent}{Continent name (matched by `continent =` arguments).}
#'   \item{region}{Sub-regional grouping.}
#'   \item{iso_3166_1_}{ISO 3166-1 numeric code.}
#'   \item{french_shor}{Short country name in French.}
#' }
#' @source Bundled global country/territory boundary polygons.
#' @seealso [mask.sxts()], [zonal_stats()]
#' @examples
#' # Continents available for `continent =` masking
#' sort(unique(world_data$continent))
#'
#' # Exact country name used by `country =` arguments
#' grep("Great Britain", world_data$name, value = TRUE)
"world_data"


#' Precomputed L-moment lookup tables (internal)
#'
#' Reference tables consumed internally by the distribution-fitting routines.
#' The `*_LSpace` objects delimit the attainable (L-CV, L-skewness) region of
#' a distribution over the shape bounds of its fitter. They back [LRatio_check()]
#' and decide the `"outside_domain"` status of [fitlm_burr()], [fitlm_dagum()],
#' [fitlm_gengamma()] and [fitlm_expweibull()]. They are built by
#' `data-raw/make_LSpace.R`, which maps a 120 x 120 log grid over the shape box
#' through the exact L-moment engine and unions the grid-cell images. The
#' `*_InitValues` tables hold unit-scale shape pairs with their L-moments. They
#' seed the numerical (L-BFGS-B) fit of these four families by nearest-neighbour
#' matching on the scale-free (L-CV, \eqn{\tau_3}, \eqn{\tau_4}), and cover the
#' fitter boxes (`data-raw/make_*_InitValues.R`).
#'
#' @format
#' \itemize{
#'   \item `LMom_LSpace`: a data frame (17994 x 3) with columns `L-CV`,
#'     `L-Skewness`, `Dist`.
#'   \item `BurrXII_LSpace`, `Dagum_LSpace`: [sp::SpatialPolygonsDataFrame]
#'     objects delimiting each distribution's attainable L-ratio region (the
#'     generalised gamma and exponentiated Weibull polygons are internal).
#'   \item `Burr_InitValues`, `Dagum_InitValues`, `ExpWeibull_InitValues`,
#'     `GG_InitValues`: data frames with 8 columns (`lambda_1`, `lambda_2`,
#'     `tau_3`, `tau_4`, `lcv`, `scale`, `shape1`, `shape2`), at `scale = 1`.
#' }
#' @keywords internal
#' @name lmoment-tables
#' @aliases LMom_LSpace BurrXII_LSpace Dagum_LSpace Burr_InitValues Dagum_InitValues ExpWeibull_InitValues GG_InitValues
NULL
