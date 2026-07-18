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
#' The `*_LSpace` objects delimit the feasible L-moment ratio support region of
#' a distribution and back [LRatio_check()]. The `*_InitValues` tables map
#' sampled L-moment ratios to good starting parameters for the numerical
#' (L-BFGS-B) fit of the distributions that lack closed-form L-moments
#' (Burr XII, Dagum, Exponentiated Weibull, Generalised Gamma), as used by
#' [fitlm_nc()] and its relatives.
#'
#' @format
#' \itemize{
#'   \item `LMom_LSpace`: a data frame (17994 x 3) with columns `L-CV`,
#'     `L-Skewness`, `Dist`.
#'   \item `BurrXII_LSpace`, `Dagum_LSpace`: [sp::SpatialPolygonsDataFrame]
#'     objects delimiting each distribution's L-ratio support polygon.
#'   \item `Burr_InitValues`, `Dagum_InitValues`, `ExpWeibull_InitValues`,
#'     `GG_InitValues`: data frames with 8 columns (`lambda_1`, `lambda_2`,
#'     `tau_3`, `tau_4`, `lcv`, `scale`, `shape1`, `shape2`).
#' }
#' @keywords internal
#' @name lmoment-tables
#' @aliases LMom_LSpace BurrXII_LSpace Dagum_LSpace Burr_InitValues Dagum_InitValues ExpWeibull_InitValues GG_InitValues
NULL
