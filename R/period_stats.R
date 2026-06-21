#' Period-based summary statistics
#'
#' @description
#' Computes a comprehensive suite of period-based summary statistics for
#' one or more time series. Statistics include count, missing-data
#' diagnostics, min, max, mean, variance, standard deviation, coefficient
#' of variation, third moment, skewness, kurtosis, L-moments (mean, scale,
#' L3, L4, L-CV, L-skewness, L-kurtosis), quantiles (5, 25, 50, 75, 95),
#' and inter-quartile range. Each statistic is computed column-wise per
#' period using \code{matrixStats} with
#' \code{na.rm = TRUE}; L-moments are computed independently via
#' \code{lmom::samlmu} since no column-wise matrixStats equivalent
#' exists. The result is a named list of xts objects, one per statistic,
#' each with one row per period and one column per input series. Accepts
#' any period supported by \code{xts::endpoints} (e.g. \code{"months"},
#' \code{"years"}) with an optional multiplier for custom period lengths.
#'
#' @param ts An xts object containing the time series data.
#' @param period A period string passed to \code{endpoints} (default \code{"months"}).
#' @param period_multiplier Integer multiplier for custom period lengths (default 1).
#'
#' @return A named list with one element per statistic (\code{NumofData},
#'   \code{NumofMisData}, \code{PercOfMissingData}, \code{Min}, \code{Max},
#'   \code{Mean}, \code{Var}, \code{StDev}, \code{Variation}, \code{Mom3},
#'   \code{Skewness}, \code{Kurtosis}, \code{LMean}, \code{LScale}, \code{L3},
#'   \code{L4}, \code{LVariation}, \code{LSkewness}, \code{LKurtosis}, \code{Q5},
#'   \code{Q25}, \code{Q50}, \code{Q75}, \code{Q95}, \code{IQR}). Each element is an
#'   xts with one row per period and one column per input series.
#'
#' @examples
#' # Synthetic xts
#' set.seed(123)
#' dates <- seq(as.POSIXct("2000-01-01"), as.POSIXct("2000-12-31"), by = "day")
#' vals <- matrix(rnorm(length(dates) * 3, 10, 5), nrow = length(dates), ncol = 3)
#' colnames(vals) <- c("A", "B", "C")
#' ts <- xts::xts(vals, order.by = dates)
#'
#' # Statistics per 3-month periods
#' pstats <- period_stats(ts, period = "months", period_multiplier = 3)
#' pstats$Mean
#' pstats$Q95
#'
#' @importFrom xts endpoints period.apply xts
#' @importFrom zoo coredata index
#' @importFrom matrixStats colMins colMaxs colMeans2 colVars colSds colSums2 colQuantiles
#'
#' @export
#'

period_stats <-function(ts, period = 'months', period_multiplier = 1){
  spec_period <- endpoints(ts, on = period, k = period_multiplier)
  cn <- colnames(ts)

  # Apply a column-wise reducer to every period, returning a named xts
  # [periods x series]. matrixStats reducers drop names, so reattach them.
  reduce <- function(reducer){
    out <- xts::period.apply(ts, INDEX = spec_period, FUN = reducer)
    colnames(out) <- cn
    out
  }

  NumofData         <- reduce(function(s) rep(nrow(s), ncol(s)))
  NumofMisData      <- reduce(function(s) matrixStats::colSums2(is.na(coredata(s))))
  PercOfMissingData <- (NumofMisData / NumofData) * 100

  Min   <- reduce(function(s) suppressWarnings(matrixStats::colMins(coredata(s), na.rm = TRUE)))
  Max   <- reduce(function(s) suppressWarnings(matrixStats::colMaxs(coredata(s), na.rm = TRUE)))
  Mean  <- reduce(function(s) matrixStats::colMeans2(coredata(s), na.rm = TRUE))
  Var   <- reduce(function(s) matrixStats::colVars(coredata(s), na.rm = TRUE))
  StDev <- reduce(function(s) matrixStats::colSds(coredata(s), na.rm = TRUE))
  Variation <- StDev / Mean

  Skewness <- reduce(function(s) .col_moment(s, 'skewness'))
  Kurtosis <- reduce(function(s) .col_moment(s, 'kurtosis'))
  Mom3     <- reduce(function(s) .col_moment(s, 'mom3'))

  Q5  <- reduce(function(s) matrixStats::colQuantiles(coredata(s), probs = 0.05, na.rm = TRUE))
  Q25 <- reduce(function(s) matrixStats::colQuantiles(coredata(s), probs = 0.25, na.rm = TRUE))
  Q50 <- reduce(function(s) matrixStats::colQuantiles(coredata(s), probs = 0.50, na.rm = TRUE))
  Q75 <- reduce(function(s) matrixStats::colQuantiles(coredata(s), probs = 0.75, na.rm = TRUE))
  Q95 <- reduce(function(s) matrixStats::colQuantiles(coredata(s), probs = 0.95, na.rm = TRUE))
  IQR <- abs(Q75 - Q25)

  # L-moments have no matrixStats equivalent: per-period, per-column samlmu.
  L <- .period_lmoments(ts, spec_period, order.by = index(Mean))

  list(NumofData = NumofData, NumofMisData = NumofMisData,
       PercOfMissingData = PercOfMissingData,
       Min = Min, Max = Max, Mean = Mean, Var = Var, StDev = StDev,
       Variation = Variation, Mom3 = Mom3, Skewness = Skewness, Kurtosis = Kurtosis,
       LMean = L$LMean, LScale = L$LScale, L3 = L$L3, L4 = L$L4,
       LVariation = L$LVariation, LSkewness = L$LSkewness, LKurtosis = L$LKurtosis,
       Q5 = Q5, Q25 = Q25, Q50 = Q50, Q75 = Q75, Q95 = Q95, IQR = IQR)
}

#' Column-wise moment statistics
#'
#' @description
#' Computes column-wise skewness, kurtosis, or third central moment (Mom3)
#' for a period slice, NA-robust. Switches on the \code{which} argument:
#' \code{"skewness"} returns the population skewness,
#' \code{"kurtosis"} the population (non-excess) kurtosis, and any other
#' value returns Mom3 = \eqn{\sum(x-\mu)^3 / (n-1)}{sum((x-m)^3)/(n-1)}.
#'
#' @noRd
# Column-wise central-moment statistic for a period slice (NA-robust). Reproduces
# moments::skewness / moments::kurtosis (population / non-excess) and the
# original Mom3 = sum((x-mean)^3)/(n-1).
.col_moment <- function(s, which){
  s  <- coredata(s)
  mu <- matrixStats::colMeans2(s, na.rm = TRUE)
  xc <- sweep(s, 2, mu, FUN = '-')
  m2 <- matrixStats::colMeans2(xc^2, na.rm = TRUE)
  if (which == 'skewness') return(matrixStats::colMeans2(xc^3, na.rm = TRUE) / m2^1.5)
  if (which == 'kurtosis') return(matrixStats::colMeans2(xc^4, na.rm = TRUE) / m2^2)
  nobs <- matrixStats::colSums2(!is.na(s))
  matrixStats::colSums2(xc^3, na.rm = TRUE) / (nobs - 1)
}

#' Per-period L-moments
#'
#' @description
#' Computes per-period, per-column L-moments via \code{lmom::samlmu}, returning
#' one xts per L-statistic aligned to the period timestamps in
#' \code{order.by}. Handles both single-column and multi-column xts inputs.
#' Re-uses \code{.basic_stats_lmom} from \code{basic_stats.R}.
#'
#' @noRd
# Per-period, per-column L-moments, returned as one xts per L-statistic aligned
# to the period timestamps in order.by. Reuses .basic_stats_lmom (basic_stats.R).
.period_lmoments <- function(ts, ep, order.by){
  m <- coredata(ts)
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)
  k <- ncol(m)
  nper <- length(ep) - 1
  labs <- c('LMean', 'LScale', 'L3', 'L4', 'LVariation', 'LSkewness', 'LKurtosis')
  mats <- stats::setNames(lapply(labs, function(x) matrix(NA_real_, nper, k)), labs)
  for (p in seq_len(nper)){
    rows <- (ep[p] + 1):ep[p + 1]
    for (j in seq_len(k)){
      v <- .basic_stats_lmom(m[rows, j])
      for (s in seq_along(labs)) mats[[s]][p, j] <- v[s]
    }
  }
  cn <- colnames(ts)
  lapply(mats, function(M){ colnames(M) <- cn; xts::xts(M, order.by = order.by) })
}
