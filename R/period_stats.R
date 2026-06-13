#' @title period_stats
#'
#' @description Function to calculate period based basic statistics. Specifically, the number of data points,
#' number of missing data, percentage of missing data, min, max, mean, variance,
#' standard deviation, variation, 3rd moment, skewness, kurtosis, l-statistics, i.e. mean, scale,
#' 3rd and 4th order l-moments, l-moment variation, l-moment skewness, l-moment kurtosis,
#' quantiles of 5, 25, 50, 75 and 95, and the inter quartile range.
#'
#' Accepts a multi-column (wide) xts: every statistic is computed column-wise for
#' each period using single-pass \code{matrixStats} reducers (NA-robust, i.e.
#' \code{na.rm = TRUE}). The result is one xts per statistic.
#'
#' @param ts a xts object containing the time series data
#' @param period a period for which to calculate statistics. Default is 'months'
#' @param period_multiplier a multiplier to define arbitrary periods based on the period argument.
#' E.g. 2 months, 3 years etc. Default is 1.
#'
#' @return A named list with one element per statistic (\code{NumofData},
#' \code{NumofMisData}, \code{PercOfMissingData}, \code{Min}, \code{Max},
#' \code{Mean}, \code{Var}, \code{StDev}, \code{Variation}, \code{Mom3},
#' \code{Skewness}, \code{Kurtosis}, \code{LMean}, \code{LScale}, \code{L3},
#' \code{L4}, \code{LVariation}, \code{LSkewness}, \code{LKurtosis}, \code{Q5},
#' \code{Q25}, \code{Q50}, \code{Q75}, \code{Q95}, \code{IQR}). Each element is an
#' xts with one row per period and one column per input series (columns named
#' after the input).
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' # Wide input: statistics computed per column, per period
#' pstats <- period_stats(data[, 1:4], period = "months", period_multiplier = 3)
#' pstats$Mean   # xts of period means, one column per series
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
