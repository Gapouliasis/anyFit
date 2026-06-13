library(testthat)

# ---------------------------------------------------------------------------
# Fixtures & helpers
# ---------------------------------------------------------------------------
make_period_ts <- function() {
  set.seed(20)
  n     <- 365 * 2
  vals  <- matrix(rnorm(n * 2, mean = 5, sd = 2), ncol = 2,
                  dimnames = list(NULL, c("A", "B")))
  vals[sample(n, 30), 1] <- NA                      # scattered NAs in A
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  xts::xts(vals, order.by = dates)
}

expect_parity <- function(actual, expected, tol = 1e-8) {
  a <- as.numeric(actual); e <- as.numeric(expected)
  fa <- is.finite(a); fe <- is.finite(e)
  expect_equal(fa, fe)
  if (any(fa)) expect_equal(a[fa], e[fe], tolerance = tol)
}

# Independent per-period, per-column oracle using base apply().
period_oracle <- function(ts, period, FUN) {
  ep <- xts::endpoints(ts, on = period)
  m  <- zoo::coredata(ts)
  t(vapply(seq_len(length(ep) - 1), function(p) {
    s <- m[(ep[p] + 1):ep[p + 1], , drop = FALSE]
    apply(s, 2, FUN)
  }, numeric(ncol(ts))))
}

STAT_NAMES <- c("NumofData", "NumofMisData", "PercOfMissingData", "Min", "Max",
                "Mean", "Var", "StDev", "Variation", "Mom3", "Skewness",
                "Kurtosis", "LMean", "LScale", "L3", "L4", "LVariation",
                "LSkewness", "LKurtosis", "Q5", "Q25", "Q50", "Q75", "Q95", "IQR")

# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------
test_that("period_stats returns one named xts per statistic", {
  ts  <- make_period_ts()
  res <- period_stats(ts, period = "months")

  expect_type(res, "list")
  expect_named(res, STAT_NAMES, ignore.order = TRUE)

  ep   <- xts::endpoints(ts, on = "months")
  nper <- length(ep) - 1
  for (nm in names(res)) {
    expect_true(xts::is.xts(res[[nm]]))
    expect_equal(ncol(res[[nm]]), ncol(ts))
    expect_equal(colnames(res[[nm]]), colnames(ts))
    expect_equal(nrow(res[[nm]]), nper)
  }
})

# ---------------------------------------------------------------------------
# Numerical parity against independent per-period oracles
# ---------------------------------------------------------------------------
test_that("period_stats Mean / Min / Max match an apply() oracle (NA-robust)", {
  ts  <- make_period_ts()
  res <- period_stats(ts, period = "months")

  expect_parity(zoo::coredata(res$Mean),
                period_oracle(ts, "months", function(c) mean(c, na.rm = TRUE)))
  expect_parity(zoo::coredata(res$Min),
                period_oracle(ts, "months", function(c) suppressWarnings(min(c, na.rm = TRUE))))
  expect_parity(zoo::coredata(res$Max),
                period_oracle(ts, "months", function(c) suppressWarnings(max(c, na.rm = TRUE))))
  expect_parity(zoo::coredata(res$StDev),
                period_oracle(ts, "months", function(c) sd(c, na.rm = TRUE)))
})

test_that("period_stats Q95 matches stats::quantile(type = 7) per period", {
  ts  <- make_period_ts()
  res <- period_stats(ts, period = "months")
  expect_parity(
    zoo::coredata(res$Q95),
    period_oracle(ts, "months",
                  function(c) stats::quantile(c, 0.95, na.rm = TRUE, names = FALSE, type = 7)))
})

test_that("period_stats Skewness/Kurtosis match the moments package per period", {
  ts  <- make_period_ts()
  res <- period_stats(ts, period = "months")
  expect_parity(zoo::coredata(res$Skewness),
                period_oracle(ts, "months", function(c) moments::skewness(c, na.rm = TRUE)))
  expect_parity(zoo::coredata(res$Kurtosis),
                period_oracle(ts, "months", function(c) moments::kurtosis(c, na.rm = TRUE)))
})

test_that("period_stats PercOfMissingData is consistent with the NA counts", {
  ts  <- make_period_ts()
  res <- period_stats(ts, period = "months")
  expect_parity(zoo::coredata(res$PercOfMissingData),
                zoo::coredata(res$NumofMisData) / zoo::coredata(res$NumofData) * 100)
})

# ---------------------------------------------------------------------------
# Edge cases & wide-equivalence
# ---------------------------------------------------------------------------
test_that("period_stats survives a fully-NA period without error", {
  ts <- make_period_ts()
  ts["2015-03", 1] <- NA                            # wipe one whole month in A
  expect_no_error(res <- period_stats(ts, period = "months"))
  march <- which(format(zoo::index(res$Mean), "%Y-%m") == "2015-03")
  expect_false(is.finite(as.numeric(res$Mean[march, 1])))
})

test_that("period_stats column k equals the single-series call on column k", {
  ts <- make_period_ts()
  wide <- period_stats(ts, period = "months")
  for (j in seq_len(ncol(ts))) {
    single <- period_stats(ts[, j], period = "months")
    for (nm in c("Mean", "StDev", "Q50", "Skewness", "LScale")) {
      expect_parity(zoo::coredata(wide[[nm]])[, j], zoo::coredata(single[[nm]]))
    }
  }
})
