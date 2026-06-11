library(testthat)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
make_plain_ts <- function() {
  set.seed(10)
  n     <- 365 * 3
  vals  <- rnorm(n, mean = 5, sd = 2)
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 3600
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "V1")), order.by = dates)
}

make_zeros_ts <- function() {
  set.seed(11)
  ts  <- make_plain_ts()
  idx <- sample(nrow(ts), size = round(nrow(ts) * 0.20))
  ts[idx, ] <- 0
  ts
}

make_na_ts <- function() {
  set.seed(12)
  ts  <- make_plain_ts()
  idx <- sample(nrow(ts), size = round(nrow(ts) * 0.10))
  ts[idx, ] <- NA
  ts
}

make_multi_ts <- function() {
  set.seed(13)
  n     <- 365 * 3
  vals  <- matrix(rnorm(n * 2, mean = 5, sd = 2), ncol = 2,
                  dimnames = list(NULL, c("A", "B")))
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 3600
  xts::xts(vals, order.by = dates)
}

# ---------------------------------------------------------------------------
# basic_stats — return structure
# ---------------------------------------------------------------------------
test_that("basic_stats() returns list with 'plot' and 'stats_table'", {
  ts  <- make_plain_ts()
  res <- basic_stats(ts)
  expect_type(res, "list")
  expect_named(res, c("plot", "stats_table"), ignore.order = TRUE)
})

test_that("basic_stats() stats_table has 34 rows when ignore_zeros=FALSE", {
  ts  <- make_plain_ts()
  res <- basic_stats(ts)
  expect_equal(nrow(res$stats_table), 34L)
})

test_that("basic_stats() stats_table contains all expected row names", {
  ts       <- make_plain_ts()
  res      <- basic_stats(ts)
  expected <- c("NumofData", "NumofMisData", "PercOfMissingData",
                "Min", "Max", "Mean", "Var", "StDev",
                "Variation", "Mom3", "Skewness", "Kurtosis",
                "Lmean", "LScale", "L3", "L4",
                "LVariation", "LSkewness", "LKurtosis",
                "Pdr", "Q5", "Q25", "Q50", "Q75", "Q95", "IQR",
                "MeanDAfterZero", "VarDAfterZero",
                "MeanDBeforeZero", "VarDBeforeZero",
                "MeanDAfterD", "VarDAfterD", "ProbDD", "ProbNDND")
  expect_true(all(expected %in% rownames(res$stats_table)))
})

test_that("basic_stats() NumofData matches nrow(ts)", {
  ts  <- make_plain_ts()
  res <- basic_stats(ts)
  expect_equal(res$stats_table["NumofData", "Value"], nrow(ts))
})

test_that("basic_stats() Mean is close to mean(ts)", {
  ts  <- make_plain_ts()
  res <- basic_stats(ts)
  expect_equal(res$stats_table["Mean", "Value"],
               round(mean(ts, na.rm = TRUE), 2),
               tolerance = 1e-2)
})

test_that("basic_stats() stats_table has 26 rows when ignore_zeros=TRUE", {
  ts  <- make_zeros_ts()
  res <- basic_stats(ts, ignore_zeros = TRUE)
  expect_equal(nrow(res$stats_table), 26L)
})

test_that("basic_stats() Pdr is positive when series contains zeros", {
  ts  <- make_zeros_ts()
  res <- basic_stats(ts)
  expect_gt(res$stats_table["Pdr", "Value"], 0)
})

test_that("basic_stats() NumofMisData is positive when series contains NAs", {
  ts  <- make_na_ts()
  res <- basic_stats(ts)
  expect_gt(res$stats_table["NumofMisData", "Value"], 0)
})

test_that("basic_stats() pstart/pend restricts the plotted period without error", {
  ts  <- make_plain_ts()
  expect_no_error(
    basic_stats(ts, pstart = "2016-01-01", pend = "2016-12-31")
  )
})

# ---------------------------------------------------------------------------
# lmom_stats — return structure
# ---------------------------------------------------------------------------
test_that("lmom_stats() returns a data.frame with 5 rows", {
  ts  <- make_plain_ts()
  res <- lmom_stats(ts)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5L)
})

test_that("lmom_stats() row names are the five L-moment ratios", {
  ts       <- make_plain_ts()
  res      <- lmom_stats(ts)
  expected <- c("L-Mean", "L-Scale", "L-Skew", "L-Kurtosis", "L-CV")
  expect_true(all(expected %in% rownames(res)))
})

test_that("lmom_stats() has one column per variable in ts", {
  ts  <- make_multi_ts()
  res <- lmom_stats(ts)
  expect_equal(ncol(res), ncol(ts))
  expect_equal(colnames(res), colnames(ts))
})

test_that("lmom_stats() L-Mean is close to the sample mean for Normal data", {
  ts  <- make_plain_ts()
  res <- lmom_stats(ts)
  expect_equal(res["L-Mean", 1], round(mean(ts, na.rm = TRUE), 2),
               tolerance = 0.5)
})

test_that("lmom_stats() runs without error when ignore_zeros=TRUE", {
  ts <- make_zeros_ts()
  expect_no_error(lmom_stats(ts, ignore_zeros = TRUE))
})

test_that("lmom_stats() values are all numeric and finite", {
  ts  <- make_plain_ts()
  res <- lmom_stats(ts)
  expect_true(all(is.numeric(unlist(res))))
  expect_true(all(is.finite(unlist(res))))
})
