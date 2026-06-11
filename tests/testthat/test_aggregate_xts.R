library(testthat)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
make_hourly_ts <- function() {
  set.seed(30)
  n     <- 24 * 365 * 3          # 3 years of hourly data
  vals  <- abs(rnorm(n, 5, 2))
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 3600
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "P")),
           order.by = dates)
}

make_na_ts <- function() {
  set.seed(31)
  ts  <- make_hourly_ts()
  idx <- sample(nrow(ts), size = round(nrow(ts) * 0.05))
  ts[idx, ] <- NA
  ts
}

# ---------------------------------------------------------------------------
# aggregate_xts — basic return structure
# ---------------------------------------------------------------------------
test_that("aggregate_xts() returns list with Combined_Plot and list_days", {
  ts  <- make_hourly_ts()
  res <- aggregate_xts(ts, periods = "days")
  expect_type(res, "list")
  expect_true("Combined_Plot" %in% names(res))
  expect_true("list_days" %in% names(res))
})

test_that("aggregate_xts() aggregated daily xts has fewer rows than hourly input", {
  ts   <- make_hourly_ts()
  res  <- aggregate_xts(ts, periods = "days")
  agg  <- res$list_days$aggregated
  expect_true(inherits(agg, "xts"))
  expect_lt(nrow(agg), nrow(ts))
})

test_that("aggregate_xts() creates one list element per requested period", {
  ts  <- make_hourly_ts()
  res <- aggregate_xts(ts, periods = c("days", "months", "years"))
  expect_true("list_days"   %in% names(res))
  expect_true("list_months" %in% names(res))
  expect_true("list_years"  %in% names(res))
})

test_that("aggregate_xts() each period element has 'aggregated' and 'figure'", {
  ts  <- make_hourly_ts()
  res <- aggregate_xts(ts, periods = "months")
  expect_true("aggregated" %in% names(res$list_months))
  expect_true("figure"     %in% names(res$list_months))
})

# ---------------------------------------------------------------------------
# aggregate_xts — FUN argument
# ---------------------------------------------------------------------------
test_that("aggregate_xts() with FUN='sum' gives sums close to manual calculation", {
  set.seed(32)
  n     <- 24 * 7   # one week of hourly data
  vals  <- rep(1, n)
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, n - 1) * 3600
  ts    <- xts::xts(matrix(vals, ncol = 1), order.by = dates)
  res   <- aggregate_xts(ts, periods = "days", FUN = "sum")
  agg   <- res$list_days$aggregated
  # Each day should sum to 24 (24 × 1)
  expect_true(all(abs(as.numeric(agg) - 24) < 1e-6))
})

# ---------------------------------------------------------------------------
# aggregate_xts — period_multiplier error
# ---------------------------------------------------------------------------
test_that("aggregate_xts() errors when period_multiplier length mismatches periods", {
  ts <- make_hourly_ts()
  expect_error(
    aggregate_xts(ts, periods = c("days", "months"), period_multiplier = c(2)),
    "same length as periods"
  )
})

test_that("aggregate_xts() accepts period_multiplier with same length as periods", {
  ts <- make_hourly_ts()
  expect_no_error(
    aggregate_xts(ts, periods = c("days", "months"), period_multiplier = c(1, 1))
  )
})

# ---------------------------------------------------------------------------
# aggregate_xts — NAs in input
# ---------------------------------------------------------------------------
test_that("aggregate_xts() runs without error when input contains NAs", {
  ts <- make_na_ts()
  expect_no_error(aggregate_xts(ts, periods = "months", FUN = "mean"))
})

# ---------------------------------------------------------------------------
# period_stats — return type and columns
# ---------------------------------------------------------------------------
test_that("period_stats() returns an xts object", {
  ts  <- make_hourly_ts()
  res <- period_stats(ts, period = "months")
  expect_true(inherits(res, "xts"))
})

test_that("period_stats() contains core statistic columns", {
  ts       <- make_hourly_ts()
  res      <- period_stats(ts, period = "months")
  required <- c("NumofData", "PercOfMissingData", "Mean", "Var", "Q50")
  expect_true(all(required %in% colnames(res)))
})

test_that("period_stats() monthly aggregation has correct row count", {
  ts  <- make_hourly_ts()
  res <- period_stats(ts, period = "months")
  # 3 years × 12 months
  expect_equal(nrow(res), 36L)
})

test_that("period_stats() yearly aggregation has correct row count", {
  ts  <- make_hourly_ts()
  res <- period_stats(ts, period = "years")
  expect_equal(nrow(res), 3L)
})

test_that("period_stats() PercOfMissingData is positive when NAs are present", {
  ts  <- make_na_ts()
  res <- period_stats(ts, period = "years")
  expect_gt(sum(res[, "PercOfMissingData"]), 0)
})

test_that("period_stats() statistic values are all numeric", {
  ts  <- make_hourly_ts()
  res <- period_stats(ts, period = "months")
  expect_true(all(sapply(as.data.frame(res), is.numeric)))
})

# ---------------------------------------------------------------------------
# period_apply_nc — tested on sxts
# ---------------------------------------------------------------------------
test_that("period_apply_nc() on sxts returns sxts with fewer rows than input", {
  set.seed(33)
  n_time <- 24 * 365          # 1 year hourly
  n_pts  <- 4
  mat    <- matrix(runif(n_time * n_pts, 0, 10), nrow = n_time, ncol = n_pts)
  dates  <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, n_time - 1) * 3600
  coords <- data.frame(x = c(1, 2, 1, 2), y = c(10, 10, 11, 11))
  obj    <- sxts(mat, order.by = dates, coords = coords,
                 projection = "+proj=longlat +datum=WGS84")
  res    <- period_apply_nc(data = obj, period = "months", FUN = "mean")
  expect_true(inherits(res, "sxts"))
  expect_lt(nrow(res), nrow(obj))
})
