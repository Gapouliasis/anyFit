library(testthat)

# ---------------------------------------------------------------------------
# Shared fixtures — 5-year daily xts (single and two-column) with some zeros
# ---------------------------------------------------------------------------
make_daily_ts <- function() {
  set.seed(40)
  n     <- 365 * 5
  vals  <- abs(rnorm(n, 5, 2))
  # ~10% zeros
  zeros <- sample(n, size = round(n * 0.10))
  vals[zeros] <- 0
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "P")),
           order.by = dates)
}

make_multi_ts <- function() {
  set.seed(41)
  n     <- 365 * 5
  vals  <- matrix(abs(rnorm(n * 2, 5, 2)), ncol = 2,
                  dimnames = list(NULL, c("P1", "P2")))
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  xts::xts(vals, order.by = dates)
}

# ---------------------------------------------------------------------------
# monthly_stats — base scale (aggregated=FALSE)
# ---------------------------------------------------------------------------
test_that("monthly_stats() returns list with 'base_stats' and 'fbase'", {
  ts  <- make_daily_ts()
  res <- monthly_stats(ts)
  expect_type(res, "list")
  expect_named(res, c("base_stats", "fbase"), ignore.order = TRUE)
})

test_that("monthly_stats() base_stats has 12 columns (one per month)", {
  ts  <- make_daily_ts()
  res <- monthly_stats(ts)
  expect_equal(ncol(res$base_stats), 12L)
})

test_that("monthly_stats() base_stats column names are month names", {
  ts  <- make_daily_ts()
  res <- monthly_stats(ts)
  expect_equal(colnames(res$base_stats), month.name)
})

test_that("monthly_stats() fbase is a ggplot/patchwork object", {
  ts  <- make_daily_ts()
  res <- monthly_stats(ts)
  expect_true(inherits(res$fbase, "gg") || inherits(res$fbase, "patchwork"))
})

test_that("monthly_stats() runs without error when ignore_zeros=TRUE", {
  ts <- make_daily_ts()
  expect_no_error(monthly_stats(ts, ignore_zeros = TRUE))
})

# ---------------------------------------------------------------------------
# monthly_stats — aggregated scale (aggregated=TRUE)
# ---------------------------------------------------------------------------
test_that("monthly_stats() with aggregated=TRUE returns agg_stats, faggre, lag1", {
  ts  <- make_daily_ts()
  res <- monthly_stats(ts, aggregated = TRUE)
  expect_named(res, c("agg_stats", "faggre", "lag1"), ignore.order = TRUE)
})

test_that("monthly_stats() agg_stats has 12 columns for aggregated=TRUE", {
  ts  <- make_daily_ts()
  res <- monthly_stats(ts, aggregated = TRUE)
  expect_equal(ncol(res$agg_stats), 12L)
})

test_that("monthly_stats() lag1 is a data.frame with 12 rows", {
  ts  <- make_daily_ts()
  res <- monthly_stats(ts, aggregated = TRUE)
  expect_s3_class(res$lag1, "data.frame")
  expect_equal(nrow(res$lag1), 12L)
})

# ---------------------------------------------------------------------------
# monthly_boxplots
# ---------------------------------------------------------------------------
test_that("monthly_boxplots() returns a ggplot object for single-column ts", {
  ts  <- make_daily_ts()
  res <- monthly_boxplots(ts)
  expect_true(inherits(res, "gg"))
})

test_that("monthly_boxplots() returns a ggplot object for two-column ts", {
  ts  <- make_multi_ts()
  res <- monthly_boxplots(ts)
  expect_true(inherits(res, "gg"))
})

test_that("monthly_boxplots() runs without error when ignore_zeros=TRUE", {
  ts <- make_daily_ts()
  expect_no_error(monthly_boxplots(ts, ignore_zeros = TRUE))
})

# ---------------------------------------------------------------------------
# monthly_ecdf
# ---------------------------------------------------------------------------
test_that("monthly_ecdf() returns a ggplot or patchwork object", {
  ts  <- make_daily_ts()
  res <- monthly_ecdf(ts)
  expect_true(inherits(res, "gg") || inherits(res, "patchwork"))
})

test_that("monthly_ecdf() runs without error when ignore_zeros=TRUE", {
  ts <- make_daily_ts()
  expect_no_error(monthly_ecdf(ts, ignore_zeros = TRUE))
})

# ---------------------------------------------------------------------------
# monthly_violins
# ---------------------------------------------------------------------------
test_that("monthly_violins() returns a ggplot object for single-column ts", {
  ts  <- make_daily_ts()
  res <- monthly_violins(ts)
  expect_true(inherits(res, "gg"))
})

test_that("monthly_violins() returns a ggplot object for two-column ts", {
  ts  <- make_multi_ts()
  res <- monthly_violins(ts)
  expect_true(inherits(res, "gg"))
})

test_that("monthly_violins() runs without error when ignore_zeros=TRUE", {
  ts <- make_daily_ts()
  expect_no_error(monthly_violins(ts, ignore_zeros = TRUE))
})

# ---------------------------------------------------------------------------
# ridge_plots
# ---------------------------------------------------------------------------
test_that("ridge_plots() returns list with 'plot_all' and 'plot_monthly'", {
  ts  <- make_daily_ts()
  res <- ridge_plots(ts)
  expect_type(res, "list")
  expect_named(res, c("plot_all", "plot_monthly"), ignore.order = TRUE)
})

test_that("ridge_plots() plot_all is a ggplot object", {
  ts  <- make_daily_ts()
  res <- ridge_plots(ts)
  expect_true(inherits(res$plot_all, "gg"))
})

test_that("ridge_plots() plot_monthly is a list with one entry per column", {
  ts  <- make_multi_ts()
  res <- ridge_plots(ts)
  expect_type(res$plot_monthly, "list")
  expect_equal(length(res$plot_monthly), ncol(ts))
})

test_that("ridge_plots() each plot_monthly element is a ggplot", {
  ts  <- make_daily_ts()
  res <- ridge_plots(ts)
  expect_true(all(sapply(res$plot_monthly, inherits, "gg")))
})

test_that("ridge_plots() runs without error when ignore_zeros=TRUE", {
  ts <- make_daily_ts()
  expect_no_error(ridge_plots(ts, ignore_zeros = TRUE))
})
