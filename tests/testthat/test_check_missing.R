library(testthat)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
make_daily_ts <- function() {
  set.seed(20)
  n     <- 365 * 2
  vals  <- rnorm(n, 5, 2)
  dates <- as.POSIXct("2018-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  ts    <- xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "Var1")),
                    order.by = dates)
  # Inject ~15% NAs
  na_idx      <- sample(n, size = round(n * 0.15))
  ts[na_idx, ] <- NA
  ts
}

make_multi_ts <- function() {
  set.seed(21)
  n     <- 365 * 2
  vals  <- matrix(rnorm(n * 2, 5, 2), ncol = 2,
                  dimnames = list(NULL, c("A", "B")))
  dates <- as.POSIXct("2018-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  ts    <- xts::xts(vals, order.by = dates)
  na_idx       <- sample(n, size = round(n * 0.10))
  ts[na_idx, ] <- NA
  ts
}

# ---------------------------------------------------------------------------
# check_missing — overall percentage
# ---------------------------------------------------------------------------
test_that("check_missing() returns list with 'prct_missing' element", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = "months")
  expect_type(res, "list")
  expect_true("prct_missing" %in% names(res))
})

test_that("check_missing() prct_missing values are between 0 and 100", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = "months")
  pct <- res$prct_missing
  expect_true(all(pct >= 0 & pct <= 100))
})

test_that("check_missing() prct_missing is positive when NAs are present", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = "months")
  expect_gt(res$prct_missing[[1]], 0)
})

# ---------------------------------------------------------------------------
# check_missing — period list elements
# ---------------------------------------------------------------------------
test_that("check_missing() creates list_months element for 'months' period", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = "months")
  expect_true("list_months" %in% names(res))
})

test_that("check_missing() list_months$prct_missing is an xts object", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = "months")
  expect_true(inherits(res$list_months$prct_missing, "xts"))
})

test_that("check_missing() creates list_years element for 'years' period", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = "years")
  expect_true("list_years" %in% names(res))
})

test_that("check_missing() handles multiple periods simultaneously", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = c("months", "years"))
  expect_true("list_months" %in% names(res))
  expect_true("list_years" %in% names(res))
})

# ---------------------------------------------------------------------------
# check_missing — group_months
# ---------------------------------------------------------------------------
test_that("check_missing() with group_months=TRUE adds grouped_months matrix", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = "months", group_months = TRUE)
  gm  <- res$list_months$prct_missing$grouped_months
  expect_true(is.matrix(gm))
  expect_equal(nrow(gm), 12L)
})

test_that("check_missing() grouped_months values are between 0 and 100", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = "months", group_months = TRUE)
  gm  <- res$list_months$prct_missing$grouped_months
  expect_true(all(gm >= 0 & gm <= 100, na.rm = TRUE))
})

test_that("check_missing() grouped_months ncol matches ncol of input", {
  ts  <- make_multi_ts()
  res <- check_missing(ts, periods = "months", group_months = TRUE)
  gm  <- res$list_months$prct_missing$grouped_months
  expect_equal(ncol(gm), ncol(ts))
})

# ---------------------------------------------------------------------------
# check_missing — plot argument
# ---------------------------------------------------------------------------
test_that("check_missing() with plot=TRUE includes figure element", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = "months", plot = TRUE)
  expect_true("figure" %in% names(res$list_months))
  expect_true(inherits(res$list_months$figure, "gg"))
})

test_that("check_missing() with plot=FALSE omits figure element", {
  ts  <- make_daily_ts()
  res <- check_missing(ts, periods = "months", plot = FALSE)
  expect_false("figure" %in% names(res$list_months))
})

# ---------------------------------------------------------------------------
# check_missing — multi-column input
# ---------------------------------------------------------------------------
test_that("check_missing() works on two-column xts", {
  ts  <- make_multi_ts()
  expect_no_error(check_missing(ts, periods = "months"))
})

# ---------------------------------------------------------------------------
# plot_missing
# ---------------------------------------------------------------------------
test_that("plot_missing() returns a ggplot object", {
  ts  <- make_daily_ts()
  res <- plot_missing(ts)
  expect_true(inherits(res, "gg"))
})

test_that("plot_missing() runs without error when series contains NAs", {
  ts <- make_daily_ts()
  expect_no_error(plot_missing(ts))
})

test_that("plot_missing() runs without error when series has no NAs", {
  set.seed(99)
  n     <- 100
  vals  <- rnorm(n)
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  ts    <- xts::xts(matrix(vals, ncol = 1), order.by = dates)
  expect_no_error(plot_missing(ts))
})
