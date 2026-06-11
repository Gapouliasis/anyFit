library(testthat)

# ---------------------------------------------------------------------------
# Shared fixture — 5-year daily AR(1) time series
# ---------------------------------------------------------------------------
make_ar1_ts <- function() {
  set.seed(1)
  n    <- 365 * 5
  vals <- numeric(n)
  vals[1] <- rnorm(1, 5, 1)
  for (i in 2:n) vals[i] <- 0.6 * vals[i - 1] + rnorm(1, 2, 1)
  dates <- as.POSIXct("2010-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  xts::xts(matrix(vals, ncol = 1), order.by = dates)
}

make_zeros_ts <- function() {
  set.seed(2)
  ts  <- make_ar1_ts()
  idx <- sample(nrow(ts), size = 200)
  ts[idx, ] <- 0
  ts
}

# ---------------------------------------------------------------------------
# CAS_ACF
# ---------------------------------------------------------------------------
test_that("CAS_ACF() returns data.frame with lag and ACF columns", {
  res <- CAS_ACF(kappa = 1, beta = 0.5)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("lag", "ACF") %in% names(res)))
})

test_that("CAS_ACF() lag 0 value equals 1", {
  res <- CAS_ACF(kappa = 1, beta = 0.5)
  expect_equal(res$ACF[res$lag == 0], 1)
})

test_that("CAS_ACF() respects lag_max and returns lag_max+1 rows", {
  res <- CAS_ACF(kappa = 1, beta = 0.5, lag_max = 5)
  expect_equal(nrow(res), 6L)
  expect_equal(max(res$lag), 5L)
})

test_that("CAS_ACF() produces monotonically decreasing ACF values", {
  res <- CAS_ACF(kappa = 1, beta = 0.5, lag_max = 10)
  expect_true(all(diff(res$ACF) < 0))
})

# ---------------------------------------------------------------------------
# HK_ACF
# ---------------------------------------------------------------------------
test_that("HK_ACF() returns data.frame with lag and ACF columns", {
  res <- HK_ACF(H = 0.7)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("lag", "ACF") %in% names(res)))
})

test_that("HK_ACF() lag 0 value equals 1", {
  res <- HK_ACF(H = 0.7)
  expect_equal(res$ACF[res$lag == 0], 1)
})

test_that("HK_ACF() respects lag_max and returns lag_max+1 rows", {
  res <- HK_ACF(H = 0.7, lag_max = 7)
  expect_equal(nrow(res), 8L)
})

test_that("HK_ACF() with H=0.5 gives SRD-like rapid decay", {
  res_hk  <- HK_ACF(H = 0.5, lag_max = 5)
  res_hk2 <- HK_ACF(H = 0.9, lag_max = 5)
  # H=0.5 decays faster than H=0.9 (LRD)
  expect_lt(res_hk$ACF[res_hk$lag == 5], res_hk2$ACF[res_hk2$lag == 5])
})

# ---------------------------------------------------------------------------
# SRD_ACF
# ---------------------------------------------------------------------------
test_that("SRD_ACF() returns data.frame with lag and ACF columns", {
  res <- SRD_ACF(kappa = 0.3)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("lag", "ACF") %in% names(res)))
})

test_that("SRD_ACF() lag 0 value equals 1", {
  res <- SRD_ACF(kappa = 0.3)
  expect_equal(res$ACF[res$lag == 0], 1)
})

test_that("SRD_ACF() respects lag_max and returns lag_max+1 rows", {
  res <- SRD_ACF(kappa = 0.5, lag_max = 8)
  expect_equal(nrow(res), 9L)
})

test_that("SRD_ACF() larger kappa gives faster decay", {
  slow <- SRD_ACF(kappa = 0.1, lag_max = 5)
  fast <- SRD_ACF(kappa = 1.0, lag_max = 5)
  expect_lt(fast$ACF[fast$lag == 5], slow$ACF[slow$lag == 5])
})

# ---------------------------------------------------------------------------
# fit_ACF
# ---------------------------------------------------------------------------
test_that("fit_ACF() returns list with ACF_params, ACF_fitted, ACF_plot", {
  ts  <- make_ar1_ts()
  res <- fit_ACF(ts, lag_max = 5)
  expect_type(res, "list")
  expect_named(res, c("ACF_params", "ACF_fitted", "ACF_plot"), ignore.order = TRUE)
})

test_that("fit_ACF() ACF_params contains an entry per fitted type", {
  ts  <- make_ar1_ts()
  res <- fit_ACF(ts, lag_max = 5, type = list("CAS", "HK", "SRD"))
  expect_equal(length(res$ACF_params), 3L)
  expect_named(res$ACF_params, c("CAS", "HK", "SRD"), ignore.order = TRUE)
})

test_that("fit_ACF() fitted parameters are numeric and finite", {
  ts  <- make_ar1_ts()
  res <- fit_ACF(ts, lag_max = 5)
  all_params <- unlist(res$ACF_params)
  expect_true(all(is.numeric(all_params)))
  expect_true(all(is.finite(all_params)))
})

test_that("fit_ACF() with type=list('SRD') returns only SRD in ACF_params", {
  ts  <- make_ar1_ts()
  res <- fit_ACF(ts, lag_max = 5, type = list("SRD"))
  expect_equal(length(res$ACF_params), 1L)
  expect_named(res$ACF_params, "SRD")
})

test_that("fit_ACF() ACF_fitted has lag column and one column per fitted type", {
  ts  <- make_ar1_ts()
  res <- fit_ACF(ts, lag_max = 5, type = list("CAS", "SRD"))
  expect_true("lag" %in% names(res$ACF_fitted))
  expect_true("CAS" %in% names(res$ACF_fitted))
  expect_true("SRD" %in% names(res$ACF_fitted))
})

test_that("fit_ACF() runs without error when ignore_zeros=TRUE", {
  ts <- make_zeros_ts()
  expect_no_error(fit_ACF(ts, lag_max = 5, ignore_zeros = TRUE))
})

# ---------------------------------------------------------------------------
# fit_ACF_monthly
# ---------------------------------------------------------------------------
test_that("fit_ACF_monthly() returns list with ACF_params_monthly and ACF_monthly_plot", {
  ts  <- make_ar1_ts()
  res <- fit_ACF_monthly(ts, lag = 5, type = list("CAS", "SRD"))
  expect_type(res, "list")
  expect_named(res, c("ACF_params_monthly", "ACF_monthly_plot"), ignore.order = TRUE)
})

test_that("fit_ACF_monthly() ACF_params_monthly has one row per month present in data", {
  ts  <- make_ar1_ts()
  res <- fit_ACF_monthly(ts, lag = 5)
  expect_true(nrow(res$ACF_params_monthly) <= 12L)
  expect_true(nrow(res$ACF_params_monthly) >= 1L)
})

test_that("fit_ACF_monthly() parameter values are all numeric", {
  ts  <- make_ar1_ts()
  res <- fit_ACF_monthly(ts, lag = 5)
  expect_true(all(is.numeric(unlist(res$ACF_params_monthly))))
})
