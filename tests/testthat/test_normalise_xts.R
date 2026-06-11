library(testthat)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
make_daily_ts <- function() {
  set.seed(50)
  n     <- 365 * 5
  vals  <- abs(rnorm(n, 5, 2))
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "P")),
           order.by = dates)
}

make_multi_ts <- function() {
  set.seed(51)
  n     <- 365 * 5
  vals  <- matrix(abs(rnorm(n * 2, 5, 2)), ncol = 2,
                  dimnames = list(NULL, c("P1", "P2")))
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  xts::xts(vals, order.by = dates)
}

make_zeros_ts <- function() {
  ts  <- make_daily_ts()
  set.seed(52)
  idx <- sample(nrow(ts), size = round(nrow(ts) * 0.15))
  ts[idx, ] <- 0
  ts
}

# ---------------------------------------------------------------------------
# normalise_xts — return structure
# ---------------------------------------------------------------------------
test_that("normalise_xts() returns a list", {
  ts  <- make_daily_ts()
  res <- normalise_xts(ts)
  expect_type(res, "list")
})

test_that("normalise_xts() list length equals ncol(ts)", {
  ts  <- make_multi_ts()
  res <- normalise_xts(ts)
  expect_equal(length(res), ncol(ts))
})

test_that("normalise_xts() list names match colnames(ts)", {
  ts  <- make_multi_ts()
  res <- normalise_xts(ts)
  expect_equal(names(res), colnames(ts))
})

test_that("normalise_xts() each element is an xts object", {
  ts  <- make_daily_ts()
  res <- normalise_xts(ts)
  expect_true(inherits(res[[1]], "xts"))
})

# ---------------------------------------------------------------------------
# normalise_xts — monthly normalization
# ---------------------------------------------------------------------------
test_that("normalise_xts() monthly mode returns same row count as input", {
  ts  <- make_daily_ts()
  res <- normalise_xts(ts, dist_period = "monthly")
  expect_equal(nrow(res[[1]]), nrow(ts))
})

test_that("normalise_xts() monthly mode returns numeric values", {
  ts  <- make_daily_ts()
  res <- normalise_xts(ts, dist_period = "monthly")
  vals <- as.numeric(res[[1]])
  expect_true(all(is.numeric(vals)))
})

test_that("normalise_xts() monthly mode values are not all NA", {
  ts  <- make_daily_ts()
  res <- normalise_xts(ts, dist_period = "monthly")
  expect_false(all(is.na(as.numeric(res[[1]]))))
})

# ---------------------------------------------------------------------------
# normalise_xts — global normalization
# ---------------------------------------------------------------------------
test_that("normalise_xts() with dist_period=NA runs without error", {
  ts <- make_daily_ts()
  expect_no_error(normalise_xts(ts, dist_period = NA))
})

test_that("normalise_xts() global mode returns same row count as input", {
  ts  <- make_daily_ts()
  res <- normalise_xts(ts, dist_period = NA)
  expect_equal(nrow(res[[1]]), nrow(ts))
})

test_that("normalise_xts() global mode values are not all NA", {
  ts  <- make_daily_ts()
  res <- normalise_xts(ts, dist_period = NA)
  expect_false(all(is.na(as.numeric(res[[1]]))))
})

# ---------------------------------------------------------------------------
# normalise_xts — multi-column
# ---------------------------------------------------------------------------
test_that("normalise_xts() produces independent normalizations per column", {
  ts   <- make_multi_ts()
  res  <- normalise_xts(ts)
  # Mean of normalized series should be near 0 for each column
  m1   <- mean(as.numeric(res[[1]]), na.rm = TRUE)
  m2   <- mean(as.numeric(res[[2]]), na.rm = TRUE)
  expect_lt(abs(m1), 0.5)
  expect_lt(abs(m2), 0.5)
})

# ---------------------------------------------------------------------------
# normalise_xts — ignore_zeros
# ---------------------------------------------------------------------------
test_that("normalise_xts() runs without error when ignore_zeros=TRUE", {
  ts <- make_zeros_ts()
  expect_no_error(normalise_xts(ts, ignore_zeros = TRUE))
})

test_that("normalise_xts() ignore_zeros=TRUE returns fewer rows than original", {
  ts  <- make_zeros_ts()
  res <- normalise_xts(ts, ignore_zeros = TRUE)
  expect_lt(nrow(res[[1]]), nrow(ts))
})
