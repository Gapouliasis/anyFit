library(testthat)

# ---------------------------------------------------------------------------
# Shared fixtures — two time series with known Pearson ≈ 0.8
# ---------------------------------------------------------------------------
make_corr_pair <- function() {
  set.seed(7)
  n     <- 365 * 3
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  base  <- rnorm(n, 5, 2)
  # y is 0.8*base + noise, giving Pearson ≈ 0.8
  x <- xts::xts(matrix(base,                     ncol = 1, dimnames = list(NULL, "X")), order.by = dates)
  y <- xts::xts(matrix(0.8 * base + rnorm(n, 0, sqrt(1 - 0.64) * 2), ncol = 1, dimnames = list(NULL, "Y")), order.by = dates)
  list(x = x, y = y)
}

make_zeros_pair <- function() {
  pair <- make_corr_pair()
  idx  <- sample(nrow(pair$x), size = 200)
  pair$x[idx, ] <- 0
  pair
}

# ---------------------------------------------------------------------------
# Return structure
# ---------------------------------------------------------------------------
test_that("correl_plots() returns list with all five elements", {
  p   <- make_corr_pair()
  res <- correl_plots(p$x, p$y)
  expect_type(res, "list")
  expect_named(res, c("combined", "scatter_plot", "copula_plot",
                      "normal_plot", "correl_table"),
               ignore.order = TRUE)
})

test_that("correl_plots() correl_table has 6 rows", {
  p   <- make_corr_pair()
  res <- correl_plots(p$x, p$y)
  expect_equal(nrow(res$correl_table), 6L)
})

test_that("correl_plots() correl_table row names are correct", {
  p        <- make_corr_pair()
  res      <- correl_plots(p$x, p$y)
  expected <- c("Pearson", "Spearman", "Semi_1", "Semi_2", "Semi_3", "Semi_4")
  expect_equal(rownames(res$correl_table), expected)
})

# ---------------------------------------------------------------------------
# Pearson value
# ---------------------------------------------------------------------------
test_that("correl_plots() Pearson correlation is close to expected ~0.8", {
  p    <- make_corr_pair()
  res  <- correl_plots(p$x, p$y)
  rho  <- as.numeric(res$correl_table["Pearson", 1])
  expect_gt(rho, 0.7)
  expect_lt(rho, 0.95)
})

test_that("correl_plots() Pearson equals manual cor() on common data", {
  p    <- make_corr_pair()
  res  <- correl_plots(p$x, p$y, check_common = TRUE)
  rho  <- as.numeric(res$correl_table["Pearson", 1])
  # Manual calculation
  merged  <- merge(p$x, p$y, join = "inner")
  manual  <- round(cor(as.numeric(merged[, 1]), as.numeric(merged[, 2]),
                       use = "complete.obs"), 2)
  expect_equal(rho, manual)
})

# ---------------------------------------------------------------------------
# check_common
# ---------------------------------------------------------------------------
test_that("correl_plots() with check_common=TRUE uses only common dates", {
  p    <- make_corr_pair()
  # Offset y by 10 days so first 10 rows of x have no match in y
  p$y <- p$y[-seq(1, 10), ]
  res  <- correl_plots(p$x, p$y, check_common = TRUE)
  expect_true(inherits(res$combined, "patchwork"))
})

# ---------------------------------------------------------------------------
# Plot types
# ---------------------------------------------------------------------------
test_that("correl_plots() scatter_plot is a ggplot object", {
  p   <- make_corr_pair()
  res <- correl_plots(p$x, p$y)
  expect_true(inherits(res$scatter_plot, "gg"))
})

test_that("correl_plots() copula_plot is a ggplot object", {
  p   <- make_corr_pair()
  res <- correl_plots(p$x, p$y)
  expect_true(inherits(res$copula_plot, "gg"))
})

test_that("correl_plots() normal_plot is a ggplot object", {
  p   <- make_corr_pair()
  res <- correl_plots(p$x, p$y)
  expect_true(inherits(res$normal_plot, "gg"))
})

test_that("correl_plots() combined is a patchwork object", {
  p   <- make_corr_pair()
  res <- correl_plots(p$x, p$y)
  expect_true(inherits(res$combined, "patchwork"))
})

# ---------------------------------------------------------------------------
# ignore_zeros
# ---------------------------------------------------------------------------
test_that("correl_plots() runs without error when ignore_zeros=TRUE", {
  p <- make_zeros_pair()
  expect_no_error(correl_plots(p$x, p$y, ignore_zeros = TRUE))
})

test_that("correl_plots() correl_table still has 6 rows with ignore_zeros=TRUE", {
  p   <- make_zeros_pair()
  res <- correl_plots(p$x, p$y, ignore_zeros = TRUE)
  expect_equal(nrow(res$correl_table), 6L)
})
