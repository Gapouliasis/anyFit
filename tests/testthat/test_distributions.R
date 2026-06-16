library(testthat)

# ===========================================================================
# Unit tests for every distribution in R/Distributions_Anyfit.R
#
#   Section A  d/p/q/r consistency for each distribution family
#   Section B  fitlm_* L-Moment fitters
#   Section C  cross-distribution helpers (multi / monthly / diagnostics)
# ===========================================================================

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
make_norm_ts <- function() {
  set.seed(42)
  vals  <- rnorm(500, mean = 5, sd = 2)
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, 499) * 86400
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "X")), order.by = dates)
}

make_gamma_ts <- function() {
  set.seed(43)
  vals  <- rgamma(500, shape = 2, scale = 3)
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, 499) * 86400
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "X")), order.by = dates)
}

make_zeros_ts <- function() {
  set.seed(44)
  ts  <- make_norm_ts()
  idx <- sample(500, size = 100)
  ts[idx, ] <- 0
  ts
}

make_multi_ts <- function() {
  set.seed(45)
  n     <- 500
  vals  <- matrix(abs(rnorm(n * 2, 5, 2)), ncol = 2,
                  dimnames = list(NULL, c("A", "B")))
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  xts::xts(vals, order.by = dates)
}

make_all_zeros_ts <- function() {
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, 99) * 86400
  xts::xts(matrix(0, nrow = 100, ncol = 1,
                  dimnames = list(NULL, "Z")), order.by = dates)
}

make_mixed_zeros_ts <- function() {
  n     <- 100
  vals  <- matrix(c(abs(rnorm(n, 5, 2)), rep(0, n)), ncol = 2,
                  dimnames = list(NULL, c("Good", "AllZero")))
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, n - 1) * 86400
  xts::xts(vals, order.by = dates)
}

# Helper: check that a fitlm result has the standard structure
check_fitlm_structure <- function(res) {
  expect_type(res, "list")
  expect_true(all(c("Distribution", "Param", "TheorLMom", "DataLMom", "GoF") %in% names(res)))
  gof_names <- c("MLE", "CM", "KS", "MSEquant", "DiffOfMax", "MeanDiffOf10Max")
  expect_true(all(gof_names %in% names(res$GoF)))
  params <- unlist(res$Param)
  expect_true(all(is.numeric(params)))
  expect_true(all(is.finite(params)))
}

# ===========================================================================
# Section A — d/p/q/r consistency checks
# ===========================================================================

# --- Exponential (location, scale) ---
test_that("dexp() returns non-negative density values", {
  x  <- seq(0.1, 5, by = 0.1)
  fx <- dexp(x, location = 0, scale = 2)
  expect_true(all(fx >= 0))
  expect_equal(length(fx), length(x))
})

test_that("pexp() is monotone increasing and bounded [0, 1]", {
  q  <- seq(0, 10, by = 0.5)
  px <- pexp(q, location = 0, scale = 2)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qexp(pexp(x)) round-trip returns original values", {
  x   <- seq(0.5, 5, by = 0.5)
  loc <- 0; sc <- 2
  expect_equal(qexp(pexp(x, loc, sc), loc, sc), x, tolerance = 1e-6)
})

test_that("rexp() returns n samples all >= location", {
  samp <- rexp(100, location = 1, scale = 2)
  expect_equal(length(samp), 100L)
  expect_true(all(samp >= 1))
})

# --- Rayleigh (location, scale) ---
test_that("drayleigh() returns non-negative density values", {
  x  <- seq(0.1, 5, by = 0.1)
  fx <- drayleigh(x, location = 0, scale = 1)
  expect_true(all(fx >= 0))
  expect_equal(length(fx), length(x))
})

test_that("prayleigh() is monotone and bounded [0, 1]", {
  q  <- seq(0, 5, by = 0.5)
  px <- prayleigh(q, location = 0, scale = 1)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qrayleigh(prayleigh(x)) round-trip returns original values", {
  x <- seq(0.5, 4, by = 0.5)
  expect_equal(qrayleigh(prayleigh(x, 0, 1), 0, 1), x, tolerance = 1e-5)
})

test_that("rrayleigh() returns n samples >= location", {
  samp <- rrayleigh(100, location = 0, scale = 1)
  expect_equal(length(samp), 100L)
  expect_true(all(samp >= 0))
})

# --- Gamma (scale, shape) ---
test_that("dgamma() returns non-negative density values", {
  x  <- seq(0.1, 10, by = 0.1)
  fx <- dgamma(x, scale = 3, shape = 2)
  expect_true(all(fx >= 0))
  expect_equal(length(fx), length(x))
})

test_that("pgamma() is monotone and bounded [0, 1]", {
  q  <- seq(0, 20, by = 0.5)
  px <- pgamma(q, scale = 3, shape = 2)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qgamma(pgamma(x)) round-trip returns original values", {
  x <- seq(0.5, 12, by = 0.5)
  expect_equal(qgamma(pgamma(x, 3, 2), 3, 2), x, tolerance = 1e-5)
})

test_that("rgamma() returns n positive samples", {
  samp <- rgamma(100, scale = 3, shape = 2)
  expect_equal(length(samp), 100L)
  expect_true(all(samp > 0))
})

# --- Gamma3 (location, scale, shape) ---
test_that("dgamma3() returns non-negative density values", {
  x  <- seq(2, 10, by = 0.5)
  fx <- dgamma3(x, location = 1, scale = 2, shape = 3)
  expect_true(all(fx >= 0))
})

test_that("pgamma3() is monotone increasing and bounded [0, 1]", {
  q  <- seq(2, 10, by = 0.5)
  px <- pgamma3(q, location = 1, scale = 2, shape = 3)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qgamma3(pgamma3(x)) round-trip returns original values", {
  x   <- seq(2, 8, by = 0.5)
  loc <- 1; sc <- 2; sh <- 3
  expect_equal(qgamma3(pgamma3(x, loc, sc, sh), loc, sc, sh), x, tolerance = 1e-4)
})

test_that("rgamma3() returns n samples", {
  samp <- rgamma3(100, location = 1, scale = 2, shape = 3)
  expect_equal(length(samp), 100L)
})

# --- Generalised Logistic (location, scale, shape) — d/p/q only, no r ---
test_that("dgenlogi() returns non-negative density values", {
  x  <- seq(-5, 5, by = 0.25)
  fx <- dgenlogi(x, location = 0, scale = 1, shape = -0.5)
  expect_true(all(fx >= 0))
  expect_equal(length(fx), length(x))
})

test_that("pgenlogi() is monotone and bounded [0, 1]", {
  q  <- seq(-5, 5, by = 0.25)
  px <- pgenlogi(q, location = 0, scale = 1, shape = -0.5)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qgenlogi(pgenlogi(x)) round-trip returns original values", {
  # shape = -0.5 bounds the GLO support below at -2, so stay strictly inside it
  x <- seq(-1.5, 3, by = 0.5)
  expect_equal(qgenlogi(pgenlogi(x, 0, 1, -0.5), 0, 1, -0.5), x, tolerance = 1e-5)
})

# --- Normal (mean, sd) ---
test_that("dnorm() returns non-negative density values", {
  x  <- seq(-5, 15, by = 0.25)
  fx <- dnorm(x, mean = 5, sd = 2)
  expect_true(all(fx >= 0))
  expect_equal(length(fx), length(x))
})

test_that("pnorm() is monotone and bounded [0, 1]", {
  q  <- seq(-5, 15, by = 0.25)
  px <- pnorm(q, mean = 5, sd = 2)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qnorm(pnorm(x)) round-trip returns original values", {
  x <- seq(-2, 12, by = 0.5)
  expect_equal(qnorm(pnorm(x, 5, 2), 5, 2), x, tolerance = 1e-6)
})

test_that("rnorm() returns n samples", {
  samp <- rnorm(100, mean = 5, sd = 2)
  expect_equal(length(samp), 100L)
})

# --- Weibull3 (location, scale, shape) ---
test_that("dweibull() returns non-negative density values", {
  x  <- seq(0.1, 10, by = 0.1)
  fx <- dweibull(x, location = 0, scale = 5, shape = 1.5)
  expect_true(all(fx >= 0))
  expect_equal(length(fx), length(x))
})

test_that("pweibull() is monotone and bounded [0, 1]", {
  q  <- seq(0, 15, by = 0.5)
  px <- pweibull(q, location = 0, scale = 5, shape = 1.5)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qweibull(pweibull(x)) round-trip returns original values", {
  x <- seq(0.5, 10, by = 0.5)
  expect_equal(qweibull(pweibull(x, 0, 5, 1.5), 0, 5, 1.5), x, tolerance = 1e-5)
})

test_that("rweibull() returns n samples >= location", {
  samp <- rweibull(100, location = 0, scale = 5, shape = 1.5)
  expect_equal(length(samp), 100L)
  expect_true(all(samp >= 0))
})

# --- Gumbel (location, scale) ---
test_that("dgumbel() returns non-negative density values", {
  x  <- seq(0, 10, by = 0.5)
  fx <- dgumbel(x, location = 3, scale = 2)
  expect_true(all(fx >= 0))
})

test_that("pgumbel() is monotone and bounded [0, 1]", {
  q  <- seq(0, 10, by = 0.5)
  px <- pgumbel(q, location = 3, scale = 2)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qgumbel(pgumbel(x)) round-trip returns original values", {
  x <- seq(1, 8, by = 0.5)
  expect_equal(qgumbel(pgumbel(x, 3, 2), 3, 2), x, tolerance = 1e-5)
})

test_that("rgumbel() returns n samples", {
  samp <- rgumbel(100, location = 3, scale = 2)
  expect_equal(length(samp), 100L)
})

# --- Log-Normal (location, scale, shape) ---
test_that("dlognorm() returns non-negative density values", {
  x  <- seq(0.1, 5, by = 0.1)
  fx <- dlognorm(x, location = 0, scale = 1, shape = 0.5)
  expect_true(all(fx >= 0))
})

test_that("plognorm() is monotone and bounded [0, 1]", {
  q  <- seq(0.1, 10, by = 0.5)
  px <- plognorm(q, location = 0, scale = 1, shape = 0.5)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qlognorm(plognorm(x)) round-trip returns original values", {
  x <- seq(0.5, 5, by = 0.5)
  expect_equal(qlognorm(plognorm(x, 0, 1, 0.5), 0, 1, 0.5), x, tolerance = 1e-5)
})

test_that("rlognorm() returns n positive samples", {
  samp <- rlognorm(100, location = 0, scale = 1, shape = 0.5)
  expect_equal(length(samp), 100L)
  expect_true(all(samp > 0))
})

# --- GEV (location, scale, shape) ; note qgev(p, shape, scale, location) ---
test_that("dgev() returns non-negative density values", {
  x  <- seq(0, 10, by = 0.5)
  fx <- dgev(x, location = 5, scale = 2, shape = 0.1)
  expect_true(all(fx >= 0))
})

test_that("pgev() is monotone and bounded [0, 1]", {
  q  <- seq(0, 10, by = 0.5)
  px <- pgev(q, location = 5, scale = 2, shape = 0.1)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qgev(pgev(x)) round-trip returns original values", {
  x <- seq(2, 9, by = 0.5)
  # qgev signature: qgev(p, shape, scale, location) — different arg order
  expect_equal(qgev(pgev(x, 5, 2, 0.1), 0.1, 2, 5), x, tolerance = 1e-5)
})

test_that("rgev() returns n samples", {
  samp <- rgev(100, location = 5, scale = 2, shape = 0.1)
  expect_equal(length(samp), 100L)
})

# --- Generalized Pareto (location, scale, shape) ---
# lmom GPA: shape = 0.2 bounds support above at location + scale/shape = 10
test_that("dgpd() returns non-negative density values", {
  x  <- seq(0, 9, by = 0.5)
  fx <- dgpd(x, location = 0, scale = 2, shape = 0.2)
  expect_true(all(fx >= 0))
  expect_equal(length(fx), length(x))
})

test_that("pgpd() is monotone and bounded [0, 1]", {
  q  <- seq(0, 9, by = 0.5)
  px <- pgpd(q, location = 0, scale = 2, shape = 0.2)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qgpd(pgpd(x)) round-trip returns original values", {
  x <- seq(0.5, 9, by = 0.5)
  expect_equal(qgpd(pgpd(x, 0, 2, 0.2), 0, 2, 0.2), x, tolerance = 1e-5)
})

test_that("rgpd() returns n samples >= location", {
  samp <- rgpd(100, location = 0, scale = 2, shape = 0.2)
  expect_equal(length(samp), 100L)
  expect_true(all(samp >= 0))
})

# --- Generalised Gamma (scale, shape1, shape2) ---
test_that("dgengamma() returns non-negative density values", {
  skip_if_not_installed("VGAM")
  x  <- seq(0.1, 20, by = 0.2)
  fx <- dgengamma(x, scale = 5, shape1 = 2, shape2 = 1.5)
  expect_true(all(fx >= 0))
  expect_equal(length(fx), length(x))
})

test_that("pgengamma() is monotone and bounded [0, 1]", {
  skip_if_not_installed("VGAM")
  q  <- seq(0.1, 30, by = 0.5)
  px <- pgengamma(q, scale = 5, shape1 = 2, shape2 = 1.5)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qgengamma(pgengamma(x)) round-trip returns original values", {
  skip_if_not_installed("VGAM")
  x <- seq(1, 15, by = 0.5)
  expect_equal(qgengamma(pgengamma(x, 5, 2, 1.5), 5, 2, 1.5), x, tolerance = 1e-4)
})

test_that("rgengamma() returns n positive samples", {
  skip_if_not_installed("VGAM")
  samp <- rgengamma(100, scale = 5, shape1 = 2, shape2 = 1.5)
  expect_equal(length(samp), 100L)
  expect_true(all(samp > 0))
})

# --- Generalised Gamma with location (location, scale, shape1, shape2) ---
test_that("dgengamma_loc() returns non-negative density values", {
  skip_if_not_installed("VGAM")
  x  <- seq(2.1, 22, by = 0.2)
  fx <- dgengamma_loc(x, location = 2, scale = 5, shape1 = 2, shape2 = 1.5)
  expect_true(all(fx >= 0))
  expect_equal(length(fx), length(x))
})

test_that("pgengamma_loc() is monotone and bounded [0, 1]", {
  skip_if_not_installed("VGAM")
  q  <- seq(2.1, 32, by = 0.5)
  px <- pgengamma_loc(q, location = 2, scale = 5, shape1 = 2, shape2 = 1.5)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qgengamma_loc(pgengamma_loc(x)) round-trip returns original values", {
  skip_if_not_installed("VGAM")
  x <- seq(3, 17, by = 0.5)
  expect_equal(
    qgengamma_loc(pgengamma_loc(x, 2, 5, 2, 1.5), 2, 5, 2, 1.5),
    x, tolerance = 1e-4
  )
})

test_that("rgengamma_loc() returns n samples >= location", {
  skip_if_not_installed("VGAM")
  samp <- rgengamma_loc(100, location = 2, scale = 5, shape1 = 2, shape2 = 1.5)
  expect_equal(length(samp), 100L)
  expect_true(all(samp >= 2))
})

# --- Burr XII (scale, shape1, shape2, PW=1) — no location param ---
test_that("dburr() returns non-negative density values", {
  x  <- seq(0.1, 5, by = 0.1)
  fx <- dburr(x, scale = 2, shape1 = 1.5, shape2 = 0.5)
  expect_true(all(fx >= 0))
})

test_that("pburr() is monotone and bounded [0, 1]", {
  q  <- seq(0.1, 5, by = 0.1)
  px <- pburr(q, scale = 2, shape1 = 1.5, shape2 = 0.5)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qburr(pburr(x)) round-trip returns original values", {
  x <- seq(0.5, 4, by = 0.5)
  expect_equal(qburr(pburr(x, 2, 1.5, 0.5), 2, 1.5, 0.5), x, tolerance = 1e-5)
})

test_that("rburr() returns n positive samples", {
  samp <- rburr(100, scale = 2, shape1 = 1.5, shape2 = 0.5)
  expect_equal(length(samp), 100L)
  expect_true(all(samp > 0))
})

# --- Dagum (scale, shape1, shape2, PW=1) ---
test_that("ddagum() returns non-negative density values", {
  x  <- seq(0.1, 10, by = 0.1)
  fx <- ddagum(x, scale = 2, shape1 = 1.5, shape2 = 0.3)
  expect_true(all(fx >= 0))
})

test_that("pdagum() is monotone and bounded [0, 1]", {
  q  <- seq(0.1, 20, by = 0.2)
  px <- pdagum(q, scale = 2, shape1 = 1.5, shape2 = 0.3)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qdagum(pdagum(x)) round-trip returns original values", {
  x <- seq(0.5, 8, by = 0.5)
  expect_equal(qdagum(pdagum(x, 2, 1.5, 0.3), 2, 1.5, 0.3), x, tolerance = 1e-4)
})

test_that("rdagum() returns n positive samples", {
  samp <- rdagum(100, scale = 2, shape1 = 1.5, shape2 = 0.3)
  expect_equal(length(samp), 100L)
  expect_true(all(samp > 0))
})

# --- Exponentiated Weibull (scale, shape1, shape2) ---
test_that("dexpweibull() returns non-negative density values", {
  x  <- seq(0.1, 10, by = 0.1)
  fx <- dexpweibull(x, scale = 5, shape1 = 1.5, shape2 = 2)
  expect_true(all(fx >= 0))
  expect_equal(length(fx), length(x))
})

test_that("pexpweibull() is monotone and bounded [0, 1]", {
  q  <- seq(0, 15, by = 0.5)
  px <- pexpweibull(q, scale = 5, shape1 = 1.5, shape2 = 2)
  expect_true(all(diff(px) >= 0))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("qexpweibull(pexpweibull(x)) round-trip returns original values", {
  x <- seq(0.5, 10, by = 0.5)
  expect_equal(
    qexpweibull(pexpweibull(x, 5, 1.5, 2), 5, 1.5, 2), x, tolerance = 1e-5
  )
})

test_that("rexpweibull() returns n positive samples", {
  samp <- rexpweibull(100, scale = 5, shape1 = 1.5, shape2 = 2)
  expect_equal(length(samp), 100L)
  expect_true(all(samp > 0))
})

# ===========================================================================
# Section B — fitlm_* single-series fitters
# ===========================================================================

test_that("fitlm_norm() returns standard structure with finite params", {
  res <- fitlm_norm(make_norm_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_norm() runs without error with ignore_zeros=TRUE", {
  expect_no_error(fitlm_norm(make_zeros_ts(), ignore_zeros = TRUE))
})

test_that("fitlm_exp() returns standard structure with finite params", {
  res <- fitlm_exp(make_gamma_ts())   # positive values needed for exponential
  check_fitlm_structure(res)
})

test_that("fitlm_exp() runs without error with ignore_zeros=TRUE", {
  expect_no_error(fitlm_exp(make_zeros_ts(), ignore_zeros = TRUE))
})

test_that("fitlm_gamma() returns standard structure with finite params", {
  res <- fitlm_gamma(make_gamma_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_gamma3() returns standard structure with finite params", {
  res <- fitlm_gamma3(make_gamma_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_genlogi() returns standard structure with finite params", {
  res <- fitlm_genlogi(make_norm_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_rayleigh() returns standard structure with finite params", {
  res <- fitlm_rayleigh(make_gamma_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_weibull() returns standard structure with finite params", {
  res <- fitlm_weibull(make_gamma_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_gumbel() returns standard structure with finite params", {
  res <- fitlm_gumbel(make_norm_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_lognorm() returns standard structure with finite params", {
  res <- fitlm_lognorm(make_gamma_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_gev() returns standard structure with finite params", {
  res <- fitlm_gev(make_norm_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_gengamma() returns standard structure with finite params", {
  skip_if_not_installed("VGAM")
  res <- fitlm_gengamma(make_gamma_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_burr() returns standard structure with finite params", {
  res <- fitlm_burr(make_gamma_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_dagum() returns standard structure with finite params", {
  res <- fitlm_dagum(make_gamma_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_expweibull() returns standard structure with finite params", {
  res <- fitlm_expweibull(make_gamma_ts())
  check_fitlm_structure(res)
})

# GPD now ships self-contained dgpd/pgpd/qgpd/rgpd wrappers (lmom GPA params).
test_that("fitlm_GPD() returns standard structure with finite params", {
  res <- fitlm_GPD(make_gamma_ts())
  check_fitlm_structure(res)
})

test_that("fitlm_gengamma_loc() returns standard structure with finite params", {
  skip_if_not_installed("VGAM")
  res <- fitlm_gengamma_loc(make_gamma_ts(), location = 0)
  check_fitlm_structure(res)
})

# ===========================================================================
# Section C — cross-distribution helpers
# ===========================================================================

test_that("fitlm_multi() with diagnostic_plots=TRUE returns 5 list elements", {
  res <- fitlm_multi(make_norm_ts(), candidates = list("norm", "gamma"))
  expect_type(res, "list")
  expect_true(all(c("parameter_list", "GoF_summary", "diagnostics",
                    "QQplot", "PPplot") %in% names(res)))
})

test_that("fitlm_multi() GoF_summary has one column per candidate", {
  res <- fitlm_multi(make_norm_ts(), candidates = list("norm", "gamma"))
  expect_equal(ncol(res$GoF_summary), 2L)
  expect_equal(colnames(res$GoF_summary), c("norm", "gamma"))
})

test_that("fitlm_multi() with diagnostic_plots=FALSE returns only params and GoF", {
  res <- fitlm_multi(make_norm_ts(), candidates = list("norm"), diagnostic_plots = FALSE)
  expect_true("parameter_list" %in% names(res))
  expect_true("GoF_summary"    %in% names(res))
  expect_false("QQplot" %in% names(res))
})

test_that("fitlm_multi() with ignore_zeros=TRUE handles all-zero data without error", {
  expect_no_error({
    res <- fitlm_multi(make_all_zeros_ts(), candidates = list("gev"),
                       ignore_zeros = TRUE, diagnostic_plots = FALSE)
  })
  expect_type(res, "list")
  expect_true(all(c("parameter_list", "GoF_summary") %in% names(res)))
  expect_true(is.na(res$parameter_list$gev$Param))
  expect_true(all(is.na(res$parameter_list$gev$TheorLMom)))
})

test_that("fitlm_multi() with ignore_zeros=TRUE and diagnostic_plots=TRUE returns 5 elements for empty data", {
  res <- fitlm_multi(make_all_zeros_ts(), candidates = list("gev"),
                     ignore_zeros = TRUE, diagnostic_plots = TRUE)
  expect_type(res, "list")
  expect_true(all(c("parameter_list", "GoF_summary", "diagnostics",
                    "QQplot", "PPplot") %in% names(res)))
  expect_null(res$diagnostics)
  expect_null(res$QQplot)
  expect_null(res$PPplot)
})

test_that("fitlm_multi() with ignore_zeros=FALSE on all-zero data errors as expected (degenerate data)", {
  expect_error(
    fitlm_multi(make_all_zeros_ts(), candidates = list("gev"),
                ignore_zeros = FALSE, diagnostic_plots = FALSE)
  )
})

test_that("fitlm_nxts() returns list with params, diagnostic_plots, QQ_plots, PP_plots", {
  res <- fitlm_nxts(make_multi_ts(), candidates = list("norm"), diagnostic_plots = TRUE)
  expect_true(all(c("params", "diagnostic_plots", "QQ_plots", "PP_plots") %in% names(res)))
})

test_that("fitlm_nxts() params has one entry per column of ts", {
  ts  <- make_multi_ts()
  res <- fitlm_nxts(ts, candidates = list("norm"), diagnostic_plots = FALSE)
  expect_equal(length(res$params), ncol(ts))
  expect_equal(names(res$params), colnames(ts))
})

test_that("fitlm_nxts() handles mixed valid and all-zero columns with ignore_zeros=TRUE", {
  ts <- make_mixed_zeros_ts()
  expect_no_error({
    res <- fitlm_nxts(ts, candidates = list("gev"), ignore_zeros = TRUE,
                      diagnostic_plots = FALSE)
  })
  expect_equal(length(res$params), ncol(ts))
  expect_equal(names(res$params), c("Good", "AllZero"))
  expect_true(all(!is.na(unlist(res$params$Good$params$gev$Param))))
  expect_true(is.na(res$params$AllZero$params$gev$Param))
})

test_that("fitlm_monthly() returns params_monthly, GoF_monthly, monthly_QQplot, monthly_PPplot", {
  res <- fitlm_monthly(make_norm_ts(), candidates = list("norm"))
  expect_true(all(c("params_monthly", "GoF_monthly",
                    "monthly_QQplot", "monthly_PPplot") %in% names(res)))
})

test_that("fitlm_monthly() params_monthly is a list with one entry per candidate", {
  res <- fitlm_monthly(make_norm_ts(), candidates = list("norm", "gamma"))
  expect_equal(length(res$params_monthly), 2L)
  expect_equal(names(res$params_monthly), c("norm", "gamma"))
})

test_that("fitlm_monthly() params_monthly columns are month names", {
  res <- fitlm_monthly(make_norm_ts(), candidates = list("norm"))
  expect_true(any(colnames(res$params_monthly[[1]]) %in% month.name))
})

test_that("fit_diagnostics() returns Diagnostic_Plots, GoF, QQplot, PPplot", {
  ts  <- make_norm_ts()
  fit <- fitlm_norm(ts)
  res <- fit_diagnostics(ts, dist = "norm", params = fit$Param)
  expect_type(res, "list")
  expect_named(res, c("Diagnostic_Plots", "GoF", "QQplot", "PPplot"),
               ignore.order = TRUE)
})

test_that("fit_diagnostics() GoF contains CramerVonMises and KolmogorovSmirnov", {
  ts  <- make_norm_ts()
  fit <- fitlm_norm(ts)
  res <- fit_diagnostics(ts, dist = "norm", params = fit$Param)
  expect_true(all(c("CramerVonMises", "KolmogorovSmirnov") %in% names(res$GoF)))
})

test_that("cvm_two_sample / ks_two_sample match CDFt reference values (with ties)", {
  set.seed(42)
  S1 <- round(stats::rgamma(500, shape = 2,   scale = 3), 1)   # rainfall-like, heavy ties
  S2 <- round(stats::rgamma(500, shape = 2.3, scale = 3), 1)
  # Reference values computed once from CDFt::CramerVonMisesTwoSamples /
  # CDFt::KolmogorovSmirnov; CDFt has since been dropped as a dependency and these
  # vectorized helpers replace it (verified bit-identical, incl. ties).
  expect_equal(cvm_two_sample(S1, S2), 2.5041039999999839, tolerance = 1e-12)
  expect_equal(ks_two_sample(S1, S2),  0.15600000000000003, tolerance = 1e-12)
})

test_that("fit_diagnostics() QQplot and PPplot are ggplot objects", {
  ts  <- make_norm_ts()
  fit <- fitlm_norm(ts)
  res <- fit_diagnostics(ts, dist = "norm", params = fit$Param)
  expect_true(inherits(res$QQplot, "gg"))
  expect_true(inherits(res$PPplot, "gg"))
})

test_that("LRatio_check() returns list with 'distributions' and 'multi_plots'", {
  lmoms <- lmom_stats(make_multi_ts())
  res   <- LRatio_check(lmoms)
  expect_type(res, "list")
  expect_named(res, c("distributions", "multi_plots"), ignore.order = TRUE)
})

test_that("LRatio_check() distributions contains entries for the four tested distributions", {
  lm  <- lmom_stats(make_multi_ts())
  res <- LRatio_check(lm)
  dist_names <- names(res$distributions)
  expect_true(all(c("Dagum", "GGamma", "ExpWeibull", "BurrXII") %in% dist_names))
})
