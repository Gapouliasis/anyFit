library(testthat)

# ===========================================================================
# Quadrature L-moment engine lmom_gengamma (R/lmom_closedform.R) and its use
# inside fitlm_gengamma / fitlm_gengamma_loc.
# ===========================================================================

# Independent oracle: beta_r = integral_0^1 q(F) F^r dF via the package quantile,
# then PWM -> L-moments. A different computation path from the quadrature engine.
gg_pwm_to_lmom <- function(beta) {
  rmax <- length(beta); lambda <- numeric(rmax)
  for (r in 1:rmax) { j <- 0:(r - 1)
    lambda[r] <- sum((-1)^(r - 1 - j) * choose(r - 1, j) * choose(r - 1 + j, j) * beta[j + 1]) }
  c(lambda[1], lambda[2], lambda[3] / lambda[2], lambda[4] / lambda[2], lambda[5] / lambda[2])
}
gg_oracle <- function(scale, shape1, shape2, rmax = 5) {
  beta <- vapply(0:(rmax - 1), function(r)
    integrate(function(F) anyFit::qgengamma(F, scale, shape1, shape2) * F^r, 0, 1,
              rel.tol = 1e-10, subdivisions = 4000L)$value, numeric(1))
  gg_pwm_to_lmom(beta)
}

make_pos_ts <- function(vals) {
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq_along(vals) * 86400
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "X")), order.by = dates)
}

# ---------------------------------------------------------------------------
# Engine shape / layout
# ---------------------------------------------------------------------------
test_that("lmom_gengamma returns a length-5 named vector for order 1:5", {
  g <- anyFit:::lmom_gengamma(1:5, scale = 2, shape1 = 1.5, shape2 = 0.8)
  expect_length(g, 5L)
  expect_equal(names(g), c("lambda_1", "lambda_2", "tau_3", "tau_4", "tau_5"))
  expect_true(all(is.finite(g)))
})

test_that("lmom_gengamma returns NA-named vector for degenerate shapes", {
  g <- anyFit:::lmom_gengamma(1:5, scale = 2, shape1 = -1, shape2 = 0.8)
  expect_true(all(is.na(g)))
})

# ---------------------------------------------------------------------------
# Special cases: Weibull (k=1, i.e. shape1==shape2) and Gamma (d=1, shape2==1)
# ---------------------------------------------------------------------------
test_that("lmom_gengamma reproduces the Weibull and Gamma special cases", {
  # Weibull k=1: lambda1 = scale*Gamma(1+1/d), lambda2 = scale*Gamma(1+1/d)*(1-2^(-1/d))
  d <- 1.7; sc <- 3
  v <- anyFit:::lmom_gengamma(1:5, sc, d, d)
  expect_equal(unname(v[1]), sc * gamma(1 + 1/d), tolerance = 1e-9)
  expect_equal(unname(v[2]), sc * gamma(1 + 1/d) * (1 - 2^(-1/d)), tolerance = 1e-6)
  # Gamma d=1: lambda1 = k*scale, lambda2 = scale*Gamma(k+0.5)/(sqrt(pi)*Gamma(k))
  k <- 2.3; sc <- 5
  v <- anyFit:::lmom_gengamma(1:5, sc, k, 1)
  expect_equal(unname(v[1]), k * sc, tolerance = 1e-9)
  expect_equal(unname(v[2]), sc * exp(lgamma(k + 0.5) - lgamma(k)) / sqrt(pi), tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# Numerical correctness vs the independent qgengamma-integration oracle
# ---------------------------------------------------------------------------
test_that("lmom_gengamma matches the independent qgengamma-integration oracle (1e-5)", {
  skip_if_not_installed("VGAM")
  for (p in list(c(2, 1.5, 0.8), c(10, 3, 2), c(1, 2, 1.3), c(5, 0.8, 0.6))) {
    cf  <- as.numeric(anyFit:::lmom_gengamma(1:5, p[1], p[2], p[3]))
    orc <- gg_oracle(p[1], p[2], p[3])
    expect_equal(cf, orc, tolerance = 1e-5,
                 info = sprintf("gengamma (%.1f,%.1f,%.2f)", p[1], p[2], p[3]))
  }
})

# ---------------------------------------------------------------------------
# fitlm_gengamma is wired to lmom_gengamma: TheorLMom equals the engine at the
# fitted params, and the fitted params lie in the optim box.
# ---------------------------------------------------------------------------
test_that("fitlm_gengamma is wired to lmom_gengamma and stays in the optim box", {
  skip_if_not_installed("VGAM")
  set.seed(1)
  vals <- anyFit::rgengamma(2000, scale = 8, shape1 = 3, shape2 = 1.4)
  res  <- fitlm_gengamma(make_pos_ts(vals))
  p    <- res$Param
  expect_true(p$scale >= 0.001 && p$scale <= 500)
  expect_true(p$shape1 >= 0.005 && p$shape1 <= 2000)
  expect_true(p$shape2 >= 0.05 && p$shape2 <= 200)
  expect_true(all(is.finite(as.numeric(res$TheorLMom))))
  expect_equal(as.numeric(res$TheorLMom),
               as.numeric(anyFit:::lmom_gengamma(1:5, p$scale, p$shape1, p$shape2)),
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Recovery: a GenGamma sample is fit closely (the two-step seed avoids degenerate
# basins; the quad engine replaces the slow lmom::lmrq).
# ---------------------------------------------------------------------------
test_that("fitlm_gengamma recovers a GenGamma sample with a close fit", {
  skip_if_not_installed("VGAM")
  set.seed(42)
  vals <- anyFit::rgengamma(2000, scale = 8, shape1 = 3, shape2 = 1.4)
  res  <- fitlm_gengamma(make_pos_ts(vals))
  expect_lt(res$GoF$KS, 0.1)
  expect_true(res$Param$shape1 < 1000)   # no degenerate-basin blow-up
})

test_that("fitlm_gengamma_loc fits with a location shift (matrix-indexing safe)", {
  skip_if_not_installed("VGAM")
  set.seed(7)
  vals <- anyFit::rgengamma(2000, scale = 8, shape1 = 3, shape2 = 1.4) + 2
  res  <- fitlm_gengamma_loc(make_pos_ts(vals), location = 1.5)
  expect_equal(res$Param$location, 1.5)
  expect_true(all(is.finite(as.numeric(res$TheorLMom))))
  expect_lt(res$GoF$KS, 0.1)
})
