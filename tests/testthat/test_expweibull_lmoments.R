library(testthat)

# ===========================================================================
# tanh-sinh L-moment engine lmom_expweibull (R/lmom_closedform.R) and its use
# inside fitlm_expweibull (two-step seed over ExpWeibull_InitValues).
# ===========================================================================

# Independent oracle: beta_r = integral_0^1 q(F) F^r dF via the package quantile,
# then PWM -> L-moments. A different computation path from the tanh-sinh engine.
ew_pwm_to_lmom <- function(beta) {
  rmax <- length(beta); lambda <- numeric(rmax)
  for (r in 1:rmax) { j <- 0:(r - 1)
    lambda[r] <- sum((-1)^(r - 1 - j) * choose(r - 1, j) * choose(r - 1 + j, j) * beta[j + 1]) }
  c(lambda[1], lambda[2], lambda[3] / lambda[2], lambda[4] / lambda[2], lambda[5] / lambda[2])
}
ew_oracle <- function(scale, shape1, shape2, rmax = 5) {
  beta <- vapply(0:(rmax - 1), function(r)
    integrate(function(F) anyFit::qexpweibull(F, scale, shape1, shape2) * F^r, 0, 1,
              rel.tol = 1e-10, subdivisions = 4000L)$value, numeric(1))
  ew_pwm_to_lmom(beta)
}

make_pos_ts <- function(vals) {
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq_along(vals) * 86400
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "X")), order.by = dates)
}

# ---------------------------------------------------------------------------
# Engine shape / layout
# ---------------------------------------------------------------------------
test_that("lmom_expweibull returns a length-5 named vector for order 1:5", {
  g <- anyFit:::lmom_expweibull(1:5, scale = 2, shape1 = 1.5, shape2 = 0.8)
  expect_length(g, 5L)
  expect_equal(names(g), c("lambda_1", "lambda_2", "tau_3", "tau_4", "tau_5"))
  expect_true(all(is.finite(g)))
})

# ---------------------------------------------------------------------------
# Special cases: exponential (a=b=1) and ordinary Weibull (b=1)
# ---------------------------------------------------------------------------
test_that("lmom_expweibull reproduces the exponential special case", {
  # exponential a=b=1: (lambda1, lambda2, tau3, tau4, tau5) = (s, s/2, 1/3, 1/6, 1/10)
  s <- 3
  v <- as.numeric(anyFit:::lmom_expweibull(1:5, s, 1, 1, h = 1/64, X = 5))
  expect_equal(v, c(s, s/2, 1/3, 1/6, 1/10), tolerance = 1e-8)
})

test_that("lmom_expweibull reproduces the Weibull (b=1) special case", {
  # Weibull b=1: lambda1 = s*Gamma(1+1/a), lambda2 = s*Gamma(1+1/a)*(1-2^(-1/a))
  a <- 2; s <- 10
  v <- anyFit:::lmom_expweibull(1:2, s, a, 1, h = 1/64, X = 5)
  expect_equal(unname(v[1]), s * gamma(1 + 1/a), tolerance = 1e-9)
  expect_equal(unname(v[2]), s * gamma(1 + 1/a) * (1 - 2^(-1/a)), tolerance = 1e-9)
  # lambda1..tau4 vs lmom::lmrwei (exact Weibull). NB lmrwei tau5 is sign-wrong, so only 1:4 here.
  e4  <- as.numeric(anyFit:::lmom_expweibull(1:4, s, a, 1, h = 1/64, X = 5))
  ex4 <- as.numeric(lmom::lmrwei(c(0, s, a), nmom = 4))
  expect_equal(e4, ex4, tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# Numerical correctness vs the independent qexpweibull-integration oracle
# (also covers tau5, which lmom::lmrwei gets sign-wrong).
# ---------------------------------------------------------------------------
test_that("lmom_expweibull matches the independent qexpweibull-integration oracle (1e-6)", {
  for (p in list(c(2, 1.5, 0.8), c(10, 0.8, 2), c(1, 2, 1.3), c(5, 0.6, 0.5))) {
    cf  <- as.numeric(anyFit:::lmom_expweibull(1:5, p[1], p[2], p[3], h = 1/64, X = 5))
    orc <- ew_oracle(p[1], p[2], p[3])
    expect_equal(cf, orc, tolerance = 1e-6,
                 info = sprintf("expweibull (%.1f,%.1f,%.2f)", p[1], p[2], p[3]))
  }
})

# ---------------------------------------------------------------------------
# fitlm_expweibull is wired to lmom_expweibull: TheorLMom equals the engine at the
# fitted params, and the fitted params lie in the optim box.
# ---------------------------------------------------------------------------
test_that("fitlm_expweibull is wired to lmom_expweibull and stays in the optim box", {
  set.seed(1)
  vals <- anyFit::rexpweibull(2000, scale = 8, shape1 = 1.5, shape2 = 1.2)
  res  <- fitlm_expweibull(make_pos_ts(vals))
  p    <- res$Param
  expect_true(p$scale >= 0.001 && p$scale <= 1e6)
  expect_true(p$shape1 >= 0.05 && p$shape1 <= 50)
  expect_true(p$shape2 >= 0.02 && p$shape2 <= 50)
  expect_true(all(is.finite(as.numeric(res$TheorLMom))))
  expect_equal(as.numeric(res$TheorLMom),
               as.numeric(anyFit:::lmom_expweibull(1:5, p$scale, p$shape1, p$shape2)),
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Recovery: an ExpWeibull sample is fit closely without a shape blow-up
# (the two-step seed + tanh-sinh engine replace the slow lmom::pelq path).
# ---------------------------------------------------------------------------
test_that("fitlm_expweibull recovers an ExpWeibull sample without a shape blow-up", {
  set.seed(42)
  vals <- anyFit::rexpweibull(2000, scale = 10, shape1 = 2, shape2 = 0.5)
  res  <- fitlm_expweibull(make_pos_ts(vals))
  expect_true(res$Param$shape1 < 50)
  expect_true(res$Param$shape2 < 50)
  expect_lt(res$GoF$KS, 0.1)
})
