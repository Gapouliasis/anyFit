library(testthat)

# ===========================================================================
# Closed-form L-moment engines lmom_burr / lmom_dagum (R/lmom_closedform.R)
# and their use inside fitlm_burr / fitlm_dagum.
# ===========================================================================

# Independent oracle: beta_r = integral_0^1 q(F) F^r dF via the package quantile,
# then PWM -> L-moments. A different computation path from the closed Beta form.
pwm_to_lmom <- function(beta) {
  rmax <- length(beta); lambda <- numeric(rmax)
  for (r in 1:rmax) { j <- 0:(r - 1)
    lambda[r] <- sum((-1)^(r - 1 - j) * choose(r - 1, j) * choose(r - 1 + j, j) * beta[j + 1]) }
  c(lambda[1], lambda[2], lambda[3] / lambda[2], lambda[4] / lambda[2], lambda[5] / lambda[2])
}
oracle_lmom <- function(qfun, scale, shape1, shape2, rmax = 5) {
  beta <- vapply(0:(rmax - 1), function(r)
    integrate(function(F) qfun(F, scale, shape1, shape2) * F^r, 0, 1,
              rel.tol = 1e-10, subdivisions = 4000L)$value, numeric(1))
  pwm_to_lmom(beta)
}

make_pos_ts <- function(seed = 7, shape = 2, scale = 4, n = 600) {
  set.seed(seed)
  vals  <- rgamma(n, shape = shape, scale = scale)
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq_len(n) * 86400
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "X")), order.by = dates)
}

# ---------------------------------------------------------------------------
# Engine shape / layout
# ---------------------------------------------------------------------------
test_that("lmom_burr/lmom_dagum return a length-5 named vector for order 1:5", {
  b <- anyFit:::lmom_burr(1:5, scale = 2, shape1 = 1.5, shape2 = 0.5)
  d <- anyFit:::lmom_dagum(1:5, scale = 2, shape1 = 1.5, shape2 = 0.3)
  expect_length(b, 5L); expect_length(d, 5L)
  expect_equal(names(b), c("lambda_1", "lambda_2", "tau_3", "tau_4", "tau_5"))
  expect_equal(names(d), c("lambda_1", "lambda_2", "tau_3", "tau_4", "tau_5"))
  expect_true(all(is.finite(b))); expect_true(all(is.finite(d)))
})

# ---------------------------------------------------------------------------
# Numerical correctness vs the independent quantile-integration oracle
# ---------------------------------------------------------------------------
test_that("lmom_burr matches the independent qburr-integration oracle (1e-6)", {
  for (p in list(c(2, 1.5, 0.5), c(5, 2, 0.3), c(1, 1, 0.4))) {
    cf  <- as.numeric(anyFit:::lmom_burr(1:5, p[1], p[2], p[3]))
    orc <- oracle_lmom(anyFit::qburr, p[1], p[2], p[3])
    expect_equal(cf, orc, tolerance = 1e-6,
                 info = sprintf("burr (%.1f,%.1f,%.2f)", p[1], p[2], p[3]))
  }
})

test_that("lmom_dagum matches the independent qdagum-integration oracle (1e-6)", {
  for (p in list(c(2, 1.5, 0.3), c(10, 3, 0.2), c(1, 0.8, 0.4))) {
    cf  <- as.numeric(anyFit:::lmom_dagum(1:5, p[1], p[2], p[3]))
    orc <- oracle_lmom(anyFit::qdagum, p[1], p[2], p[3])
    expect_equal(cf, orc, tolerance = 1e-6,
                 info = sprintf("dagum (%.1f,%.1f,%.2f)", p[1], p[2], p[3]))
  }
})

# ---------------------------------------------------------------------------
# Domain guard: mean (hence every L-moment) diverges for shape2 >= 1
# ---------------------------------------------------------------------------
test_that("engines return Inf when shape2 >= 1", {
  expect_true(all(is.infinite(anyFit:::lmom_burr(1:5,  2, 1.5, 1.0))))
  expect_true(all(is.infinite(anyFit:::lmom_burr(1:5,  2, 1.5, 1.2))))
  expect_true(all(is.infinite(anyFit:::lmom_dagum(1:5, 2, 1.5, 1.0))))
  expect_true(all(is.infinite(anyFit:::lmom_dagum(1:5, 2, 1.5, 1.5))))
})

# ---------------------------------------------------------------------------
# fitlm_* use the closed engine: TheorLMom consistency, params in the optim box,
# and the fit actually matches the sample L-moments
# ---------------------------------------------------------------------------
# The wiring check: fitlm_*$TheorLMom is produced by the closed engine (not lmrp), so it
# equals lmom_*(1:5, fitted params) exactly, and the fitted params lie in the optim box.
# (Fit *quality* is model-dependent and is covered by the end-to-end gate, not here.)
test_that("fitlm_burr is wired to lmom_burr and stays in the optim box", {
  res <- fitlm_burr(make_pos_ts())
  p   <- res$Param
  expect_true(p$shape2 > 0 && p$shape2 <= 1)
  expect_true(p$shape1 >= 0.5 && p$shape1 <= 50)
  expect_true(all(is.finite(as.numeric(res$TheorLMom))))
  expect_equal(as.numeric(res$TheorLMom),
               as.numeric(anyFit:::lmom_burr(1:5, p$scale, p$shape1, p$shape2)),
               tolerance = 1e-10)
})

test_that("fitlm_dagum is wired to lmom_dagum and stays in the optim box", {
  res <- fitlm_dagum(make_pos_ts())
  p   <- res$Param
  expect_true(p$shape2 > 0 && p$shape2 <= 0.5)
  expect_true(all(is.finite(as.numeric(res$TheorLMom))))
  expect_equal(as.numeric(res$TheorLMom),
               as.numeric(anyFit:::lmom_dagum(1:5, p$scale, p$shape1, p$shape2)),
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Regression guard for the min-max-NN start (R/Distributions_Anyfit.R): on data
# drawn from a Dagum inside its L-space, the fit must recover a sensible shape1
# (the old scale-free-ratio start trapped optim at shape1 >> 100) and fit well.
# ---------------------------------------------------------------------------
test_that("fitlm_dagum recovers a Dagum sample without the shape1 blow-up", {
  set.seed(42)
  vals  <- anyFit::rdagum(2000, scale = 10, shape1 = 2, shape2 = 0.3)
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq_along(vals) * 86400
  ts    <- xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "X")), order.by = dates)
  res   <- fitlm_dagum(ts)
  expect_true(res$Param$shape1 < 50)     # bad start gave shape1 >> 100
  expect_lt(res$GoF$KS, 0.1)             # close fit to its own family
})
