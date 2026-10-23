library(testthat)

# ===========================================================================
# Shared method-of-L-moments solver (R/fitlm_lmom_solver.R) behind fitlm_burr,
# fitlm_dagum, fitlm_gengamma(_loc) and fitlm_expweibull: log-shape L-BFGS-B on
# the scale-free L-ratios, analytic scale, restart on stall, and the fit status
# returned in Res$Diag. Also the log-space qburr / dburr the widened Burr box
# needs.
# ===========================================================================

make_pos_ts <- function(vals) {
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq_along(vals) * 86400
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "X")), order.by = dates)
}

# One sample per family, drawn well inside its attainable L-ratio region.
make_family_samples <- function(n = 20000) {
  set.seed(101)
  list(burr       = list(x = anyFit::rburr(n, scale = 10, shape1 = 2, shape2 = 0.3),
                         par = c(10, 2, 0.3)),
       dagum      = list(x = anyFit::rdagum(n, scale = 10, shape1 = 2, shape2 = 0.3),
                         par = c(10, 2, 0.3)),
       gengamma   = list(x = anyFit::rgengamma(n, scale = 8, shape1 = 3, shape2 = 1.4),
                         par = c(8, 3, 1.4)),
       expweibull = list(x = anyFit::rexpweibull(n, scale = 3, shape1 = 1.5, shape2 = 2),
                         par = c(3, 1.5, 2)))
}

STATUS_LEVELS <- c("matched", "matched_at_bound", "outside_domain", "not_converged")

# ---------------------------------------------------------------------------
# Diag structure
# ---------------------------------------------------------------------------
test_that("the four fitters and fitlm_gengamma_loc return a Diag with the documented fields", {
  skip_if_not_installed("VGAM")
  s <- make_family_samples(2000)
  for (fam in names(s)) {
    res <- get(paste0("fitlm_", fam))(make_pos_ts(s[[fam]]$x))
    expect_named(res$Diag, c("status", "resid", "convergence", "evals", "restarts"))
    expect_true(res$Diag$status %in% STATUS_LEVELS)
    expect_true(is.finite(res$Diag$resid) && res$Diag$resid >= 0)
    expect_true(res$Diag$restarts %in% 0:1)
    expect_true(res$Diag$evals > 0)
  }
  res <- fitlm_gengamma_loc(make_pos_ts(s$gengamma$x + 2), location = 1.5)
  expect_named(res$Diag, c("status", "resid", "convergence", "evals", "restarts"))
})

# ---------------------------------------------------------------------------
# Recovery: an exact (order 1:3) match on a large sample from each family
# ---------------------------------------------------------------------------
test_that("each family is matched to its own large sample (order 1:3)", {
  skip_if_not_installed("VGAM")
  s <- make_family_samples()
  for (fam in names(s)) {
    res <- get(paste0("fitlm_", fam))(make_pos_ts(s[[fam]]$x), order = 1:3)
    expect_equal(res$Diag$status, "matched", info = fam)
    expect_lt(res$Diag$resid, 1e-3)
    # fitted L-ratios reproduce the sample ones
    th <- as.numeric(res$TheorLMom); dl <- as.numeric(unlist(res$DataLMom))
    expect_equal(c(th[2] / th[1], th[3]), c(dl[2] / dl[1], dl[3]), tolerance = 1e-4, info = fam)
    expect_equal(unname(unlist(res$Param)), s[[fam]]$par, tolerance = 0.15, info = fam)
  }
})

test_that("order 1:5 is a least-squares fit whose status follows the optimiser", {
  s <- make_family_samples(5000)
  res <- fitlm_burr(make_pos_ts(s$burr$x), order = 1:5)
  expect_equal(res$Diag$convergence, 0)
  expect_equal(res$Diag$status, "matched")
})

test_that("the exponentiated Weibull reaches shape2 beyond the old ceiling of 50", {
  set.seed(5)
  vals <- anyFit::rexpweibull(20000, scale = 1, shape1 = 4, shape2 = 300)
  res  <- fitlm_expweibull(make_pos_ts(vals))
  expect_equal(res$Diag$status, "matched")
  expect_gt(res$Param$shape2, 50)
  expect_lte(res$Param$shape2, 1e4)
})

test_that("a sample more skewed than the lognormal is outside the generalised gamma", {
  skip_if_not_installed("VGAM")
  set.seed(9)
  # population (L-CV, tau_3) = (0.33, 0.40): above the lognormal edge of the GG family
  vals <- anyFit::rburr(5000, scale = 1, shape1 = 4, shape2 = 0.4)
  res  <- fitlm_gengamma(make_pos_ts(vals), order = 1:3)
  expect_equal(res$Diag$status, "outside_domain")
  expect_gte(res$Diag$resid, 1e-3)
  expect_true(all(is.finite(unlist(res$Param))))
  # the same sample is inside the Burr XII family
  expect_equal(fitlm_burr(make_pos_ts(vals), order = 1:3)$Diag$status, "matched")
})

test_that("ignore_zeros drops the zeros before solving", {
  set.seed(13)
  wet  <- anyFit::rburr(3000, scale = 10, shape1 = 2, shape2 = 0.3)
  vals <- c(wet, rep(0, 1500))[sample(4500)]
  res_zero <- fitlm_burr(make_pos_ts(vals), ignore_zeros = TRUE, order = 1:3)
  res_wet  <- fitlm_burr(make_pos_ts(vals[vals > 0.01]), order = 1:3)
  expect_equal(res_zero$Param, res_wet$Param)
  expect_equal(res_zero$Diag$status, "matched")
})

# ---------------------------------------------------------------------------
# order -> L-ratio mapping
# ---------------------------------------------------------------------------
test_that("order selects L-ratios; lambda_1 always sets the scale", {
  s <- make_family_samples(3000)
  x <- make_pos_ts(s$dagum$x)
  expect_equal(fitlm_dagum(x, order = c(2, 4))$Param, fitlm_dagum(x, order = c(1, 2, 4))$Param)
  expect_false(isTRUE(all.equal(fitlm_dagum(x, order = c(2, 4))$Param,
                                fitlm_dagum(x, order = 1:3)$Param)))
})

test_that("an order with fewer than two L-ratios is an error", {
  s <- make_family_samples(1000)
  x <- make_pos_ts(s$burr$x)
  expect_error(fitlm_burr(x, order = 1:2), "at least two")
  expect_error(fitlm_expweibull(x, order = c(1, 3)), "at least two")
})

# ---------------------------------------------------------------------------
# Log-space qburr / dburr
# ---------------------------------------------------------------------------
test_that("qburr and dburr agree with the direct formulas and stay finite for large shape1", {
  q_direct <- function(p, scale, a, b) scale * ((((1 - p))^(-a * b) - 1) / (a * b))^(1 / a)
  d_direct <- function(x, scale, a, b) a * scale^(-a) * x^(a - 1) * (a * b * (x / scale)^a + 1)^(-1 / (a * b) - 1)
  p <- c(0.001, 0.1, 0.5, 0.9, 0.999)
  expect_equal(qburr(p, 3, 2, 0.4), q_direct(p, 3, 2, 0.4), tolerance = 1e-10)
  x <- c(0.1, 1, 3, 10)
  expect_equal(dburr(x, 3, 2, 0.4), d_direct(x, 3, 2, 0.4), tolerance = 1e-10)
  # shape1 * shape2 = 5000: the direct quantile overflows to Inf, the log form does not
  expect_false(is.finite(q_direct(0.99, 1, 1e4, 0.5)))
  expect_true(all(is.finite(qburr(p, 1, 1e4, 0.5))))
  expect_true(all(is.finite(dburr(c(0.5, 1, 2), 1, 1e4, 0.5))))
  expect_equal(qburr(0, 1, 2, 0.4), 0)
  expect_equal(qburr(c(0.2, 0.5), 1, 2, 0.4, PW = 0.7)[1], 0)   # p <= 1 - PW is a zero
})

# ---------------------------------------------------------------------------
# Attainable-region polygons
# ---------------------------------------------------------------------------
test_that("the *_LSpace polygons are single-ring SpatialPolygonsDataFrames used by LRatio_check", {
  for (nm in c("BurrXII_LSpace", "Dagum_LSpace", "GGamma_LSpace", "ExpWeibull_LSpace")) {
    p <- get(nm, envir = asNamespace("anyFit"))
    expect_s4_class(p, "SpatialPolygonsDataFrame")
    expect_gt(nrow(p@polygons[[1]]@Polygons[[1]]@coords), 100)
  }
  set.seed(42)
  ts <- xts::xts(cbind(a = rgamma(2000, shape = 0.8, scale = 5), b = rgamma(2000, shape = 2, scale = 7)),
                 order.by = as.Date("2000-01-01") + 0:1999)
  lc <- LRatio_check(lmom_stats(ts, ignore_zeros = TRUE))
  expect_named(lc, c("distributions", "multi_plots"))
  expect_equal(length(lc$distributions$ExpWeibull$check), 2L)
})
