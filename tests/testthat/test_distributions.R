library(testthat)

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

# ---------------------------------------------------------------------------
# Section A — d/p/q/r consistency checks
# ---------------------------------------------------------------------------

# --- Exponential ---
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

# --- Gamma3 ---
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

# --- Rayleigh ---
test_that("drayleigh() returns non-negative density values", {
  x  <- seq(0.1, 5, by = 0.1)
  fx <- drayleigh(x, location = 0, scale = 1)
  expect_true(all(fx >= 0))
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

# --- GEV ---
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
  # qgev signature: qgev(p, shape, scale, location) — note different arg order
  expect_equal(qgev(pgev(x, 5, 2, 0.1), 0.1, 2, 5), x, tolerance = 1e-5)
})

test_that("rgev() returns n samples", {
  samp <- rgev(100, location = 5, scale = 2, shape = 0.1)
  expect_equal(length(samp), 100L)
})

# --- Gumbel ---
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

# --- Burr ---
# dburr signature: dburr(x, scale, shape1, shape2, PW=1) — no location param
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

# --- Lognormal ---
# dlognorm signature: dlognorm(x, location=0, scale, shape)
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

# ---------------------------------------------------------------------------
# Section B — fitlm_* single-series functions
# ---------------------------------------------------------------------------

test_that("fitlm_norm() returns standard structure with finite params", {
  ts  <- make_norm_ts()
  res <- fitlm_norm(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_norm() runs without error with ignore_zeros=TRUE", {
  ts <- make_zeros_ts()
  expect_no_error(fitlm_norm(ts, ignore_zeros = TRUE))
})

test_that("fitlm_exp() returns standard structure with finite params", {
  ts  <- make_gamma_ts()   # positive values needed for exponential
  res <- fitlm_exp(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_exp() runs without error with ignore_zeros=TRUE", {
  ts <- make_zeros_ts()
  expect_no_error(fitlm_exp(ts, ignore_zeros = TRUE))
})

test_that("fitlm_gamma() returns standard structure with finite params", {
  ts  <- make_gamma_ts()
  res <- fitlm_gamma(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_gamma3() returns standard structure with finite params", {
  ts  <- make_gamma_ts()
  res <- fitlm_gamma3(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_gev() returns standard structure with finite params", {
  ts  <- make_norm_ts()
  res <- fitlm_gev(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_lognorm() returns standard structure with finite params", {
  ts  <- make_gamma_ts()
  res <- fitlm_lognorm(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_weibull() returns standard structure with finite params", {
  ts  <- make_gamma_ts()
  res <- fitlm_weibull(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_rayleigh() returns standard structure with finite params", {
  ts  <- make_gamma_ts()
  res <- fitlm_rayleigh(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_burr() returns standard structure with finite params", {
  ts  <- make_gamma_ts()
  res <- fitlm_burr(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_dagum() returns standard structure with finite params", {
  ts  <- make_gamma_ts()
  res <- fitlm_dagum(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_gengamma() returns standard structure with finite params", {
  ts  <- make_gamma_ts()
  res <- fitlm_gengamma(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_expweibull() returns standard structure with finite params", {
  ts  <- make_gamma_ts()
  res <- fitlm_expweibull(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_GPD() returns standard structure with finite params", {
  ts  <- make_gamma_ts()
  res <- fitlm_GPD(ts)
  check_fitlm_structure(res)
})

test_that("fitlm_genlogi() returns standard structure with finite params", {
  ts  <- make_norm_ts()
  res <- fitlm_genlogi(ts)
  check_fitlm_structure(res)
})

# ---------------------------------------------------------------------------
# Section C — multi-distribution helpers
# ---------------------------------------------------------------------------

test_that("fitlm_multi() with diagnostic_plots=TRUE returns 5 list elements", {
  ts  <- make_norm_ts()
  res <- fitlm_multi(ts, candidates = list("norm", "gamma"))
  expect_type(res, "list")
  expect_true(all(c("parameter_list", "GoF_summary", "diagnostics",
                    "QQplot", "PPplot") %in% names(res)))
})

test_that("fitlm_multi() GoF_summary has one column per candidate", {
  ts  <- make_norm_ts()
  res <- fitlm_multi(ts, candidates = list("norm", "gamma"))
  expect_equal(ncol(res$GoF_summary), 2L)
  expect_equal(colnames(res$GoF_summary), c("norm", "gamma"))
})

test_that("fitlm_multi() with diagnostic_plots=FALSE returns only params and GoF", {
  ts  <- make_norm_ts()
  res <- fitlm_multi(ts, candidates = list("norm"), diagnostic_plots = FALSE)
  expect_true("parameter_list" %in% names(res))
  expect_true("GoF_summary"    %in% names(res))
  expect_false("QQplot" %in% names(res))
})

test_that("fitlm_nxts() returns list with params, diagnostic_plots, QQ_plots, PP_plots", {
  ts  <- make_multi_ts()
  res <- fitlm_nxts(ts, candidates = list("norm"), diagnostic_plots = TRUE)
  expect_true(all(c("params", "diagnostic_plots", "QQ_plots", "PP_plots") %in% names(res)))
})

test_that("fitlm_nxts() params has one entry per column of ts", {
  ts  <- make_multi_ts()
  res <- fitlm_nxts(ts, candidates = list("norm"), diagnostic_plots = FALSE)
  expect_equal(length(res$params), ncol(ts))
  expect_equal(names(res$params), colnames(ts))
})

test_that("fitlm_monthly() returns params_monthly, GoF_monthly, monthly_QQplot, monthly_PPplot", {
  ts  <- make_norm_ts()
  res <- fitlm_monthly(ts, candidates = list("norm"))
  expect_true(all(c("params_monthly", "GoF_monthly",
                    "monthly_QQplot", "monthly_PPplot") %in% names(res)))
})

test_that("fitlm_monthly() params_monthly is a list with one entry per candidate", {
  ts  <- make_norm_ts()
  res <- fitlm_monthly(ts, candidates = list("norm", "gamma"))
  expect_equal(length(res$params_monthly), 2L)
  expect_equal(names(res$params_monthly), c("norm", "gamma"))
})

test_that("fitlm_monthly() params_monthly columns are month names", {
  ts  <- make_norm_ts()
  res <- fitlm_monthly(ts, candidates = list("norm"))
  # At least some months should appear as column names
  expect_true(any(colnames(res$params_monthly[[1]]) %in% month.name))
})

# ---------------------------------------------------------------------------
# Section D — diagnostic helpers
# ---------------------------------------------------------------------------

test_that("fit_diagnostics() returns Diagnostic_Plots, GoF, QQplot, PPplot", {
  ts     <- make_norm_ts()
  fit    <- fitlm_norm(ts)
  res    <- fit_diagnostics(ts, dist = "norm", params = fit$Param)
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

test_that("fit_diagnostics() QQplot and PPplot are ggplot objects", {
  ts  <- make_norm_ts()
  fit <- fitlm_norm(ts)
  res <- fit_diagnostics(ts, dist = "norm", params = fit$Param)
  expect_true(inherits(res$QQplot, "gg"))
  expect_true(inherits(res$PPplot, "gg"))
})

test_that("LRatio_check() returns list with 'distributions' and 'multi_plots'", {
  ts      <- make_multi_ts()
  lmoms   <- lmom_stats(ts)
  res     <- LRatio_check(lmoms)
  expect_type(res, "list")
  expect_named(res, c("distributions", "multi_plots"), ignore.order = TRUE)
})

test_that("LRatio_check() distributions contains entries for the four tested distributions", {
  ts    <- make_multi_ts()
  lm    <- lmom_stats(ts)
  res   <- LRatio_check(lm)
  dist_names <- names(res$distributions)
  expect_true(all(c("Dagum", "GGamma", "ExpWeibull", "BurrXII") %in% dist_names))
})
