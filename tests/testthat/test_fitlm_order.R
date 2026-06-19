library(testthat)

# ===========================================================================
# Tests for the per-candidate `order` argument threaded through the fitlm
# wrapper functions (fitlm_multi -> fitlm_nxts -> fitlm_nc -> fitlm_monthly_nc,
# and fitlm_monthly).
#
# `order` is a named list mapping a candidate to the L-moment orders matched by
# its optimiser, e.g. list(gengamma = 1:5). Only the numerically-fitted
# distributions (gengamma, gengamma_loc, burr, dagum, expweibull) accept it.
#
# The optimiser-based fitters carry NO RNG in the fit (the seed is data-derived
# and deterministic), so the strongest assertion is equivalence: a wrapper call
# with a given `order` must produce the SAME parameters as calling the
# underlying fitlm_<dist> directly with that same `order`. gengamma is the probe
# distribution (its default order is 1:5).
# ===========================================================================

# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

make_gengamma_ts <- function() {
  set.seed(1234)
  vals  <- anyFit::rgengamma(400, scale = 3, shape1 = 2, shape2 = 1.5)
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, 399) * 86400
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "X")), order.by = dates)
}

# Multi-column xts grid (positive gengamma draws) for fitlm_nxts.
make_ggamma_grid_xts <- function() {
  set.seed(1235)
  n_time <- 400
  n_col  <- 3
  mat    <- matrix(anyFit::rgengamma(n_time * n_col, scale = 3, shape1 = 2, shape2 = 1.5),
                   nrow = n_time, ncol = n_col,
                   dimnames = list(NULL, c("A", "B", "C")))
  dates  <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, n_time - 1) * 86400
  xts::xts(mat, order.by = dates)
}

# 2x2 regular grid sxts, ~2 years daily so all 12 months are present, for the
# nc / monthly_nc wrappers (pattern from test_fitlm_monthly_nc.R).
make_ggamma_grid_sxts <- function() {
  set.seed(1236)
  dates  <- seq(as.POSIXct("2018-01-01", tz = "UTC"),
                as.POSIXct("2019-12-31", tz = "UTC"), by = "1 day")
  n_time <- length(dates)
  n_pts  <- 4
  mat    <- matrix(anyFit::rgengamma(n_time * n_pts, scale = 3, shape1 = 2, shape2 = 1.5),
                   nrow = n_time, ncol = n_pts)
  coords <- data.frame(x = c(1.0, 2.0, 1.0, 2.0),
                       y = c(10.0, 10.0, 11.0, 11.0))
  sxts(mat, order.by = dates, coords = coords,
       projection = "+proj=longlat +datum=WGS84")
}

# Probe (from test_parallel.R): can future workers load anyFit?
future_parallel_ok <- function() {
  isTRUE(tryCatch({
    old <- future::plan(future::multisession, workers = 2)
    on.exit(future::plan(old), add = TRUE)
    res <- future.apply::future_lapply(1:2, function(i)
      exists("fitlm_multi", where = asNamespace("anyFit")),
      future.packages = "anyFit")
    all(unlist(res))
  }, error = function(e) FALSE))
}

skip_unless_parallel <- function() {
  skip_on_cran()
  if (!future_parallel_ok())
    skip("future multisession workers cannot load anyFit (install the package)")
}

ord  <- c(2, 4)          # discontinuous order to match
ord3 <- c(2, 3, 4)       # contiguous order to match

# ---------------------------------------------------------------------------
# Group 1 - fitlm_multi (the dispatch point)
# ---------------------------------------------------------------------------

test_that("order has an observable effect on the fitted parameters", {
  x <- make_gengamma_ts()
  default_fit <- fitlm_multi(x, list("gengamma"))$parameter_list$gengamma$Param
  custom_fit  <- fitlm_multi(x, list("gengamma"),
                             order = list(gengamma = ord3))$parameter_list$gengamma$Param
  expect_false(isTRUE(all.equal(default_fit, custom_fit)))
})

test_that("order matches a direct fitlm_gengamma call (discontinuous order)", {
  x <- make_gengamma_ts()
  via_multi <- fitlm_multi(x, list("gengamma"),
                           order = list(gengamma = ord))$parameter_list$gengamma$Param
  direct    <- fitlm_gengamma(x, order = ord)$Param
  expect_equal(via_multi, direct)
})

test_that("order = NULL is a pure no-op (defaults preserved)", {
  x <- make_gengamma_ts()
  via_default  <- fitlm_multi(x, list("gengamma"))$parameter_list$gengamma$Param
  direct       <- fitlm_gengamma(x)$Param
  via_explicit <- fitlm_multi(x, list("gengamma"),
                              order = list(gengamma = 1:5))$parameter_list$gengamma$Param
  expect_equal(via_default, direct)
  expect_equal(via_default, via_explicit)   # 1:5 IS the gengamma default
})

test_that("mixed candidates with a partial order list", {
  x   <- make_gengamma_ts()
  res <- fitlm_multi(x, list("gengamma", "norm"), order = list(gengamma = ord))
  # norm is not in the order list -> untouched
  expect_equal(res$parameter_list$norm$Param, fitlm_norm(x)$Param)
  # gengamma uses the requested order
  expect_equal(res$parameter_list$gengamma$Param, fitlm_gengamma(x, order = ord)$Param)
})

test_that("listing a non-order candidate is silently ignored", {
  x <- make_gengamma_ts()
  expect_silent(
    res <- fitlm_multi(x, list("norm"), order = list(norm = ord))
  )
  expect_equal(res$parameter_list$norm$Param, fitlm_norm(x)$Param)
})

# ---------------------------------------------------------------------------
# Group 2 - fitlm_monthly
# ---------------------------------------------------------------------------

test_that("fitlm_monthly threads order into each month's fit", {
  ts  <- make_ggamma_grid_sxts()[, 1]          # one cell, all 12 months
  res <- fitlm_monthly(ts, list("gengamma"), order = list(gengamma = ord))

  jan_slice <- ts[lubridate::month(ts) == 1]
  direct    <- unname(unlist(fitlm_gengamma(jan_slice, order = ord)$Param))
  via_month <- unname(as.numeric(res$params_monthly$gengamma[, "January"]))
  expect_equal(via_month, direct)
})

# ---------------------------------------------------------------------------
# Group 3 - fitlm_nxts (covers the .make_col_worker factory)
# ---------------------------------------------------------------------------

test_that("fitlm_nxts matches per-column direct fits with order set", {
  grid <- make_ggamma_grid_xts()
  res  <- fitlm_nxts(grid, list("gengamma"), order = list(gengamma = ord),
                     diagnostic_plots = FALSE)
  for (j in seq_len(ncol(grid))) {
    direct <- fitlm_gengamma(grid[, j], order = ord)$Param
    via_nxts <- res$params[[colnames(grid)[j]]]$params$gengamma$Param
    expect_equal(via_nxts, direct)
  }
})

test_that("fitlm_nxts: parallel == serial with order set", {
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  grid <- make_ggamma_grid_xts()
  ser  <- fitlm_nxts(grid, list("gengamma"), order = list(gengamma = ord),
                     diagnostic_plots = FALSE)
  par  <- fitlm_nxts(grid, list("gengamma"), order = list(gengamma = ord),
                     diagnostic_plots = FALSE, parallel = TRUE)
  sig <- function(r) unlist(lapply(r$params, function(v) v$params$gengamma$Param))
  expect_equal(sig(par), sig(ser))
})

# ---------------------------------------------------------------------------
# Group 4 - fitlm_nc
# ---------------------------------------------------------------------------

test_that("fitlm_nc threads order through to the rasters", {
  skip_if_not_installed("raster")
  g <- make_ggamma_grid_sxts()
  default_run <- fitlm_nc(data = g, candidates = "gengamma")
  custom_run  <- fitlm_nc(data = g, candidates = "gengamma",
                          order = list(gengamma = ord))
  explicit5   <- fitlm_nc(data = g, candidates = "gengamma",
                          order = list(gengamma = 1:5))

  v <- function(res) raster::values(res$fit_results$gengamma$raster_params)
  expect_false(isTRUE(all.equal(v(default_run), v(custom_run))))  # effect reaches raster
  expect_equal(v(default_run), v(explicit5))                      # 1:5 == default no-op
})

# ---------------------------------------------------------------------------
# Group 5 - fitlm_monthly_nc
# ---------------------------------------------------------------------------

test_that("fitlm_monthly_nc threads order to the deepest wrapper", {
  skip_if_not_installed("raster")
  g   <- make_ggamma_grid_sxts()
  res <- fitlm_monthly_nc(data = g, candidates = "gengamma",
                          order = list(gengamma = ord))

  jan_slice <- g[lubridate::month(zoo::index(g)) == 1, ]
  direct    <- fitlm_nc(data = jan_slice, candidates = "gengamma",
                        order = list(gengamma = ord))

  v <- function(r) raster::values(r$fit_results$gengamma$raster_params)
  expect_equal(v(res[["January"]]), v(direct))
})
