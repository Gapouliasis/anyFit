library(testthat)

# ===========================================================================
# Tests for fitlm_monthly_nc: monthly L-moment fits on gridded data.
#
# Structure/correctness tests run everywhere. Parallelism tests need real
# future workers that can library(anyFit); they skip_on_cran and skip when
# that is not possible (e.g. an uninstalled / stale package build), reusing
# the guard pattern from test_parallel.R.
# ===========================================================================

# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

# 2x2 regular grid (4 cells), ~2 years of daily dates so all 12 months are
# present with enough positive observations per month to fit norm/gamma.
make_monthly_grid_sxts <- function() {
  set.seed(70)
  dates  <- seq(as.POSIXct("2018-01-01", tz = "UTC"),
                as.POSIXct("2019-12-31", tz = "UTC"), by = "1 day")
  n_time <- length(dates)
  n_pts  <- 4
  mat    <- matrix(abs(rnorm(n_time * n_pts, mean = 8, sd = 2)) + 0.5,
                   nrow = n_time, ncol = n_pts)
  coords <- data.frame(x = c(1.0, 2.0, 1.0, 2.0),
                       y = c(10.0, 10.0, 11.0, 11.0))
  sxts(mat, order.by = dates, coords = coords,
       projection = "+proj=longlat +datum=WGS84")
}

# Same grid but covering only months 1-3, for the missing-month case.
make_partial_months_sxts <- function() {
  set.seed(71)
  dates  <- seq(as.POSIXct("2018-01-01", tz = "UTC"),
                as.POSIXct("2018-03-31", tz = "UTC"), by = "1 day")
  n_time <- length(dates)
  n_pts  <- 4
  mat    <- matrix(abs(rnorm(n_time * n_pts, mean = 8, sd = 2)) + 0.5,
                   nrow = n_time, ncol = n_pts)
  coords <- data.frame(x = c(1.0, 2.0, 1.0, 2.0),
                       y = c(10.0, 10.0, 11.0, 11.0))
  sxts(mat, order.by = dates, coords = coords,
       projection = "+proj=longlat +datum=WGS84")
}

# Numeric signature over every month / candidate / raster component, so
# serial and parallel runs can be compared for exact equality.
sig_rasters <- function(res) {
  comps <- c("raster_params", "raster_TheorLMom", "raster_DataLMom", "raster_GoF")
  unlist(lapply(res, function(month) {
    lapply(month$fit_results, function(cand) {
      lapply(comps, function(k) as.numeric(raster::values(cand[[k]])))
    })
  }))
}

# Probe (copied from test_parallel.R): can future workers load anyFit?
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

cands <- c("norm", "gamma")

# ---------------------------------------------------------------------------
# Structure / correctness tests (always run)
# ---------------------------------------------------------------------------

test_that("returns a list named by month.name, one per month present", {
  skip_if_not_installed("raster")
  res <- fitlm_monthly_nc(data = make_monthly_grid_sxts(), candidates = "norm")
  expect_type(res, "list")
  expect_length(res, 12L)
  expect_identical(names(res), month.name)
})

test_that("each month mirrors fitlm_nc output", {
  skip_if_not_installed("raster")
  res <- fitlm_monthly_nc(data = make_monthly_grid_sxts(), candidates = "norm")
  expect_named(res[["January"]], c("fit_results", "gof_plots"), ignore.order = TRUE)
})

test_that("one fit_results entry per candidate, full-grid rasters", {
  skip_if_not_installed("raster")
  res <- fitlm_monthly_nc(data = make_monthly_grid_sxts(), candidates = cands)
  jan <- res[["January"]]
  expect_equal(names(jan$fit_results), cands)
  params <- jan$fit_results[["norm"]]$raster_params
  expect_true(inherits(params, "Raster"))
  expect_equal(raster::ncell(params), 4L)
})

test_that("missing months yield exactly the months present", {
  skip_if_not_installed("raster")
  res <- fitlm_monthly_nc(data = make_partial_months_sxts(), candidates = "norm")
  expect_length(res, 3L)
  expect_identical(names(res), month.name[1:3])
})

test_that("ignore_zeros keeps grid extent with NA at dropped cells", {
  skip_if_not_installed("raster")
  obj <- make_monthly_grid_sxts()
  obj[, c(1, 3)] <- 0          # empty the x = 1 column (cells 1 and 3)

  res    <- fitlm_monthly_nc(data = obj, candidates = "norm", ignore_zeros = TRUE)
  params <- res[["January"]]$fit_results[["norm"]]$raster_params

  expect_equal(raster::ncell(params), 4L)            # full 2x2 grid preserved
  pts <- raster::rasterToPoints(params)              # non-NA cells only
  expect_equal(nrow(pts), 2L)                        # only the x = 2 column fitted
  expect_true(all(pts[, "x"] == 2))                  # dropped x = 1 cells are NA
})

test_that("parallel_by rejects an unknown value", {
  skip_if_not_installed("raster")
  expect_error(
    fitlm_monthly_nc(data = make_partial_months_sxts(), candidates = "norm",
                     parallel_by = "bogus"),
    "should be one of"
  )
})

# ---------------------------------------------------------------------------
# Parallelism tests
# ---------------------------------------------------------------------------

test_that("cells axis == serial (raster values identical)", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  g   <- make_monthly_grid_sxts()
  ser <- fitlm_monthly_nc(data = g, candidates = cands, ignore_zeros = TRUE)
  par <- fitlm_monthly_nc(data = g, candidates = cands, ignore_zeros = TRUE,
                          parallel = TRUE, parallel_by = "cells")
  expect_equal(sig_rasters(par), sig_rasters(ser))
})

test_that("months axis == serial, and == cells axis", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  g     <- make_monthly_grid_sxts()
  ser   <- fitlm_monthly_nc(data = g, candidates = cands, ignore_zeros = TRUE)
  cells <- fitlm_monthly_nc(data = g, candidates = cands, ignore_zeros = TRUE,
                            parallel = TRUE, parallel_by = "cells")
  mons  <- fitlm_monthly_nc(data = g, candidates = cands, ignore_zeros = TRUE,
                            parallel = TRUE, parallel_by = "months")
  expect_equal(sig_rasters(mons), sig_rasters(ser))
  expect_equal(sig_rasters(mons), sig_rasters(cells))
})

test_that("months axis: user-set plan is respected and left intact", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  before <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(before), add = TRUE)
  cur <- future::plan()
  invisible(fitlm_monthly_nc(data = make_monthly_grid_sxts(), candidates = "norm",
                             parallel = TRUE, parallel_by = "months"))
  expect_identical(class(future::plan()), class(cur))
})

test_that("months axis: result is invariant to worker count", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  g   <- make_monthly_grid_sxts()
  old <- future::plan(future::multisession, workers = 2)
  r2  <- sig_rasters(fitlm_monthly_nc(data = g, candidates = "norm",
                                      parallel = TRUE, parallel_by = "months"))
  future::plan(future::multisession, workers = 3)
  r3  <- sig_rasters(fitlm_monthly_nc(data = g, candidates = "norm",
                                      parallel = TRUE, parallel_by = "months"))
  future::plan(old)
  expect_equal(r2, r3)
})

test_that("months axis: deterministic and emits no RNG warning", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  g <- make_monthly_grid_sxts()
  expect_no_warning(
    r1 <- fitlm_monthly_nc(data = g, candidates = "norm",
                           parallel = TRUE, parallel_by = "months"))
  r2 <- fitlm_monthly_nc(data = g, candidates = "norm",
                         parallel = TRUE, parallel_by = "months")
  expect_equal(sig_rasters(r1), sig_rasters(r2))
})

test_that("months axis: restores plan and globals option (on.exit handlers)", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  # Start from a default sequential plan; the call sets and must restore it.
  before_opt <- getOption("future.globals.maxSize")
  invisible(fitlm_monthly_nc(data = make_monthly_grid_sxts(), candidates = "norm",
                             parallel = TRUE, parallel_by = "months"))
  expect_true(inherits(future::plan(), "sequential"))
  expect_identical(getOption("future.globals.maxSize"), before_opt)
})
