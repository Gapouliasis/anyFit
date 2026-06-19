library(testthat)

# ===========================================================================
# Tests for monthly_stats_nc: monthly basic statistics on gridded data.
#
# Structure/correctness tests run everywhere. Parallelism tests need real
# future workers that can library(anyFit); they skip_on_cran and skip when
# that is not possible, reusing the guard pattern from test_parallel.R.
# ===========================================================================

# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

# 2x2 regular grid (4 cells), ~2 years of daily dates so all 12 months are
# present with enough positive observations per month for the statistics.
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

# Numeric signature over every month's stat raster, so serial and parallel
# runs can be compared for exact equality.
sig_rasters <- function(res) {
  unlist(lapply(res, function(r) as.numeric(raster::values(r))))
}

# Probe (copied from test_parallel.R): can future workers load anyFit?
future_parallel_ok <- function() {
  isTRUE(tryCatch({
    old <- future::plan(future::multisession, workers = 2)
    on.exit(future::plan(old), add = TRUE)
    res <- future.apply::future_lapply(1:2, function(i)
      exists("basic_stats_nc", where = asNamespace("anyFit")),
      future.packages = "anyFit")
    all(unlist(res))
  }, error = function(e) FALSE))
}

skip_unless_parallel <- function() {
  skip_on_cran()
  if (!future_parallel_ok())
    skip("future multisession workers cannot load anyFit (install the package)")
}

# ---------------------------------------------------------------------------
# Structure / correctness tests (always run)
# ---------------------------------------------------------------------------

test_that("returns a list named by month.name, one per month present", {
  skip_if_not_installed("raster")
  res <- monthly_stats_nc(data = make_monthly_grid_sxts())
  expect_type(res, "list")
  expect_length(res, 12L)
  expect_identical(names(res), month.name)
})

test_that("each month is a full-grid stat raster", {
  skip_if_not_installed("raster")
  res <- monthly_stats_nc(data = make_monthly_grid_sxts())
  jan <- res[["January"]]
  expect_true(inherits(jan, "Raster"))
  expect_equal(raster::ncell(jan), 4L)
  expect_true("Mean" %in% names(jan))
})

test_that("ignore_zeros drops the transition layers", {
  skip_if_not_installed("raster")
  g    <- make_monthly_grid_sxts()
  full <- monthly_stats_nc(data = g, ignore_zeros = FALSE)
  noz  <- monthly_stats_nc(data = g, ignore_zeros = TRUE)
  expect_lt(raster::nlayers(noz[["January"]]),
            raster::nlayers(full[["January"]]))
})

test_that("missing months yield exactly the months present", {
  skip_if_not_installed("raster")
  res <- monthly_stats_nc(data = make_partial_months_sxts())
  expect_length(res, 3L)
  expect_identical(names(res), month.name[1:3])
})

test_that("per-month output equals basic_stats_nc on that month's subset", {
  skip_if_not_installed("raster")
  g   <- make_monthly_grid_sxts()
  res <- monthly_stats_nc(data = g, ignore_zeros = TRUE)

  jul_rows <- which(lubridate::month(zoo::index(g)) == 7)
  ref      <- basic_stats_nc(data = g[jul_rows, ], ignore_zeros = TRUE)

  expect_equal(raster::values(res[["July"]]), raster::values(ref))
})

# ---------------------------------------------------------------------------
# Parallelism tests
# ---------------------------------------------------------------------------

test_that("months axis == serial (raster values identical)", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  g   <- make_monthly_grid_sxts()
  ser <- monthly_stats_nc(data = g, ignore_zeros = TRUE)
  par <- monthly_stats_nc(data = g, ignore_zeros = TRUE, parallel = TRUE)
  expect_equal(sig_rasters(par), sig_rasters(ser))
})

test_that("user-set plan is respected and left intact", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  before <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(before), add = TRUE)
  cur <- future::plan()
  invisible(monthly_stats_nc(data = make_monthly_grid_sxts(), parallel = TRUE))
  expect_identical(class(future::plan()), class(cur))
})

test_that("result is invariant to worker count", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  g   <- make_monthly_grid_sxts()
  old <- future::plan(future::multisession, workers = 2)
  r2  <- sig_rasters(monthly_stats_nc(data = g, parallel = TRUE))
  future::plan(future::multisession, workers = 3)
  r3  <- sig_rasters(monthly_stats_nc(data = g, parallel = TRUE))
  future::plan(old)
  expect_equal(r2, r3)
})

test_that("parallel is deterministic with no parallel-RNG warning", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  g <- make_monthly_grid_sxts()

  # basic_stats_nc emits benign min/max-on-NA warnings (from rasterFromXYZ); only
  # the parallel-RNG warning, guarded by future.seed = TRUE, must never appear.
  warns <- character(0)
  r1 <- withCallingHandlers(
    monthly_stats_nc(data = g, parallel = TRUE),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("random|seed|RNG|UNRELIABLE", warns, ignore.case = TRUE)))

  r2 <- suppressWarnings(monthly_stats_nc(data = g, parallel = TRUE))
  expect_equal(sig_rasters(r1), sig_rasters(r2))
})

test_that("restores plan and globals option (on.exit handlers)", {
  skip_if_not_installed("raster")
  skip_unless_parallel()
  before_opt <- getOption("future.globals.maxSize")
  invisible(monthly_stats_nc(data = make_monthly_grid_sxts(), parallel = TRUE))
  expect_true(inherits(future::plan(), "sequential"))
  expect_identical(getOption("future.globals.maxSize"), before_opt)
})
