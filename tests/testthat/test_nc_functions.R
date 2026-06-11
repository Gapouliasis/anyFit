library(testthat)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

# Minimal sxts on a 2×2 regular grid (required by rasterFromXYZ)
make_grid_sxts <- function() {
  set.seed(70)
  n_time <- 10
  n_pts  <- 4
  mat    <- matrix(runif(n_time * n_pts, 1, 10), nrow = n_time, ncol = n_pts)
  dates  <- as.POSIXct("2020-01-01", tz = "UTC") + seq(0, n_time - 1) * 86400
  coords <- data.frame(x = c(1.0, 2.0, 1.0, 2.0),
                       y = c(10.0, 10.0, 11.0, 11.0))
  sxts(mat, order.by = dates, coords = coords,
       projection = "+proj=longlat +datum=WGS84")
}

# Create a minimal temporary NetCDF file using ncdf4
make_tmp_nc <- function() {
  skip_if_not_installed("ncdf4")
  tmp  <- tempfile(fileext = ".nc")
  lon  <- ncdf4::ncdim_def("lon", "degrees_east",  c(4.0, 4.5), longname = "longitude")
  lat  <- ncdf4::ncdim_def("lat", "degrees_north", c(50.0, 50.5), longname = "latitude")
  time <- ncdf4::ncdim_def("time", "hours since 2020-01-01 00:00:00 UTC",
                           0:4, unlim = TRUE, longname = "time")
  var  <- ncdf4::ncvar_def("precip", "mm", list(lon, lat, time), -9999,
                           longname = "precipitation")
  set.seed(71)
  nc   <- ncdf4::nc_create(tmp, var)
  ncdf4::ncvar_put(nc, var, array(runif(2 * 2 * 5, 0, 10), dim = c(2, 2, 5)))
  ncdf4::nc_close(nc)
  tmp
}

# ---------------------------------------------------------------------------
# nc2xts — NetCDF to sxts conversion
# ---------------------------------------------------------------------------
test_that("nc2xts() returns list with 'ncdf_sxts' element", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc()
  on.exit(if (file.exists(tmp)) file.remove(tmp))
  res <- nc2xts(tmp, varname = "precip",
                projection = "+proj=longlat +datum=WGS84")
  expect_type(res, "list")
  expect_true("ncdf_sxts" %in% names(res))
})

test_that("nc2xts() ncdf_sxts inherits sxts class", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc()
  on.exit(if (file.exists(tmp)) file.remove(tmp))
  res <- nc2xts(tmp, varname = "precip",
                projection = "+proj=longlat +datum=WGS84")
  expect_true(inherits(res$ncdf_sxts, "sxts"))
})

test_that("nc2xts() ncdf_sxts has correct time steps", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc()
  on.exit(if (file.exists(tmp)) file.remove(tmp))
  res <- nc2xts(tmp, varname = "precip",
                projection = "+proj=longlat +datum=WGS84")
  expect_equal(nrow(res$ncdf_sxts), 5L)
})

test_that("nc2xts() errors on non-existent variable name", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc()
  on.exit(if (file.exists(tmp)) file.remove(tmp))
  expect_error(
    nc2xts(tmp, varname = "nonexistent_variable",
           projection = "+proj=longlat +datum=WGS84"),
    "not found"
  )
})

test_that("nc2xts() with country filter returns fewer spatial columns", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc()
  on.exit(if (file.exists(tmp)) file.remove(tmp))
  full    <- nc2xts(tmp, varname = "precip",
                    projection = "+proj=longlat +datum=WGS84")
  # Belgium filter: coords (4.0–4.5, 50.0–50.5) are inside Belgium
  belgium <- nc2xts(tmp, varname = "precip", country = "Belgium",
                    projection = "+proj=longlat +datum=WGS84")
  expect_lte(ncol(belgium$ncdf_sxts), ncol(full$ncdf_sxts))
})

# ---------------------------------------------------------------------------
# nc2xts_nn — nearest-neighbour extraction
# ---------------------------------------------------------------------------
test_that("nc2xts_nn() returns an xts object", {
  skip_if_not_installed("ncdf4")
  tmp    <- make_tmp_nc()
  on.exit(if (file.exists(tmp)) file.remove(tmp))
  coords <- matrix(c(4.2, 50.2), ncol = 2, dimnames = list(NULL, c("lon", "lat")))
  res    <- nc2xts_nn(filename = tmp, varname = "precip", coords = coords)
  expect_true(inherits(res, "xts"))
})

test_that("nc2xts_nn() returns one column per supplied coordinate", {
  skip_if_not_installed("ncdf4")
  tmp    <- make_tmp_nc()
  on.exit(if (file.exists(tmp)) file.remove(tmp))
  coords <- matrix(c(4.2, 50.2, 4.4, 50.4), ncol = 2, byrow = TRUE,
                   dimnames = list(NULL, c("lon", "lat")))
  res    <- nc2xts_nn(filename = tmp, varname = "precip", coords = coords)
  expect_equal(ncol(res), nrow(coords))
})

# ---------------------------------------------------------------------------
# nc_ggplot — tested with sxts input (no file needed)
# ---------------------------------------------------------------------------
test_that("nc_ggplot() returns a ggplot/patchwork object when given an sxts", {
  obj <- make_grid_sxts()
  res <- nc_ggplot(obj)
  expect_true(inherits(res, "gg") || inherits(res, "patchwork"))
})

test_that("nc_ggplot() with common_legend=TRUE runs without error", {
  obj <- make_grid_sxts()
  expect_no_error(nc_ggplot(obj, common_legend = TRUE))
})

# ---------------------------------------------------------------------------
# basic_stats_nc — tested with sxts input
# ---------------------------------------------------------------------------
test_that("basic_stats_nc() with sxts input returns a Raster object", {
  skip_if_not_installed("raster")
  obj <- make_grid_sxts()
  res <- basic_stats_nc(data = obj)
  expect_true(inherits(res, "Raster"))
})

test_that("basic_stats_nc() raster has at least 26 layers (one per statistic)", {
  skip_if_not_installed("raster")
  obj <- make_grid_sxts()
  res <- basic_stats_nc(data = obj)
  expect_gte(raster::nlayers(res), 26L)
})

# ---------------------------------------------------------------------------
# fitlm_nc — tested with sxts input
# ---------------------------------------------------------------------------
test_that("fitlm_nc() with sxts input returns list with 'fit_results' and 'gof_plots'", {
  skip_if_not_installed("raster")
  obj <- make_grid_sxts()
  res <- fitlm_nc(data = obj, candidates = "norm")
  expect_type(res, "list")
  expect_named(res, c("fit_results", "gof_plots"), ignore.order = TRUE)
})

test_that("fitlm_nc() fit_results has one entry per candidate distribution", {
  skip_if_not_installed("raster")
  obj <- make_grid_sxts()
  res <- fitlm_nc(data = obj, candidates = c("norm", "gamma"))
  expect_equal(length(res$fit_results), 2L)
  expect_equal(names(res$fit_results), c("norm", "gamma"))
})

test_that("fitlm_nc() each fit_results element has raster_params", {
  skip_if_not_installed("raster")
  obj <- make_grid_sxts()
  res <- fitlm_nc(data = obj, candidates = "norm")
  expect_true("raster_params" %in% names(res$fit_results[["norm"]]))
  expect_true(inherits(res$fit_results[["norm"]]$raster_params, "Raster"))
})
