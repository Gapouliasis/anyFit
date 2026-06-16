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

# NetCDF on a 1/24-degree grid (like NClimGrid-Daily) with single-precision (float32)
# axis values, so the cell spacing wobbles at ~1e-5 the way the real files do.
make_tmp_nc_fine <- function() {
  skip_if_not_installed("ncdf4")
  f32 <- function(v) { con <- rawConnection(raw(0), "wb"); writeBin(as.vector(v), con, size = 4)
                       r <- rawConnectionValue(con); close(con); readBin(r, "numeric", length(v), size = 4) }
  res  <- 1 / 24
  tmp  <- tempfile(fileext = ".nc")
  lon  <- ncdf4::ncdim_def("lon", "degrees_east",  f32(seq(-100, by = res, length.out = 10)))
  lat  <- ncdf4::ncdim_def("lat", "degrees_north", f32(seq(  40, by = res, length.out =  8)))
  time <- ncdf4::ncdim_def("time", "hours since 2020-01-01 00:00:00 UTC", 0:29, unlim = TRUE)
  var  <- ncdf4::ncvar_def("precip", "mm", list(lon, lat, time), -9999)
  set.seed(73)
  nc   <- ncdf4::nc_create(tmp, var)
  ncdf4::ncvar_put(nc, var, array(runif(10 * 8 * 30, 1, 50), dim = c(10, 8, 30)))
  ncdf4::nc_close(nc)
  tmp
}

# ---------------------------------------------------------------------------
# nc2xts — NetCDF to sxts conversion
# ---------------------------------------------------------------------------
test_that("nc2xts() returns an sxts object", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc()
  on.exit(if (file.exists(tmp)) file.remove(tmp))
  res <- nc2xts(tmp, varname = "precip",
                projection = "+proj=longlat +datum=WGS84")
  expect_true(inherits(res, "sxts"))
})

test_that("nc2xts() inherits sxts class", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc()
  on.exit(if (file.exists(tmp)) file.remove(tmp))
  res <- nc2xts(tmp, varname = "precip",
                projection = "+proj=longlat +datum=WGS84")
  expect_true(inherits(res, "sxts"))
})

test_that("nc2xts() has correct time steps", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc()
  on.exit(if (file.exists(tmp)) file.remove(tmp))
  res <- nc2xts(tmp, varname = "precip",
                projection = "+proj=longlat +datum=WGS84")
  expect_equal(nrow(res), 5L)
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
  expect_lte(ncol(belgium), ncol(full))
})

# ---------------------------------------------------------------------------
# nc2xts — bounding-box prefilter (slab read) correctness
# ---------------------------------------------------------------------------

# A larger regular grid so a sub-window is a strict subset of the full grid.
make_tmp_nc_grid <- function() {
  skip_if_not_installed("ncdf4")
  tmp  <- tempfile(fileext = ".nc")
  lon  <- ncdf4::ncdim_def("lon", "degrees_east",  seq(0, 9), longname = "longitude")
  lat  <- ncdf4::ncdim_def("lat", "degrees_north", seq(40, 49), longname = "latitude")
  time <- ncdf4::ncdim_def("time", "hours since 2020-01-01 00:00:00 UTC",
                           0:3, unlim = TRUE)
  var  <- ncdf4::ncvar_def("precip", "mm", list(lon, lat, time), -9999)
  set.seed(72)
  nc   <- ncdf4::nc_create(tmp, var)
  ncdf4::ncvar_put(nc, var, array(runif(10 * 10 * 4, 0, 10), dim = c(10, 10, 4)))
  ncdf4::nc_close(nc)
  tmp
}

test_that("nc2xts() xlim/ylim slab read matches full read + bbox mask", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc_grid()
  on.exit(if (file.exists(tmp)) file.remove(tmp))

  full <- nc2xts(tmp, varname = "precip",
                 projection = "+proj=longlat +datum=WGS84")
  ref  <- mask.sxts(full, xlim = c(2, 5), ylim = c(42, 45))

  slab <- nc2xts(tmp, varname = "precip", xlim = c(2, 5), ylim = c(42, 45),
                 projection = "+proj=longlat +datum=WGS84")

  # Same locations selected (order-independent) and same values
  ref_coords  <- coords(ref)[order(coords(ref)$y,  coords(ref)$x), ]
  slab_coords <- coords(slab)[order(coords(slab)$y, coords(slab)$x), ]
  expect_equal(ref_coords$x, slab_coords$x)
  expect_equal(ref_coords$y, slab_coords$y)
  expect_equal(unname(sum(zoo::coredata(slab))),
               unname(sum(zoo::coredata(ref))))
})

test_that("nc2xts() country slab read matches full read + country mask", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc_grid()  # spans lon 0-9, lat 40-49: covers part of Europe
  on.exit(if (file.exists(tmp)) file.remove(tmp))

  full <- nc2xts(tmp, varname = "precip",
                 projection = "+proj=longlat +datum=WGS84")
  ref  <- mask.sxts(full, mask = world_data[world_data$name == "France", ])

  slab <- nc2xts(tmp, varname = "precip", country = "France",
                 projection = "+proj=longlat +datum=WGS84")

  expect_equal(ncol(slab), ncol(ref))
  expect_equal(unname(sum(zoo::coredata(slab))),
               unname(sum(zoo::coredata(ref))))
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

test_that("nc2xts_nn() matches raster::extract on a multi-cell grid", {
  skip_if_not_installed("ncdf4")
  skip_if_not_installed("raster")
  tmp <- make_tmp_nc_grid()  # lon 0-9, lat 40-49, 4 time steps
  on.exit(if (file.exists(tmp)) file.remove(tmp))

  pts <- matrix(c(2.2, 42.3, 7.4, 47.8), ncol = 2, byrow = TRUE)
  res <- nc2xts_nn(filename = tmp, varname = "precip", coords = pts)

  brick <- raster::brick(tmp, varname = "precip")
  ref   <- t(raster::extract(brick, pts, method = "simple"))

  expect_equal(unname(zoo::coredata(res)), unname(ref))
})

test_that("nc2xts_nn() returns an all-NA column for a point outside the grid", {
  skip_if_not_installed("ncdf4")
  tmp <- make_tmp_nc_grid()
  on.exit(if (file.exists(tmp)) file.remove(tmp))

  # First point inside the grid, second well outside (lon/lat range is 0-9 / 40-49)
  pts <- matrix(c(3.0, 44.0, 100.0, 100.0), ncol = 2, byrow = TRUE)
  res <- nc2xts_nn(filename = tmp, varname = "precip", coords = pts)

  expect_equal(ncol(res), 2L)
  expect_false(all(is.na(zoo::coredata(res)[, 1])))
  expect_true(all(is.na(zoo::coredata(res)[, 2])))
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

test_that("nc_ggplot() accepts a timeseries Raster* object", {
  skip_if_not_installed("raster")
  r   <- rasterFromSxts(make_grid_sxts())
  res <- nc_ggplot(r)
  expect_true(inherits(res, "gg") || inherits(res, "patchwork"))
})

test_that("nc_ggplot() accepts a non-timeseries stats raster", {
  skip_if_not_installed("raster")
  stats_raster <- basic_stats_nc(make_grid_sxts())
  expect_no_error(res <- nc_ggplot(stats_raster, title = TRUE))
  expect_true(inherits(res, "gg") || inherits(res, "patchwork"))
})

test_that("nc_ggplot() errors on unsupported input", {
  expect_error(nc_ggplot("some/path.nc"), "must be an sxts object")
  expect_error(nc_ggplot(data.frame(a = 1)), "must be an sxts object")
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

test_that("fitlm_nc() ignore_zeros keeps the full grid extent with NA at dropped cells", {
  skip_if_not_installed("raster")
  obj <- make_grid_sxts()
  obj[, c(1, 3)] <- 0          # empty the entire x = 1 column (cells 1 and 3)

  res    <- fitlm_nc(data = obj, candidates = "norm", ignore_zeros = TRUE)
  params <- res$fit_results[["norm"]]$raster_params

  expect_true(inherits(params, "Raster"))
  expect_equal(raster::ncell(params), 4L)            # full 2x2 grid preserved
  pts <- raster::rasterToPoints(params)              # non-NA cells only
  expect_equal(nrow(pts), 2L)                        # only the x = 2 column fitted
  expect_true(all(pts[, "x"] == 2))                  # dropped x = 1 cells are NA
})

test_that("fitlm_nc() works on a float-precision (1/24-deg) NetCDF grid via nc2xts", {
  skip_if_not_installed("ncdf4")
  skip_if_not_installed("raster")
  tmp <- make_tmp_nc_fine()
  on.exit(if (file.exists(tmp)) file.remove(tmp))

  obj    <- nc2xts(tmp, varname = "precip", projection = "+proj=longlat +datum=WGS84")
  res    <- fitlm_nc(data = obj, candidates = "norm", ignore_zeros = TRUE)
  params <- res$fit_results[["norm"]]$raster_params

  expect_true(inherits(params, "Raster"))
  expect_equal(raster::ncell(params), 10L * 8L)   # full grid rasterized
})
