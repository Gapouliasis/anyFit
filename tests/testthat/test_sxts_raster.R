library(testthat)

# ---------------------------------------------------------------------------
# Shared fixture — 5 time steps x 4 spatial points on a 2×2 regular grid
# (must be a regular grid for rasterFromXYZ to work)
# ---------------------------------------------------------------------------
make_obj <- function() {
  set.seed(60)
  n_time <- 5
  n_pts  <- 4
  mat    <- matrix(runif(n_time * n_pts, 0, 10), nrow = n_time, ncol = n_pts)
  dates  <- as.POSIXct("2020-01-01", tz = "UTC") + 0:(n_time - 1) * 3600
  coords <- data.frame(x = c(1.0, 2.0, 1.0, 2.0),
                       y = c(10.0, 10.0, 11.0, 11.0))
  sxts(mat, order.by = dates,
       coords = coords, projection = "+proj=longlat +datum=WGS84")
}

# ---------------------------------------------------------------------------
# rasterFromSxts
# ---------------------------------------------------------------------------
test_that("rasterFromSxts() returns a Raster object", {
  obj <- make_obj()
  res <- rasterFromSxts(obj)
  expect_true(inherits(res, "Raster"))
})

test_that("rasterFromSxts() RasterBrick has one layer per time step", {
  obj <- make_obj()
  res <- rasterFromSxts(obj)
  expect_equal(raster::nlayers(res), nrow(obj))
})

test_that("rasterFromSxts() preserves projection string", {
  obj  <- make_obj()
  res  <- rasterFromSxts(obj)
  proj <- raster::projection(res)
  expect_true(grepl("longlat", proj, ignore.case = TRUE))
})

# ---------------------------------------------------------------------------
# sxtsFromRaster — round-trip
# ---------------------------------------------------------------------------
test_that("sxtsFromRaster() returns an sxts object", {
  obj     <- make_obj()
  ras     <- rasterFromSxts(obj)
  back    <- sxtsFromRaster(ras)
  expect_true(inherits(back, "sxts"))
})

test_that("sxtsFromRaster() round-trip preserves number of spatial points", {
  obj  <- make_obj()
  ras  <- rasterFromSxts(obj)
  back <- sxtsFromRaster(ras)
  expect_equal(ncol(back), ncol(obj))
})

test_that("sxtsFromRaster() round-trip preserves number of time steps", {
  obj  <- make_obj()
  ras  <- rasterFromSxts(obj)
  back <- sxtsFromRaster(ras)
  expect_equal(nrow(back), nrow(obj))
})

# ---------------------------------------------------------------------------
# mask.sxts — error when no mask provided
# ---------------------------------------------------------------------------
test_that("mask.sxts() errors when all mask arguments are NA", {
  obj <- make_obj()
  expect_error(
    mask.sxts(obj, xlim = NA, ylim = NA, shapefile_name = NA, mask = NA),
    "Provide xlim"
  )
})

# ---------------------------------------------------------------------------
# zonal_stats — continent argument
# ---------------------------------------------------------------------------
test_that("zonal_stats() with continent='Europe' returns xts with country name column", {
  # Four points clearly inside continental Europe (western France area)
  set.seed(61)
  n_time <- 3
  n_pts  <- 4
  mat    <- matrix(runif(n_time * n_pts, 0, 10), nrow = n_time, ncol = n_pts)
  dates  <- as.POSIXct("2020-01-01", tz = "UTC") + 0:(n_time - 1) * 3600
  coords <- data.frame(x = c(2.0, 2.5, 3.0, 3.5),
                       y = c(47.0, 47.5, 48.0, 48.5))
  obj    <- sxts(mat, order.by = dates, coords = coords,
                 projection = "+proj=longlat +datum=WGS84")

  result <- zonal_stats(obj, continent = "Europe")
  expect_true(inherits(result, "xts"))
  expect_equal(ncol(result), 1L)
  expect_equal(colnames(result), "France")
  expect_equal(nrow(result), n_time)
})

# ---------------------------------------------------------------------------
# zonal_stats — country argument (already tested in test_sxts.R, verify dims)
# ---------------------------------------------------------------------------
test_that("zonal_stats() country result has correct time dimensions", {
  set.seed(62)
  n_time <- 4
  n_pts  <- 3
  mat    <- matrix(runif(n_time * n_pts, 0, 10), nrow = n_time, ncol = n_pts)
  dates  <- as.POSIXct("2020-01-01", tz = "UTC") + 0:(n_time - 1) * 3600
  coords <- data.frame(x = c(4.3, 4.5, 4.0),
                       y = c(50.8, 50.9, 50.6))
  obj    <- sxts(mat, order.by = dates, coords = coords,
                 projection = "+proj=longlat +datum=WGS84")

  result <- zonal_stats(obj, country = "Belgium")
  expect_equal(nrow(result), n_time)
})
