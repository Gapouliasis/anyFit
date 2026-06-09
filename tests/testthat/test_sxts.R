library(testthat)

# ---------------------------------------------------------------------------
# Shared fixture — 5 time steps x 4 spatial points, all in-memory
# ---------------------------------------------------------------------------
make_obj <- function() {
  set.seed(42)
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
# Constructor
# ---------------------------------------------------------------------------
test_that("sxts() creates object with correct class and attributes", {
  obj <- make_obj()
  expect_s3_class(obj, "sxts")
  expect_s3_class(obj, "xts")
  expect_equal(nrow(obj), 5L)
  expect_equal(ncol(obj), 4L)
  expect_equal(attr(obj, "elements"), 4L)
  expect_equal(attr(obj, "projection"), "+proj=longlat +datum=WGS84")
  expect_equal(nrow(attr(obj, "coords")), 4L)
})

test_that("sxts() errors when coord rows != data columns", {
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + 0:2 * 3600
  mat   <- matrix(1:12, nrow = 3, ncol = 4)
  bad_coords <- data.frame(x = 1:3, y = 1:3)   # 3 rows, but 4 cols
  expect_error(sxts(mat, order.by = dates, coords = bad_coords),
               "Number of coordinate rows")
})

test_that("sxts() errors when coords has only one column", {
  dates <- as.POSIXct("2020-01-01", tz = "UTC") + 0:2 * 3600
  mat   <- matrix(1:6, nrow = 3, ncol = 2)
  bad_coords <- data.frame(x = 1:2)
  expect_error(sxts(mat, order.by = dates, coords = bad_coords),
               "coords must have at least 2 columns")
})

# ---------------------------------------------------------------------------
# is.sxts
# ---------------------------------------------------------------------------
test_that("is.sxts() returns TRUE for sxts and FALSE for plain xts", {
  obj      <- make_obj()
  plain    <- xts::xts(matrix(1:6, 3), as.POSIXct("2020-01-01") + 0:2 * 3600)
  expect_true(is.sxts(obj))
  expect_false(is.sxts(plain))
  expect_false(is.sxts(42))
})

# ---------------------------------------------------------------------------
# Accessor: coords
# ---------------------------------------------------------------------------
test_that("coords() returns data.frame with x and y columns", {
  obj <- make_obj()
  co  <- coords(obj)
  expect_s3_class(co, "data.frame")
  expect_equal(nrow(co), 4L)
  expect_true(all(c("x", "y") %in% names(co)))
  expect_equal(co$x, c(1, 2, 1, 2))
})

# ---------------------------------------------------------------------------
# Accessor: projection
# ---------------------------------------------------------------------------
test_that("projection() returns the projection string", {
  obj <- make_obj()
  expect_equal(projection(obj), "+proj=longlat +datum=WGS84")
})

# ---------------------------------------------------------------------------
# Accessor: elements
# ---------------------------------------------------------------------------
test_that("elements() returns number of spatial locations", {
  obj <- make_obj()
  expect_equal(elements(obj), 4L)
})

# ---------------------------------------------------------------------------
# attributes.sxts
# ---------------------------------------------------------------------------
test_that("attributes.sxts() returns list with elements, coords, projection", {
  obj   <- make_obj()
  attrs <- attributes.sxts(obj)
  expect_type(attrs, "list")
  expect_named(attrs, c("elements", "coords", "projection"),
               ignore.order = TRUE)
  expect_equal(attrs$elements, 4L)
  expect_equal(attrs$projection, "+proj=longlat +datum=WGS84")
})

# ---------------------------------------------------------------------------
# print.sxts
# ---------------------------------------------------------------------------
test_that("print.sxts() runs without error and mentions 'sxts'", {
  obj <- make_obj()
  out <- capture.output(print(obj))
  expect_true(any(grepl("sxts", out, ignore.case = TRUE)))
})

# ---------------------------------------------------------------------------
# str.sxts  (was broken: used 'obj' instead of 'object' — now fixed)
# ---------------------------------------------------------------------------
test_that("str.sxts() runs without error", {
  obj <- make_obj()
  expect_no_error(capture.output(str(obj)))
})

# ---------------------------------------------------------------------------
# summary.sxts
# ---------------------------------------------------------------------------
test_that("summary.sxts() runs without error and reports element count", {
  obj <- make_obj()
  out <- capture.output(summary(obj))
  expect_true(any(grepl("Number of elements", out)))
})

# ---------------------------------------------------------------------------
# [.sxts — time subsetting
# ---------------------------------------------------------------------------
test_that("[.sxts time subsetting preserves class and coords", {
  obj <- make_obj()
  sub <- obj[1:3, ]
  expect_s3_class(sub, "sxts")
  expect_equal(nrow(sub), 3L)
  expect_equal(ncol(sub), 4L)
  expect_equal(coords(sub), coords(obj))   # coords unchanged
})

# ---------------------------------------------------------------------------
# [.sxts — spatial (column) subsetting
# ---------------------------------------------------------------------------
test_that("[.sxts spatial subsetting updates coords and elements", {
  obj <- make_obj()
  sub <- obj[, 1:2]
  expect_s3_class(sub, "sxts")
  expect_equal(ncol(sub), 2L)
  expect_equal(nrow(coords(sub)), 2L)
  expect_equal(elements(sub), 2L)
  expect_equal(coords(sub)$x, c(1, 2))
})

# ---------------------------------------------------------------------------
# lag.sxts  (needed restore_sxts to be defined)
# ---------------------------------------------------------------------------
test_that("lag.sxts() returns sxts with same dimensions", {
  obj    <- make_obj()
  lagged <- lag(obj, k = 1)
  expect_s3_class(lagged, "sxts")
  expect_equal(nrow(lagged), nrow(obj))
  expect_equal(ncol(lagged), ncol(obj))
  expect_equal(projection(lagged), projection(obj))
})

# ---------------------------------------------------------------------------
# diff.sxts  (needed restore_sxts to be defined)
# ---------------------------------------------------------------------------
test_that("diff.sxts() returns sxts with same dims and NA in first row", {
  obj  <- make_obj()
  dobj <- diff(obj)
  expect_s3_class(dobj, "sxts")
  # xts diff keeps row count; first row is NA
  expect_equal(nrow(dobj), nrow(obj))
  expect_equal(ncol(dobj), ncol(obj))
  expect_true(all(is.na(dobj[1, ])))
  expect_equal(projection(dobj), projection(obj))
})

# ---------------------------------------------------------------------------
# Ops.sxts — arithmetic  (needed restore_sxts to be defined)
# ---------------------------------------------------------------------------
test_that("Ops.sxts arithmetic preserves class and scales values", {
  obj    <- make_obj()
  scaled <- obj * 2
  expect_s3_class(scaled, "sxts")
  expect_equal(as.numeric(scaled), as.numeric(obj) * 2)
  expect_equal(projection(scaled), projection(obj))
})

# ---------------------------------------------------------------------------
# Ops.sxts — comparison
# ---------------------------------------------------------------------------
test_that("Ops.sxts comparison returns sxts of logicals", {
  obj  <- make_obj()
  mask <- obj > 5
  expect_s3_class(mask, "sxts")
  expect_type(as.numeric(mask), "double")   # xts stores logicals as 0/1
  expect_equal(nrow(mask), nrow(obj))
})

# ---------------------------------------------------------------------------
# mask.sxts — bounding-box masking
# ---------------------------------------------------------------------------
test_that("mask.sxts() bounding box retains only points inside limits", {
  obj    <- make_obj()
  # xlim [0.5, 1.5], ylim [9.5, 10.5] captures only (x=1, y=10) → column 1
  masked <- mask.sxts(obj, xlim = c(0.5, 1.5), ylim = c(9.5, 10.5))
  expect_s3_class(masked, "sxts")
  expect_equal(ncol(masked), 1L)
  expect_equal(coords(masked)$x, 1)
  expect_equal(coords(masked)$y, 10)
})

# ---------------------------------------------------------------------------
# zonal_stats — country masking using built-in world_data dataset
# ---------------------------------------------------------------------------
test_that("zonal_stats() with country returns xts with country as column name", {
  # Three points clearly inside Belgium (Brussels area)
  mat2    <- matrix(runif(3 * 3, 0, 10), nrow = 3)
  dates2  <- as.POSIXct("2020-01-01", tz = "UTC") + 0:2 * 3600
  coords2 <- data.frame(x = c(4.3, 4.5, 4.0),
                        y = c(50.8, 50.9, 50.6))
  obj2    <- sxts(mat2, order.by = dates2, coords = coords2,
                  projection = "+proj=longlat +datum=WGS84")

  result <- zonal_stats(obj2, country = "Belgium")
  expect_s3_class(result, "xts")
  expect_equal(ncol(result), 1L)
  expect_equal(colnames(result), "Belgium")
  expect_equal(nrow(result), 3L)
})
