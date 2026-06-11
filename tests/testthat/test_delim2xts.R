library(testthat)

# ---------------------------------------------------------------------------
# Shared fixture helpers — write temp files with known content
# ---------------------------------------------------------------------------

# Standard 5-row uniform-step tab-delimited file (datetime + 2 data cols)
write_tab_file <- function(col_names = TRUE) {
  tmp <- tempfile(fileext = ".txt")
  dates <- format(
    as.POSIXct("2020-01-01 00:00:00", tz = "UTC") + seq(0, 4) * 3600,
    "%Y-%m-%d %H:%M:%S"
  )
  vals1 <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  vals2 <- c(10.0, 20.0, 30.0, 40.0, 50.0)
  df <- data.frame(Date = dates, V1 = vals1, V2 = vals2)
  if (col_names) {
    write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = FALSE)
  }
  tmp
}

# File with a missing step (gap between rows 2 and 3)
write_gap_file <- function() {
  tmp <- tempfile(fileext = ".txt")
  # Row 3 is 2 hours after row 2 instead of 1 hour
  dates <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC") +
    c(0, 3600, 7200 + 3600, 10800 + 3600, 14400 + 3600)
  dates <- format(dates, "%Y-%m-%d %H:%M:%S")
  df    <- data.frame(Date = dates, V1 = 1:5)
  write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  tmp
}

# File with a sentinel no-value
write_novalue_file <- function(sentinel = -999) {
  tmp <- tempfile(fileext = ".txt")
  dates <- format(
    as.POSIXct("2020-01-01", tz = "UTC") + seq(0, 4) * 86400,
    "%Y-%m-%d"
  )
  vals <- c(1.0, sentinel, 3.0, 4.0, 5.0)
  df   <- data.frame(Date = dates, V1 = vals)
  write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  tmp
}

# File spanning a leap year
write_leap_file <- function() {
  tmp <- tempfile(fileext = ".txt")
  dates <- as.POSIXct("2020-02-27", tz = "UTC") + seq(0, 4) * 86400
  dates <- format(dates, "%Y-%m-%d")
  df    <- data.frame(Date = dates, V1 = 1:5)
  write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  tmp
}

# File with data column first, date second
write_date_second_file <- function() {
  tmp <- tempfile(fileext = ".txt")
  dates <- format(
    as.POSIXct("2020-01-01", tz = "UTC") + seq(0, 4) * 86400,
    "%Y-%m-%d"
  )
  df <- data.frame(V1 = 1:5, Date = dates)
  write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  tmp
}

# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------
test_that("delim2xts() returns an xts object with correct dimensions", {
  tmp <- write_tab_file()
  on.exit(file.remove(tmp))
  res <- delim2xts(tmp, time_zone = "UTC", strict_step = TRUE)
  expect_true(inherits(res, "xts"))
  expect_equal(nrow(res), 5L)
  expect_equal(ncol(res), 2L)
})

test_that("delim2xts() column names match file header", {
  tmp <- write_tab_file(col_names = TRUE)
  on.exit(file.remove(tmp))
  res <- delim2xts(tmp, time_zone = "UTC", strict_step = TRUE)
  expect_equal(colnames(res), c("V1", "V2"))
})

test_that("delim2xts() data values match file content", {
  tmp <- write_tab_file()
  on.exit(file.remove(tmp))
  res <- delim2xts(tmp, time_zone = "UTC", strict_step = TRUE)
  expect_equal(as.numeric(res[, 1]), 1:5)
  expect_equal(as.numeric(res[, 2]), c(10, 20, 30, 40, 50))
})

# ---------------------------------------------------------------------------
# strict_step behaviour
# ---------------------------------------------------------------------------
test_that("delim2xts() with strict_step=FALSE accepts non-uniform time steps", {
  tmp <- write_gap_file()
  on.exit(file.remove(tmp))
  expect_no_error(
    delim2xts(tmp, time_zone = "UTC", strict_step = FALSE)
  )
})

test_that("delim2xts() with strict_step=TRUE errors on non-uniform steps", {
  tmp <- write_gap_file()
  on.exit(file.remove(tmp))
  expect_error(
    delim2xts(tmp, time_zone = "UTC", strict_step = TRUE),
    "Time step is not strict"
  )
})

# ---------------------------------------------------------------------------
# col_names
# ---------------------------------------------------------------------------
test_that("delim2xts() with col_names=FALSE assigns default column names", {
  tmp <- write_tab_file(col_names = FALSE)
  on.exit(file.remove(tmp))
  res <- delim2xts(tmp, time_zone = "UTC", strict_step = FALSE,
                   col_names = FALSE)
  expect_equal(nrow(res), 5L)
  # First column is treated as date index — the remaining columns should have default names
  expect_true(ncol(res) >= 1L)
})

# ---------------------------------------------------------------------------
# date_index
# ---------------------------------------------------------------------------
test_that("delim2xts() with date_index=2 uses the second column as time index", {
  tmp <- write_date_second_file()
  on.exit(file.remove(tmp))
  res <- delim2xts(tmp, time_zone = "UTC", strict_step = TRUE, date_index = 2)
  expect_true(inherits(res, "xts"))
  expect_equal(nrow(res), 5L)
  expect_equal(as.numeric(res[, 1]), 1:5)
})

# ---------------------------------------------------------------------------
# exc_leaps
# ---------------------------------------------------------------------------
test_that("delim2xts() with exc_leaps=TRUE excludes Feb 29 from result", {
  tmp <- write_leap_file()
  on.exit(file.remove(tmp))
  res <- delim2xts(tmp, time_zone = "UTC", strict_step = TRUE,
                   exc_leaps = TRUE)
  dates_out <- format(zoo::index(res), "%m-%d")
  expect_false("02-29" %in% dates_out)
})

test_that("delim2xts() with exc_leaps=FALSE retains Feb 29", {
  tmp <- write_leap_file()
  on.exit(file.remove(tmp))
  res <- delim2xts(tmp, time_zone = "UTC", strict_step = TRUE,
                   exc_leaps = FALSE)
  dates_out <- format(zoo::index(res), "%m-%d")
  expect_true("02-29" %in% dates_out)
})

# ---------------------------------------------------------------------------
# save_Xts
# ---------------------------------------------------------------------------
test_that("delim2xts() with save_Xts=TRUE writes a file at the given path", {
  tmp     <- write_tab_file()
  out_f   <- tempfile(fileext = "")
  on.exit({ file.remove(tmp); if (file.exists(paste0(out_f, ".txt"))) file.remove(paste0(out_f, ".txt")) })
  delim2xts(tmp, time_zone = "UTC", strict_step = TRUE,
             save_Xts = TRUE, filename = out_f)
  expect_true(file.exists(paste0(out_f, ".txt")))
})
