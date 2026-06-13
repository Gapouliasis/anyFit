library(testthat)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
make_plain_ts <- function() {
  set.seed(10)
  n     <- 365 * 3
  vals  <- rnorm(n, mean = 5, sd = 2)
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 3600
  xts::xts(matrix(vals, ncol = 1, dimnames = list(NULL, "V1")), order.by = dates)
}

make_zeros_ts <- function() {
  set.seed(11)
  ts  <- make_plain_ts()
  idx <- sample(nrow(ts), size = round(nrow(ts) * 0.20))
  ts[idx, ] <- 0
  ts
}

make_na_ts <- function() {
  set.seed(12)
  ts  <- make_plain_ts()
  idx <- sample(nrow(ts), size = round(nrow(ts) * 0.10))
  ts[idx, ] <- NA
  ts
}

make_multi_ts <- function() {
  set.seed(13)
  n     <- 365 * 3
  vals  <- matrix(rnorm(n * 2, mean = 5, sd = 2), ncol = 2,
                  dimnames = list(NULL, c("A", "B")))
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 3600
  xts::xts(vals, order.by = dates)
}

# ---------------------------------------------------------------------------
# basic_stats — return structure
# ---------------------------------------------------------------------------
test_that("basic_stats() returns list with 'plot' and 'stats_table'", {
  ts  <- make_plain_ts()
  res <- basic_stats(ts)
  expect_type(res, "list")
  expect_named(res, c("plot", "stats_table"), ignore.order = TRUE)
})

test_that("basic_stats() stats_table has 34 rows when ignore_zeros=FALSE", {
  ts  <- make_plain_ts()
  res <- basic_stats(ts)
  expect_equal(nrow(res$stats_table), 34L)
})

test_that("basic_stats() stats_table contains all expected row names", {
  ts       <- make_plain_ts()
  res      <- basic_stats(ts)
  expected <- c("NumofData", "NumofMisData", "PercOfMissingData",
                "Min", "Max", "Mean", "Var", "StDev",
                "Variation", "Mom3", "Skewness", "Kurtosis",
                "Lmean", "LScale", "L3", "L4",
                "LVariation", "LSkewness", "LKurtosis",
                "Pdr", "Q5", "Q25", "Q50", "Q75", "Q95", "IQR",
                "MeanDAfterZero", "VarDAfterZero",
                "MeanDBeforeZero", "VarDBeforeZero",
                "MeanDAfterD", "VarDAfterD", "ProbDD", "ProbNDND")
  expect_true(all(expected %in% rownames(res$stats_table)))
})

test_that("basic_stats() NumofData matches nrow(ts)", {
  ts  <- make_plain_ts()
  res <- basic_stats(ts)
  expect_equal(res$stats_table["NumofData", "Value"], nrow(ts))
})

test_that("basic_stats() Mean is close to mean(ts)", {
  ts  <- make_plain_ts()
  res <- basic_stats(ts)
  expect_equal(res$stats_table["Mean", "Value"],
               round(mean(ts, na.rm = TRUE), 2),
               tolerance = 1e-2)
})

test_that("basic_stats() stats_table has 26 rows when ignore_zeros=TRUE", {
  ts  <- make_zeros_ts()
  res <- basic_stats(ts, ignore_zeros = TRUE)
  expect_equal(nrow(res$stats_table), 26L)
})

test_that("basic_stats() Pdr is positive when series contains zeros", {
  ts  <- make_zeros_ts()
  res <- basic_stats(ts)
  expect_gt(res$stats_table["Pdr", "Value"], 0)
})

test_that("basic_stats() NumofMisData is positive when series contains NAs", {
  ts  <- make_na_ts()
  res <- basic_stats(ts)
  expect_gt(res$stats_table["NumofMisData", "Value"], 0)
})

test_that("basic_stats() pstart/pend restricts the plotted period without error", {
  ts  <- make_plain_ts()
  expect_no_error(
    basic_stats(ts, pstart = "2016-01-01", pend = "2016-12-31")
  )
})

# ---------------------------------------------------------------------------
# lmom_stats — return structure
# ---------------------------------------------------------------------------
test_that("lmom_stats() returns a data.frame with 5 rows", {
  ts  <- make_plain_ts()
  res <- lmom_stats(ts)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5L)
})

test_that("lmom_stats() row names are the five L-moment ratios", {
  ts       <- make_plain_ts()
  res      <- lmom_stats(ts)
  expected <- c("L-Mean", "L-Scale", "L-Skew", "L-Kurtosis", "L-CV")
  expect_true(all(expected %in% rownames(res)))
})

test_that("lmom_stats() has one column per variable in ts", {
  ts  <- make_multi_ts()
  res <- lmom_stats(ts)
  expect_equal(ncol(res), ncol(ts))
  expect_equal(colnames(res), colnames(ts))
})

test_that("lmom_stats() L-Mean is close to the sample mean for Normal data", {
  ts  <- make_plain_ts()
  res <- lmom_stats(ts)
  expect_equal(res["L-Mean", 1], round(mean(ts, na.rm = TRUE), 2),
               tolerance = 0.5)
})

test_that("lmom_stats() runs without error when ignore_zeros=TRUE", {
  ts <- make_zeros_ts()
  expect_no_error(lmom_stats(ts, ignore_zeros = TRUE))
})

test_that("lmom_stats() values are all numeric and finite", {
  ts  <- make_plain_ts()
  res <- lmom_stats(ts)
  expect_true(all(is.numeric(unlist(res))))
  expect_true(all(is.finite(unlist(res))))
})

# ---------------------------------------------------------------------------
# Adversarial parity tests for the vectorised stats engine
# ---------------------------------------------------------------------------
# Compare finiteness pattern, then numeric values where finite. This sidesteps
# NaN-vs-NA pedantry (both are "non-finite") while still pinning exact values.
expect_parity <- function(actual, expected, tol = 1e-8) {
  a <- as.numeric(actual); e <- as.numeric(expected)
  fa <- is.finite(a); fe <- is.finite(e)
  expect_equal(fa, fe)
  if (any(fa)) expect_equal(a[fa], e[fe], tolerance = tol)
}

# Faithful transcription of the ORIGINAL single-series transition logic - the
# only valid oracle for these statistics (no external package computes them).
ref_transition <- function(v, thr = 0.01) {
  v       <- as.numeric(v)
  pos     <- 2:length(v)
  lagpos  <- pos - 1
  Ypos    <- v[pos]
  Ylagpos <- v[lagpos]

  z   <- ifelse(Ypos > thr & Ylagpos <= thr, Ypos, NA);    daz <- z[!is.na(z)]
  z2  <- ifelse(Ypos <= thr & Ylagpos > thr, Ylagpos, NA); dbz <- z2[!is.na(z2)]
  z3  <- ifelse(Ypos > thr & Ylagpos > thr, Ypos, NA);     dad <- z3[!is.na(z3)]
  zdd <- ifelse(Ypos > thr & Ylagpos > thr, Ypos, -99)
  ProbDD   <- length(which(zdd != -99 & !is.na(zdd))) / length(which(!is.na(zdd)))
  znd <- ifelse(Ypos == thr & Ylagpos == thr, Ypos, -99)
  ProbNDND <- length(which(znd != -99 & !is.na(znd))) / length(which(!is.na(znd)))

  c(MeanDAfterZero  = mean(daz, na.rm = TRUE), VarDAfterZero  = var(daz, na.rm = TRUE),
    MeanDBeforeZero = mean(dbz, na.rm = TRUE), VarDBeforeZero = var(dbz, na.rm = TRUE),
    MeanDAfterD     = mean(dad, na.rm = TRUE), VarDAfterD     = var(dad, na.rm = TRUE),
    ProbDD = ProbDD, ProbNDND = ProbNDND)
}

edge_fixtures <- function() {
  set.seed(99)
  list(
    normal      = rnorm(500, 5, 2),
    lognormal   = rlnorm(500, 1, 1),
    uniform     = runif(500),
    constant    = rep(3, 200),
    twovalue    = rep(c(1, 5), 100),
    allneg      = -abs(rnorm(300, 5, 2)),
    mixedsign   = rnorm(400, 0, 3),
    largemag    = rnorm(300, 1e8, 1e7),
    tinymag     = rnorm(300, 1e-8, 1e-9),
    n2          = c(2.5, 7.5),
    n3          = c(1, 2, 9),
    na10        = { v <- rnorm(500, 5, 2); v[sample(500, 50)]  <- NA; v },
    na90        = { v <- rnorm(500, 5, 2); v[sample(500, 450)] <- NA; v },
    naAllButOne = { v <- rep(NA_real_, 50); v[7] <- 4.2; v },
    zerosheavy  = { v <- rnorm(500, 5, 2); v[sample(500, 150)] <- 0; v }
  )
}

test_that("vectorised core matches independent oracles across edge cases", {
  fx  <- edge_fixtures()
  for (nm in names(fx)) {
    f    <- fx[[nm]]
    core <- anyFit:::.basic_stats_core(matrix(f, ncol = 1),
                                       zero_threshold = 0.01, ignore_zeros = FALSE)

    expect_parity(core["Mean", ],  mean(f, na.rm = TRUE))
    expect_parity(core["Var", ],   var(f, na.rm = TRUE))
    expect_parity(core["StDev", ], sd(f, na.rm = TRUE))
    expect_parity(core["Min", ],   suppressWarnings(min(f, na.rm = TRUE)))
    expect_parity(core["Max", ],   suppressWarnings(max(f, na.rm = TRUE)))

    # Derived central moments vs the moments package (population / non-excess)
    expect_parity(core["Skewness", ], moments::skewness(f, na.rm = TRUE))
    expect_parity(core["Kurtosis", ], moments::kurtosis(f, na.rm = TRUE))
    mu <- mean(f, na.rm = TRUE)
    expect_parity(core["Mom3", ],
                  sum((f - mu)^3, na.rm = TRUE) / (sum(!is.na(f)) - 1))

    # Quantiles vs stats::quantile(type = 7) - skip all-NA columns
    if (sum(!is.na(f)) >= 1) {
      expect_parity(core["Q5", ],  stats::quantile(f, .05, na.rm = TRUE, names = FALSE, type = 7))
      expect_parity(core["Q25", ], stats::quantile(f, .25, na.rm = TRUE, names = FALSE, type = 7))
      expect_parity(core["Q50", ], stats::quantile(f, .50, na.rm = TRUE, names = FALSE, type = 7))
      expect_parity(core["Q75", ], stats::quantile(f, .75, na.rm = TRUE, names = FALSE, type = 7))
      expect_parity(core["Q95", ], stats::quantile(f, .95, na.rm = TRUE, names = FALSE, type = 7))
    }

    expect_parity(core["Pdr", ], mean(f <= 0.01, na.rm = TRUE))
    expect_equal(unname(core["NumofMisData", ]), sum(is.na(f)))
    expect_equal(unname(core["NumofData", ]),    length(f))

    # Transition stats vs the golden per-column transcription
    ref <- ref_transition(f, 0.01)
    for (s in names(ref)) expect_parity(core[s, ], ref[[s]])
  }
})

test_that("ignore_zeros NA-mask path equals subsetting the non-zero values", {
  set.seed(7)
  f <- rnorm(500, 5, 2); f[sample(500, 150)] <- 0
  core <- anyFit:::.basic_stats_core(matrix(f, ncol = 1),
                                     zero_threshold = 0.01, ignore_zeros = TRUE)
  s <- f[f > 0.01]                                   # original ts[ts>thr] path
  expect_parity(core["Mean", ],     mean(s, na.rm = TRUE))
  expect_parity(core["Var", ],      var(s, na.rm = TRUE))
  expect_parity(core["Skewness", ], moments::skewness(s, na.rm = TRUE))
  expect_parity(core["Kurtosis", ], moments::kurtosis(s, na.rm = TRUE))
  expect_parity(core["Q95", ],      stats::quantile(s, .95, na.rm = TRUE, names = FALSE, type = 7))
  # counts and Pdr stay on the FULL series
  expect_parity(core["Pdr", ],         mean(f <= 0.01, na.rm = TRUE))
  expect_equal(unname(core["NumofData", ]), length(f))
})

test_that("basic_stats handles an all-NA column without error", {
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, 49) * 3600
  ts    <- xts::xts(matrix(NA_real_, 50, 1, dimnames = list(NULL, "V1")),
                    order.by = dates)
  expect_no_error(r <- basic_stats(ts))
  expect_equal(nrow(r$stats_table), 34L)
})

# ---------------------------------------------------------------------------
# Wide input behaviour
# ---------------------------------------------------------------------------
make_wide_na_ts <- function() {
  set.seed(21)
  n     <- 365 * 2
  vals  <- matrix(rnorm(n * 3, 5, 2), ncol = 3,
                  dimnames = list(NULL, c("S1", "S2", "S3")))
  vals[sample(n, 40), 1] <- 0      # zeros in S1
  vals[sample(n, 40), 2] <- NA     # NAs in S2
  dates <- as.POSIXct("2015-01-01", tz = "UTC") + seq(0, n - 1) * 3600
  xts::xts(vals, order.by = dates)
}

test_that("basic_stats wide stats_table columns are named after the inputs", {
  ts  <- make_multi_ts()
  res <- basic_stats(ts)
  expect_equal(colnames(res$stats_table), colnames(ts))
})

test_that("basic_stats single-column keeps the 'Value' column name", {
  res <- basic_stats(make_plain_ts())
  expect_equal(colnames(res$stats_table), "Value")
})

test_that("basic_stats wide column k equals the single-column call on column k", {
  ts   <- make_wide_na_ts()
  wide <- basic_stats(ts)$stats_table
  expect_equal(colnames(wide), colnames(ts))
  for (j in seq_len(ncol(ts))) {
    single <- basic_stats(ts[, j])$stats_table
    expect_parity(wide[, j], single[, "Value"])
  }
})

test_that("basic_stats plot flag controls plot construction", {
  ts <- make_plain_ts()
  expect_null(basic_stats(ts, plot = FALSE)$plot)

  p <- basic_stats(ts, plot = TRUE)$plot
  expect_s3_class(p, "patchwork")

  multi <- make_multi_ts()
  pm <- basic_stats(multi, plot = TRUE)$plot
  expect_true(is.list(pm) && !inherits(pm, "patchwork"))
  expect_equal(names(pm), colnames(multi))
  expect_equal(length(pm), ncol(multi))
  expect_s3_class(pm[[1]], "patchwork")
})
