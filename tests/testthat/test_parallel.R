library(testthat)

# ===========================================================================
# Failure-mode tests for the parallel fitlm_nxts backend (round 4).
#
# Pure-logic tests (F7, F8, F11, F13) run everywhere.
# Multiprocess tests (F1-F6, F9, F10, F12, F15) need real future workers that
# can library(anyFit); they skip_on_cran and skip when that is not possible
# (e.g. devtools::test() against an uninstalled / stale package build).
# ===========================================================================

make_grid <- function(K = 12, n = 200, seed = 1) {
  set.seed(seed)
  m <- vapply(seq_len(K), function(j) abs(rnorm(n, mean = 5 + j, sd = 2)) + 0.5,
              numeric(n))
  colnames(m) <- paste0("c", seq_len(K))
  xts::xts(m, order.by = as.POSIXct("2000-01-01", tz = "UTC") + seq_len(n) * 86400)
}

# Compact per-column numeric signature (also carries names -> order check).
sig_params <- function(out) {
  vapply(out$params, function(v) {
    g <- as.matrix(v$GoF_summary)
    p <- unlist(lapply(v$params, function(z) unlist(z$Param)))
    sum(as.numeric(g), na.rm = TRUE) + sum(p, na.rm = TRUE)
  }, numeric(1))
}

# Probe: can future multisession workers load anyFit and see fitlm_multi?
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
fit_serial <- function(g) fitlm_nxts(g, candidates = cands, ignore_zeros = TRUE,
                                     parallel = FALSE, diagnostic_plots = FALSE)
fit_par <- function(g, shared) fitlm_nxts(g, candidates = cands, ignore_zeros = TRUE,
                                          parallel = TRUE, shared_memory = shared,
                                          diagnostic_plots = FALSE)

# ---------------------------------------------------------------------------
# Pure-logic tests (always run)
# ---------------------------------------------------------------------------

test_that("F7: big.matrix preserves grid values incl. NA", {
  skip_if_not_installed("bigmemory")
  g <- make_grid(6, 100)
  m <- zoo::coredata(g); m[3, 2] <- NA; m[10, 5] <- NA
  bm <- bigmemory::as.big.matrix(m, type = "double")
  expect_identical(bm[, ], m)
  expect_identical(bm[, 5, drop = FALSE], m[, 5, drop = FALSE])
})

test_that("F8: splitIndices covers every column exactly, incl. nw > ncol", {
  for (ncol in c(1, 2, 5, 12, 37)) {
    for (nw in c(1, 2, 3, 8, 64)) {
      idx <- parallel::splitIndices(ncol, min(nw, ncol))
      expect_identical(sort(unlist(idx)), seq_len(ncol))
      expect_false(any(duplicated(unlist(idx))))
    }
  }
})

test_that("F11: GOF_tests sort-dedup is byte-identical incl. short series (n<10)", {
  ref_gof <- function(x, fit) {
    u <- rank(zoo::coredata(x), na.last = NA, ties.method = "average") / (length(x) + 1)
    qq <- do.call(qgamma, c(list(p = u), fit))
    ts2 <- sort(qq, decreasing = TRUE)[1:10]; ss <- sort(x, decreasing = TRUE)[1:10]
    list(MSEquant = sum((qq - x)^2) / length(x),
         DiffOfMax = ((max(qq) - max(x)) / max(x)) * 100,
         MeanDiffOf10Max = sum(abs(ts2 - ss)) / 10)
  }
  fit <- list(scale = 3, shape = 2)
  for (n in c(5, 8, 10, 50)) {
    set.seed(n); x <- round(stats::rgamma(n, shape = 2, scale = 3), 1)
    new <- GOF_tests(x, fit, "gamma"); old <- ref_gof(x, fit)
    for (k in names(old))
      expect_equal(new[[k]], old[[k]], tolerance = 1e-12, info = paste("n=", n, k))
  }
})

test_that("F13: shared_memory does not warn-fallback under a normal local plan", {
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  g <- make_grid(6, 120)
  expect_no_warning(fit_par(g, shared = TRUE))
})

# ---------------------------------------------------------------------------
# Multiprocess tests
# ---------------------------------------------------------------------------

test_that("F2/F1: serial == parallel (both layers), order preserved", {
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  g <- make_grid(12, 200)
  ser <- fit_serial(g)
  expect_equal(fit_par(g, TRUE)$params,  ser$params)   # mmap
  expect_equal(fit_par(g, FALSE)$params, ser$params)   # chunks
  expect_identical(names(fit_par(g, TRUE)$params), colnames(g))  # F1 order
})

test_that("F3: result is invariant to worker count", {
  skip_unless_parallel()
  g <- make_grid(12, 200)
  old <- future::plan(future::multisession, workers = 2); r2 <- sig_params(fit_par(g, TRUE))
  future::plan(future::multisession, workers = 3);        r3 <- sig_params(fit_par(g, TRUE))
  future::plan(old)
  expect_equal(r2, r3)
})

test_that("F4: ncores/workers > ncol does not crash and matches serial", {
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 4)
  on.exit(future::plan(old), add = TRUE)
  g <- make_grid(2, 150)
  expect_equal(fit_par(g, TRUE)$params,  fit_serial(g)$params)
  expect_equal(fit_par(g, FALSE)$params, fit_serial(g)$params)
})

test_that("F5: single-column grid works in parallel", {
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  g <- make_grid(1, 200)
  expect_equal(fit_par(g, TRUE)$params, fit_serial(g)$params)
})

test_that("F6: all-zero and all-NA columns handled at correct positions", {
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  g <- make_grid(6, 150)
  g[, 2] <- 0          # all-zero -> empty_fit under ignore_zeros
  g[, 5] <- NA_real_   # all-NA   -> empty_fit
  ser <- fit_serial(g); par <- fit_par(g, TRUE)
  expect_equal(par$params, ser$params)
  expect_true(is.na(ser$params[[2]]$GoF_summary["CM", "norm"]))
  expect_true(is.na(ser$params[[5]]$GoF_summary["CM", "norm"]))
})

test_that("F9: user-set plan is respected and left intact", {
  skip_unless_parallel()
  before <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(before), add = TRUE)
  cur <- future::plan()
  invisible(fit_par(make_grid(6, 120), TRUE))
  expect_identical(class(future::plan()), class(cur))   # not overridden
})

test_that("F10: shared_memory leaves no backing files in tempdir", {
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  before <- list.files(tempdir(), pattern = "^grid_")
  invisible(fit_par(make_grid(6, 120), TRUE))
  after <- list.files(tempdir(), pattern = "^grid_")
  expect_identical(sort(after), sort(before))
})

test_that("F12: chunks mode survives a low future.globals.maxSize", {
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  oopt <- options(future.globals.maxSize = 1024); on.exit(options(oopt), add = TRUE)
  g <- make_grid(6, 400)   # chunk >> 1024 bytes
  expect_error(fit_par(g, FALSE), NA)   # the in-function Inf override prevents error
})

test_that("F16: worker closures do not capture the grid (no giant globals)", {
  skip_unless_parallel()
  # The factory closure must be tiny regardless of any grid.
  w <- anyFit:::.make_col_worker("norm", TRUE, 0.01, FALSE)
  expect_lt(length(serialize(w, NULL)), 1e6)   # constant ~0.1 MB; a captured grid is >> 1 MB

  # End-to-end: a grid LARGER than a deliberately tiny maxSize must still run in
  # shared_memory mode (only the descriptor crosses). If the worker captured the
  # grid this errors. shared_memory mode does NOT raise maxSize itself.
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  oopt <- options(future.globals.maxSize = 2 * 1024^2)  # 2 MB
  on.exit(options(oopt), add = TRUE)
  g <- make_grid(40, 8000)                              # ~2.6 MB of doubles > 2 MB
  expect_error(fitlm_nxts(g, candidates = "norm", ignore_zeros = TRUE,
                          parallel = TRUE, shared_memory = TRUE,
                          diagnostic_plots = FALSE), NA)
})

test_that("F15: parallel is deterministic and emits no RNG warning", {
  skip_unless_parallel()
  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  g <- make_grid(8, 200)
  expect_no_warning(r1 <- fit_par(g, TRUE))
  r2 <- fit_par(g, TRUE)
  expect_equal(sig_params(r1), sig_params(r2))
})
