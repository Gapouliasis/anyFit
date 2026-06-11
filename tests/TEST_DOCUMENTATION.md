# anyFit Package — Test Suite Documentation

## Overview

This document describes every test implemented for the **anyFit** package. All tests use synthetic, in-memory data generated with `set.seed()` — no external files are required except where explicitly noted (NetCDF tests create a temporary file at runtime and delete it on exit).

**How to run:**
```r
devtools::test()                              # run all tests
devtools::test(filter = "distributions")     # run a single test file
```

**Package version at time of writing:** 0.0.0.9000

---

## Test File Index

| File | Functions covered | Test blocks |
|---|---|---|
| [test_sxts.R](#test_sxtsr) | `sxts`, `is.sxts`, `coords`, `projection`, `elements`, `attributes.sxts`, `print.sxts`, `str.sxts`, `summary.sxts`, `[.sxts`, `lag.sxts`, `diff.sxts`, `Ops.sxts`, `mask.sxts`, `zonal_stats` | 19 |
| [test_acf_functions.R](#test_acf_functionsr) | `CAS_ACF`, `HK_ACF`, `SRD_ACF`, `fit_ACF`, `fit_ACF_monthly` | 17 |
| [test_basic_stats.R](#test_basic_statsr) | `basic_stats`, `lmom_stats` | 15 |
| [test_check_missing.R](#test_check_missingr) | `check_missing`, `plot_missing` | 15 |
| [test_aggregate_xts.R](#test_aggregate_xtsr) | `aggregate_xts`, `period_stats`, `period_apply_nc` | 16 |
| [test_monthly_functions.R](#test_monthly_functionsr) | `monthly_stats`, `monthly_boxplots`, `monthly_ecdf`, `monthly_violins`, `ridge_plots` | 17 |
| [test_distributions.R](#test_distributionsr) | d/p/q/r families (11), `fitlm_*` (15 single-series + genlogi), `fitlm_multi`, `fitlm_nxts`, `fitlm_monthly`, `fit_diagnostics`, `LRatio_check` | 52 |
| [test_delim2xts.R](#test_delim2xtsr) | `delim2xts` | 12 |
| [test_correl_plots.R](#test_correl_plotsr) | `correl_plots` | 12 |
| [test_normalise_xts.R](#test_normalise_xtsr) | `normalise_xts` | 13 |
| [test_sxts_raster.R](#test_sxts_rasterr) | `rasterFromSxts`, `sxtsFromRaster`, `mask.sxts` (error path), `zonal_stats` (continent + dims) | 11 |
| [test_nc_functions.R](#test_nc_functionsr) | `nc2xts`, `nc2xts_nn`, `nc_ggplot`, `basic_stats_nc`, `fitlm_nc` | 15 |

---

## test_sxts.R

**Pre-existing file.** Tests the `sxts` S3 class constructor, accessors, subsetting, arithmetic operators, and spatial masking.

**Shared fixture:** `make_obj()` — 5 time steps × 4 spatial points, 2×2 grid.

| Test | Assertion |
|---|---|
| `sxts()` creates object with correct class and attributes | class is `sxts`/`xts`, nrow=5, ncol=4, elements=4 |
| `sxts()` errors when coord rows != data columns | `stop("Number of coordinate rows")` |
| `sxts()` errors when coords has only one column | `stop("coords must have at least 2 columns")` |
| `is.sxts()` returns TRUE/FALSE | TRUE for sxts, FALSE for xts and numeric |
| `coords()` returns data.frame with x and y columns | nrow=4, names contain x and y |
| `projection()` returns the projection string | exact string match |
| `elements()` returns number of spatial locations | 4L |
| `attributes.sxts()` returns list with elements, coords, projection | 3 named elements |
| `print.sxts()` runs without error and mentions 'sxts' | output contains "sxts" |
| `str.sxts()` runs without error | no error |
| `summary.sxts()` runs without error and reports element count | output contains "Number of elements" |
| `[.sxts` time subsetting preserves class and coords | class=sxts, nrow=3, coords unchanged |
| `[.sxts` spatial subsetting updates coords and elements | ncol=2, nrow(coords)=2, elements=2 |
| `lag.sxts()` returns sxts with same dimensions | class=sxts, same dims, same projection |
| `diff.sxts()` returns sxts with same dims and NA in first row | class=sxts, first row all NA |
| `Ops.sxts` arithmetic preserves class and scales values | class=sxts, values×2 |
| `Ops.sxts` comparison returns sxts of logicals | class=sxts, double type |
| `mask.sxts()` bounding box retains only points inside limits | ncol=1, coords match expected |
| `zonal_stats()` with country returns xts with country as column name | class=xts, ncol=1, colname="Belgium" |

---

## test_acf_functions.R

**Shared fixture:** `make_ar1_ts()` — 5-year daily AR(1) xts, `set.seed(1)`.

| Test | Assertion |
|---|---|
| `CAS_ACF()` returns data.frame with lag and ACF columns | s3 class data.frame, has lag and ACF |
| `CAS_ACF()` lag 0 value equals 1 | ACF[lag==0] == 1 |
| `CAS_ACF()` respects lag_max and returns lag_max+1 rows | nrow == lag_max+1 |
| `CAS_ACF()` produces monotonically decreasing ACF values | all(diff(ACF) < 0) |
| `HK_ACF()` returns data.frame with lag and ACF columns | s3 class data.frame, has lag and ACF |
| `HK_ACF()` lag 0 value equals 1 | ACF[lag==0] == 1 |
| `HK_ACF()` respects lag_max and returns lag_max+1 rows | nrow == lag_max+1 |
| `HK_ACF()` with H=0.5 gives SRD-like rapid decay | H=0.5 decays faster than H=0.9 at lag 5 |
| `SRD_ACF()` returns data.frame with lag and ACF columns | s3 class data.frame, has lag and ACF |
| `SRD_ACF()` lag 0 value equals 1 | ACF[lag==0] == 1 |
| `SRD_ACF()` respects lag_max and returns lag_max+1 rows | nrow == lag_max+1 |
| `SRD_ACF()` larger kappa gives faster decay | kappa=1 decays faster than kappa=0.1 at lag 5 |
| `fit_ACF()` returns list with ACF_params, ACF_fitted, ACF_plot | 3 named list elements |
| `fit_ACF()` ACF_params contains an entry per fitted type | length == 3, named CAS/HK/SRD |
| `fit_ACF()` fitted parameters are numeric and finite | all numeric, all finite |
| `fit_ACF()` with type=list('SRD') returns only SRD in ACF_params | length==1, named SRD |
| `fit_ACF()` ACF_fitted has lag column and one column per fitted type | has lag, CAS, SRD |
| `fit_ACF()` runs without error when ignore_zeros=TRUE | no error |
| `fit_ACF_monthly()` returns list with ACF_params_monthly and ACF_monthly_plot | 2 named elements |
| `fit_ACF_monthly()` ACF_params_monthly has one row per month present | 1–12 rows |
| `fit_ACF_monthly()` parameter values are all numeric | all(is.numeric(unlist(...))) |

---

## test_basic_stats.R

**Shared fixtures:** `make_plain_ts()` 3-year hourly, `make_zeros_ts()` with 20% zeros, `make_na_ts()` with 10% NAs, `make_multi_ts()` 2-column.

| Test | Assertion |
|---|---|
| `basic_stats()` returns list with 'plot' and 'stats_table' | 2 named elements |
| `basic_stats()` stats_table has 34 rows when ignore_zeros=FALSE | nrow == 34 |
| `basic_stats()` stats_table contains all expected row names | 34 specific row names present |
| `basic_stats()` NumofData matches nrow(ts) | exact match |
| `basic_stats()` Mean is close to mean(ts) | within 0.01 tolerance |
| `basic_stats()` stats_table has 26 rows when ignore_zeros=TRUE | nrow == 26 |
| `basic_stats()` Pdr is positive when series contains zeros | Pdr > 0 |
| `basic_stats()` NumofMisData is positive when series contains NAs | NumofMisData > 0 |
| `basic_stats()` pstart/pend restricts the plotted period without error | no error |
| `lmom_stats()` returns a data.frame with 5 rows | s3 class data.frame, nrow == 5 |
| `lmom_stats()` row names are the five L-moment ratios | L-Mean, L-Scale, L-Skew, L-Kurtosis, L-CV |
| `lmom_stats()` has one column per variable in ts | ncol == ncol(ts) |
| `lmom_stats()` L-Mean is close to the sample mean for Normal data | within 0.5 tolerance |
| `lmom_stats()` runs without error when ignore_zeros=TRUE | no error |
| `lmom_stats()` values are all numeric and finite | all numeric, all finite |

---

## test_check_missing.R

**Shared fixtures:** `make_daily_ts()` 2-year daily with 15% NAs; `make_multi_ts()` 2-column with 10% NAs.

| Test | Assertion |
|---|---|
| `check_missing()` returns list with 'prct_missing' element | named list, prct_missing present |
| `check_missing()` prct_missing values are between 0 and 100 | all(pct >= 0 & pct <= 100) |
| `check_missing()` prct_missing is positive when NAs are present | prct_missing[[1]] > 0 |
| `check_missing()` creates list_months element for 'months' period | list_months in names |
| `check_missing()` list_months$prct_missing is an xts object | inherits xts |
| `check_missing()` creates list_years element for 'years' period | list_years in names |
| `check_missing()` handles multiple periods simultaneously | both list_months and list_years |
| `check_missing()` with group_months=TRUE adds grouped_months matrix | is.matrix, nrow==12 |
| `check_missing()` grouped_months values are between 0 and 100 | all in [0, 100] |
| `check_missing()` grouped_months ncol matches ncol of input | ncol == ncol(ts) |
| `check_missing()` with plot=TRUE includes figure element | figure in names, inherits gg |
| `check_missing()` with plot=FALSE omits figure element | figure not in names |
| `check_missing()` works on two-column xts | no error |
| `plot_missing()` returns a ggplot object | inherits gg |
| `plot_missing()` runs without error when series contains NAs | no error |
| `plot_missing()` runs without error when series has no NAs | no error |

---

## test_aggregate_xts.R

**Shared fixture:** `make_hourly_ts()` 3-year hourly xts.

| Test | Assertion |
|---|---|
| `aggregate_xts()` returns list with Combined_Plot and list_days | both in names |
| `aggregate_xts()` aggregated daily xts has fewer rows than hourly input | nrow(agg) < nrow(ts) |
| `aggregate_xts()` creates one list element per requested period | list_days, list_months, list_years |
| `aggregate_xts()` each period element has 'aggregated' and 'figure' | both in names |
| `aggregate_xts()` with FUN='sum' gives sums close to manual calculation | all values ≈ 24 |
| `aggregate_xts()` errors when period_multiplier length mismatches periods | error "same length as periods" |
| `aggregate_xts()` accepts period_multiplier with same length as periods | no error |
| `aggregate_xts()` runs without error when input contains NAs | no error |
| `period_stats()` returns an xts object | inherits xts |
| `period_stats()` contains core statistic columns | NumofData, Mean, Var, Q50 present |
| `period_stats()` monthly aggregation has correct row count | nrow == 36 (3 years × 12 months) |
| `period_stats()` yearly aggregation has correct row count | nrow == 3 |
| `period_stats()` PercOfMissingData is positive when NAs are present | sum > 0 |
| `period_stats()` statistic values are all numeric | all(sapply(..., is.numeric)) |
| `period_apply_nc()` on sxts returns sxts with fewer rows than input | inherits sxts, nrow < nrow(input) |

---

## test_monthly_functions.R

**Shared fixtures:** `make_daily_ts()` 5-year daily with 10% zeros; `make_multi_ts()` 2-column.

| Test | Assertion |
|---|---|
| `monthly_stats()` returns list with 'base_stats' and 'fbase' | 2 named elements |
| `monthly_stats()` base_stats has 12 columns (one per month) | ncol == 12 |
| `monthly_stats()` base_stats column names are month names | colnames == month.name |
| `monthly_stats()` fbase is a ggplot/patchwork object | inherits gg or patchwork |
| `monthly_stats()` runs without error when ignore_zeros=TRUE | no error |
| `monthly_stats()` with aggregated=TRUE returns agg_stats, faggre, lag1 | 3 named elements |
| `monthly_stats()` agg_stats has 12 columns for aggregated=TRUE | ncol == 12 |
| `monthly_stats()` lag1 is a data.frame with 12 rows | s3 data.frame, nrow == 12 |
| `monthly_boxplots()` returns a ggplot object for single-column ts | inherits gg |
| `monthly_boxplots()` returns a ggplot object for two-column ts | inherits gg |
| `monthly_boxplots()` runs without error when ignore_zeros=TRUE | no error |
| `monthly_ecdf()` returns a ggplot or patchwork object | inherits gg or patchwork |
| `monthly_ecdf()` runs without error when ignore_zeros=TRUE | no error |
| `monthly_violins()` returns a ggplot object for single-column ts | inherits gg |
| `monthly_violins()` returns a ggplot object for two-column ts | inherits gg |
| `monthly_violins()` runs without error when ignore_zeros=TRUE | no error |
| `ridge_plots()` returns list with 'plot_all' and 'plot_monthly' | 2 named elements |
| `ridge_plots()` plot_all is a ggplot object | inherits gg |
| `ridge_plots()` plot_monthly is a list with one entry per column | length == ncol(ts) |
| `ridge_plots()` each plot_monthly element is a ggplot | all inherit gg |
| `ridge_plots()` runs without error when ignore_zeros=TRUE | no error |

---

## test_distributions.R

**Shared fixtures:** `make_norm_ts()` Normal(5,2) n=500; `make_gamma_ts()` Gamma(2,3) n=500; `make_zeros_ts()` with 20% zeros; `make_multi_ts()` 2-column n=500.

### Section A — d/p/q/r consistency

For each distribution family (exp, gamma3, rayleigh, GEV, gumbel, burr, lognorm):

| Test pattern | Assertion |
|---|---|
| `d*()` density | non-negative, length matches input |
| `p*()` CDF | monotone increasing, bounded [0, 1] |
| `q*(p*(x)) ≈ x` round-trip | within 1e-5 tolerance |
| `r*()` samples | length == n, all ≥ location |

### Section B — fitlm_* single-series

For each of: `fitlm_norm`, `fitlm_exp`, `fitlm_gamma`, `fitlm_gamma3`, `fitlm_gev`, `fitlm_lognorm`, `fitlm_weibull`, `fitlm_rayleigh`, `fitlm_burr`, `fitlm_dagum`, `fitlm_gengamma`, `fitlm_expweibull`, `fitlm_GPD`:

| Test | Assertion |
|---|---|
| Returns standard structure | has Distribution, Param, TheorLMom, DataLMom, GoF |
| GoF has all six metrics | MLE, CM, KS, MSEquant, DiffOfMax, MeanDiffOf10Max |
| Parameters numeric and finite | all(is.numeric), all(is.finite) |
| ignore_zeros=TRUE runs without error | no error on zeros_ts |

For `fitlm_genlogi`: returns Param, TheorLMom, DataLMom (no GoF — by design).

### Section C — multi-distribution helpers

| Test | Assertion |
|---|---|
| `fitlm_multi()` with diagnostics returns 5 elements | parameter_list, GoF_summary, diagnostics, QQplot, PPplot |
| `fitlm_multi()` GoF_summary has one column per candidate | ncol == length(candidates) |
| `fitlm_multi()` without diagnostics returns only params and GoF | no QQplot element |
| `fitlm_nxts()` returns params, diagnostic_plots, QQ_plots, PP_plots | 4 elements |
| `fitlm_nxts()` params has one entry per column | length == ncol(ts) |
| `fitlm_monthly()` returns params_monthly, GoF_monthly, monthly_QQplot, monthly_PPplot | 4 elements |
| `fitlm_monthly()` params_monthly has one entry per candidate | length == length(candidates) |
| `fitlm_monthly()` params_monthly column names are month names | colnames in month.name |
| `fit_diagnostics()` returns Diagnostic_Plots, GoF, QQplot, PPplot | 4 named elements |
| `fit_diagnostics()` GoF has CramerVonMises and KolmogorovSmirnov | 2 named elements |
| `fit_diagnostics()` plots are ggplot objects | inherits gg |
| `LRatio_check()` returns distributions and multi_plots | 2 named elements |
| `LRatio_check()` distributions contains Dagum, GGamma, ExpWeibull, BurrXII | 4 names |

---

## test_delim2xts.R

**Fixtures:** `write_tab_file()`, `write_gap_file()`, `write_novalue_file()`, `write_leap_file()`, `write_date_second_file()` — all using `tempfile()`, cleaned up with `on.exit()`.

| Test | Assertion |
|---|---|
| Returns xts with correct dimensions | inherits xts, nrow==5, ncol==2 |
| Column names match file header | colnames == c("V1", "V2") |
| Data values match file content | as.numeric matches expected |
| `strict_step=FALSE` accepts non-uniform time steps | no error |
| `strict_step=TRUE` errors on non-uniform steps | error "Time step is not strict" |
| `col_names=FALSE` assigns default column names | ncol >= 1 |
| `date_index=2` uses second column as time index | inherits xts, nrow==5, data correct |
| `exc_leaps=TRUE` excludes Feb 29 | "02-29" not in output dates |
| `exc_leaps=FALSE` retains Feb 29 | "02-29" in output dates |
| `save_Xts=TRUE` writes a file at the given path | file.exists == TRUE |

---

## test_correl_plots.R

**Shared fixture:** Two 3-year daily xts with Pearson ≈ 0.8, `set.seed(7)`.

| Test | Assertion |
|---|---|
| Returns list with all five elements | combined, scatter_plot, copula_plot, normal_plot, correl_table |
| `correl_table` has 6 rows | nrow == 6 |
| `correl_table` row names are correct | Pearson, Spearman, Semi_1–4 |
| Pearson correlation is close to expected ~0.8 | 0.7 < rho < 0.95 |
| Pearson equals manual cor() on common data | exact match after rounding |
| `check_common=TRUE` uses only common dates | no error on date-offset series |
| `scatter_plot` is a ggplot | inherits gg |
| `copula_plot` is a ggplot | inherits gg |
| `normal_plot` is a ggplot | inherits gg |
| `combined` is a patchwork | inherits patchwork |
| `ignore_zeros=TRUE` runs without error | no error |
| `correl_table` still 6 rows with ignore_zeros=TRUE | nrow == 6 |

---

## test_normalise_xts.R

**Shared fixtures:** `make_daily_ts()` 5-year daily single-column; `make_multi_ts()` 2-column; `make_zeros_ts()` with 15% zeros.

| Test | Assertion |
|---|---|
| Returns a list | is.list == TRUE |
| List length equals ncol(ts) | length(res) == ncol(ts) |
| List names match colnames(ts) | names(res) == colnames(ts) |
| Each element is an xts | inherits xts |
| Monthly mode returns same row count as input | nrow == nrow(ts) |
| Monthly mode returns numeric values | all(is.numeric) |
| Monthly mode values are not all NA | not all NA |
| `dist_period=NA` runs without error | no error |
| Global mode returns same row count as input | nrow == nrow(ts) |
| Global mode values are not all NA | not all NA |
| Independent normalizations per column (means near 0) | abs(mean) < 0.5 per column |
| `ignore_zeros=TRUE` runs without error | no error |
| `ignore_zeros=TRUE` returns fewer rows than original | nrow < nrow(ts) |

---

## test_sxts_raster.R

**Shared fixture:** `make_obj()` — 5 time steps × 4 spatial points on a 2×2 regular grid.

| Test | Assertion |
|---|---|
| `rasterFromSxts()` returns a Raster object | inherits Raster |
| `rasterFromSxts()` RasterBrick has one layer per time step | nlayers == nrow(obj) |
| `rasterFromSxts()` preserves projection string | projection contains "longlat" |
| `sxtsFromRaster()` returns an sxts object | inherits sxts |
| `sxtsFromRaster()` round-trip preserves number of spatial points | ncol == ncol(obj) |
| `sxtsFromRaster()` round-trip preserves number of time steps | nrow == nrow(obj) |
| `mask.sxts()` errors when all mask arguments are NA | error "Provide xlim" |
| `zonal_stats()` with continent='Europe' returns xts with continent column | inherits xts, ncol==1, colname=="Europe" |
| `zonal_stats()` country result has correct time dimensions | nrow == n_time |

---

## test_nc_functions.R

**Notes on fixtures:** NetCDF-dependent tests use `make_tmp_nc()` which creates a 2×2×5 temporary NetCDF file via `ncdf4::nc_create()` and cleans up with `on.exit(file.remove(tmp))`. Tests that only use sxts/raster input (`nc_ggplot`, `basic_stats_nc`, `fitlm_nc`) require no external files. All tests are skipped if required packages are not installed.

| Test | Assertion |
|---|---|
| `nc2xts()` returns list with 'ncdf_sxts' element | list, ncdf_sxts in names |
| `nc2xts()` ncdf_sxts inherits sxts class | inherits sxts |
| `nc2xts()` ncdf_sxts has correct time steps | nrow == 5 |
| `nc2xts()` errors on non-existent variable name | error "not found" |
| `nc2xts()` with country filter returns ≤ all spatial columns | ncol(filtered) <= ncol(full) |
| `nc2xts_nn()` returns an xts object | inherits xts |
| `nc2xts_nn()` returns one column per supplied coordinate | ncol == nrow(coords) |
| `nc_ggplot()` returns ggplot/patchwork when given an sxts | inherits gg or patchwork |
| `nc_ggplot()` with common_legend=TRUE runs without error | no error |
| `basic_stats_nc()` with sxts input returns a Raster object | inherits Raster |
| `basic_stats_nc()` raster has at least 26 layers | nlayers >= 26 |
| `fitlm_nc()` with sxts input returns fit_results and gof_plots | 2 named elements |
| `fitlm_nc()` fit_results has one entry per candidate | length == length(candidates) |
| `fitlm_nc()` each fit_results element has raster_params | "raster_params" in names, inherits Raster |

---

## Design Principles

- **Synthetic data only:** Every test fixture uses `runif`, `rnorm`, or arithmetic construction seeded with `set.seed()`. No files bundled with the test suite except those created at test time by `tempfile()`.
- **One assertion per `test_that` block:** Each block tests one logical behavior, making failures easy to diagnose.
- **Style consistency:** All new files match the formatting of the pre-existing `test_sxts.R`: `library(testthat)` only at the top, shared fixture functions named `make_*()`, comment-delimited section headers.
- **Error path coverage:** Every documented `stop()` call has a corresponding `expect_error()` test.
- **Ignore-zeros coverage:** Functions that accept `ignore_zeros` are tested with a zero-containing series.
