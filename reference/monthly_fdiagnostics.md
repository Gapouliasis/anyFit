# Monthly Distribution Diagnostics

Computes Q-Q and P-P diagnostic plots and Cramer-von Mises /
Kolmogorov-Smirnov goodness-of-fit statistics for a fitted distribution,
applied separately to each calendar month and assembled into 12-panel
grids.

## Usage

``` r
monthly_fdiagnostics(
  ts,
  distr,
  params,
  ignore_zeros = FALSE,
  zero_threshold = 0.01,
  nrow = 3,
  ncol = 4
)
```

## Arguments

- ts:

  An xts object containing the time series data.

- distr:

  Character string naming the distribution (e.g. `'gamma3'`).

- params:

  A named list of monthly parameter vectors, as returned by
  [`fitlm_monthly`](https://gapouliasis.github.io/anyFit/reference/fitlm_monthly.md)\$params_monthly.

- ignore_zeros:

  Logical. If `TRUE`, zeros are excluded. Default `FALSE`.

- zero_threshold:

  Numeric. Values below this threshold are treated as zero. Default
  0.01.

- nrow:

  Number of rows in the plot grid. Default 3.

- ncol:

  Number of columns in the plot grid. Default 4.

## Value

A list with elements `monthly_QQplot` (12-panel Q-Q grid),
`monthly_PPplot` (12-panel P-P grid), and `GoF_monthly` (list of
per-month CM/KS statistics).

## Examples

``` r
# Synthetic daily data
set.seed(42)
ts <- xts::xts(rgamma(3650, shape = 2, scale = 5),
          order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = 3650))
monthly_fits <- fitlm_monthly(ts, candidates = 'gamma3', ignore_zeros = TRUE)
monthly_fcheck <- monthly_fdiagnostics(ts, distr = 'gamma3',
                                       params = monthly_fits$params_monthly,
                                       ignore_zeros = TRUE)
```
