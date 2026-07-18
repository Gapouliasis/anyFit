# Period-based summary statistics

Computes a comprehensive suite of period-based summary statistics for
one or more time series. Statistics include count, missing-data
diagnostics, min, max, mean, variance, standard deviation, coefficient
of variation, third moment, skewness, kurtosis, L-moments (mean, scale,
L3, L4, L-CV, L-skewness, L-kurtosis), quantiles (5, 25, 50, 75, 95),
and inter-quartile range. Each statistic is computed column-wise per
period using `matrixStats` with `na.rm = TRUE`; L-moments are computed
independently via
[`lmom::samlmu`](https://rdrr.io/pkg/lmom/man/samlmu.html) since no
column-wise matrixStats equivalent exists. The result is a named list of
xts objects, one per statistic, each with one row per period and one
column per input series. Accepts any period supported by
[`xts::endpoints`](https://rdrr.io/pkg/xts/man/endpoints.html) (e.g.
`"months"`, `"years"`) with an optional multiplier for custom period
lengths.

## Usage

``` r
period_stats(ts, period = "months", period_multiplier = 1)
```

## Arguments

- ts:

  An xts object containing the time series data.

- period:

  A period string passed to `endpoints` (default `"months"`).

- period_multiplier:

  Integer multiplier for custom period lengths (default 1).

## Value

A named list with one element per statistic (`NumofData`,
`NumofMisData`, `PercOfMissingData`, `Min`, `Max`, `Mean`, `Var`,
`StDev`, `Variation`, `Mom3`, `Skewness`, `Kurtosis`, `LMean`, `LScale`,
`L3`, `L4`, `LVariation`, `LSkewness`, `LKurtosis`, `Q5`, `Q25`, `Q50`,
`Q75`, `Q95`, `IQR`). Each element is an xts with one row per period and
one column per input series.

## Examples

``` r
# Synthetic xts
set.seed(123)
dates <- seq(as.POSIXct("2000-01-01"), as.POSIXct("2000-12-31"), by = "day")
vals <- matrix(rnorm(length(dates) * 3, 10, 5), nrow = length(dates), ncol = 3)
colnames(vals) <- c("A", "B", "C")
ts <- xts::xts(vals, order.by = dates)

# Statistics per 3-month periods
pstats <- period_stats(ts, period = "months", period_multiplier = 3)
pstats$Mean
#>                    A         B         C
#> 2000-03-31 10.311162  9.725297 10.487451
#> 2000-06-30  9.752405 10.568295 10.615225
#> 2000-09-30 10.009313  9.285933  9.605579
#> 2000-12-31 10.548654  9.782270 10.839648
pstats$Q95
#>                   A        B        C
#> 2000-03-31 17.68795 17.33576 18.54074
#> 2000-06-30 19.38241 17.57733 19.48497
#> 2000-09-30 18.88848 17.96998 18.10591
#> 2000-12-31 18.38538 16.38306 18.99246
```
