# normalise_xts

Empirical normal-score transform of an xts time series: values are
mapped to standard normal quantiles via their empirical cumulative
distribution function (ranking with average ties). Optionally applied
per calendar month to preserve seasonal structure in cyclo-stationary
processes.

## Usage

``` r
normalise_xts(
  ts,
  dist_period = "monthly",
  ignore_zeros = FALSE,
  zero_threshold = 0.01
)
```

## Arguments

- ts:

  An xts object containing the time series data. Multi-column xts are
  supported; each column is normalised independently.

- dist_period:

  Character; the distribution period for normalisation. Use `"monthly"`
  for per-calendar-month normalisation or `NA` for the entire series.
  Default `"monthly"`.

- ignore_zeros:

  Logical; if `TRUE`, zero values (at or below `zero_threshold`) are
  excluded before normalisation. Default `FALSE`.

- zero_threshold:

  Numeric; threshold below which values are treated as zero. Default
  `0.01`.

## Value

A named list of normalised xts objects, one per column of the input.

## Examples

``` r
# Synthetic daily precipitation
set.seed(42)
n <- 365 * 3
dates <- seq(as.POSIXct("2000-01-01", tz = "UTC"), by = "day", length.out = n)
precip <- pmax(0, rnorm(n, mean = 3, sd = 5))
ts <- xts::xts(precip, order.by = dates)

# Monthly normal-score transform
ns <- normalise_xts(ts, dist_period = "monthly")
hist(zoo::coredata(ns[[1]]), main = "Normalised values")


# Full-series transform
ns_full <- normalise_xts(ts, dist_period = NA)
```
