# L-Moment Statistics

Computes the L-mean, L-scale, L-skewness, L-kurtosis, and L-CV for each
column of an xts object using
[`lmom::samlmu`](https://rdrr.io/pkg/lmom/man/samlmu.html).

## Usage

``` r
lmom_stats(ts, ignore_zeros = FALSE, zero_threshold = 0.01)
```

## Arguments

- ts:

  An xts object containing the time series data.

- ignore_zeros:

  Logical. If `TRUE`, zeros are excluded before computing L-moments.
  Default `FALSE`.

- zero_threshold:

  Numeric. Values below this threshold are treated as zero. Default
  0.01.

## Value

A data frame with rows for L-mean, L-scale, L-skewness, L-kurtosis, and
L-CV, and one column per input series.

## Examples

``` r
# Synthetic daily data
set.seed(42)
ts <- xts::xts(rgamma(3650, shape = 2, scale = 5),
          order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = 3650))
ts_lmoms <- lmom_stats(ts, ignore_zeros = TRUE)
```
