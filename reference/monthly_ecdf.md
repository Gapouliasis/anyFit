# Monthly Empirical CDF Grid

Generates a 12-panel grid of empirical cumulative distribution function
plots, one per calendar month. Supports multi-column xts input with
series distinguished by colour.

## Usage

``` r
monthly_ecdf(ts, ignore_zeros = FALSE, zero_threshold = 0.01)
```

## Arguments

- ts:

  An xts object containing the time series data.

- ignore_zeros:

  Logical. If `TRUE`, zeros are excluded. Default `FALSE`.

- zero_threshold:

  Numeric. Values below this threshold are treated as zero. Default
  0.01.

## Value

A patchwork grid of ggplot ECDF panels (4 rows by 3 columns).

## Examples

``` r
# Synthetic daily data
set.seed(42)
ts <- xts::xts(rgamma(3650, shape = 2, scale = 5),
          order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = 3650))
monthly_ecdf(ts)

```
