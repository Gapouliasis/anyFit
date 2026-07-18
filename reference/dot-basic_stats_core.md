# Column-wise statistics engine

Computes an unrounded statistics matrix with one row per statistic and
one column per input series. When `ignore_zeros = TRUE`, distributional
statistics are computed on the non-zero values only (sub-threshold
values masked to `NA`) and transition statistics are omitted.

## Usage

``` r
.basic_stats_core(m, zero_threshold = 0.01, ignore_zeros = FALSE)
```

## Arguments

- m:

  A numeric matrix (columns = series, rows = observations).

- zero_threshold:

  Numeric threshold. Default is `0.01`.

- ignore_zeros:

  Logical. Default is `FALSE`.

## Value

A numeric matrix with rownames identifying each statistic and colnames
inherited from `m`.
