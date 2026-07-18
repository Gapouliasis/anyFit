# Apply per-statistic rounding rules

Rounds each row of the core statistics matrix according to a fixed
per-statistic precision table. `PercOfMissingData` is left unrounded.

## Usage

``` r
.basic_stats_round(core, ignore_zeros = FALSE)
```

## Arguments

- core:

  A numeric matrix with rownames identifying each statistic.

- ignore_zeros:

  Logical. When `TRUE`, quantile rows are rounded to 5 decimal places;
  otherwise transition statistic rows are rounded to 5 decimal places.

## Value

The rounded numeric matrix (same dimensions as `core`).
