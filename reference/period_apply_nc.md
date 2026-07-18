# Temporal aggregation of gridded data

Aggregates gridded time series data to coarser temporal scales. Accepts
an `sxts` object, a `Raster*` object, or a NetCDF file path; when a
filename is supplied, data are imported through
[`nc2xts`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md)
first. The aggregation period is controlled by `period` and
`period_multiplier`. For common aggregation functions (`"mean"`,
`"sum"`, `"min"`, `"max"`, `"median"`, `"var"`, `"sd"`), the function
very efficient column-wise function from `matrixStats` package, which
ensures good computational efficiency. Unsupported functions (e.g.
custom functions) fall back to a per-column lapply which is
significantly slower.

## Usage

``` r
period_apply_nc(
  data = NULL,
  filename = NA,
  varname = NA,
  period = "months",
  period_multiplier = 1,
  FUN = "mean",
  ...
)
```

## Arguments

- data:

  An `sxts` object, or a `Raster*` object. Leave `NULL` when supplying
  `filename` and `varname`.

- filename:

  Optional NetCDF file path to import if `data` is not provided.

- varname:

  Optional variable name to extract from `filename`.

- period:

  Period string passed to
  [`xts::endpoints`](https://rdrr.io/pkg/xts/man/endpoints.html)
  (e.g.`"months"`, `"years"`).

- period_multiplier:

  Integer multiplier for custom period lengths (default 1).

- FUN:

  Summary function as a string (`"mean"`, `"sum"`, etc.) or a function.

- ...:

  Additional arguments passed to
  [`nc2xts`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md)
  when `filename` and `varname` are supplied.

## Value

An `sxts` object when the input is `sxts` or NetCDF; an `xts` when the
input is a `Raster*`.

## Examples

``` r
# Synthetic sxts
set.seed(123)
dates <- seq(as.POSIXct("2000-01-01"), as.POSIXct("2000-12-31"), by = "day")
coords <- data.frame(x = c(1, 2, 3, 4), y = c(1, 1, 1, 1))
vals <- matrix(rnorm(length(dates) * 4, 10, 5), nrow = length(dates), ncol = 4)
ncdf_sxts <- sxts(vals, order.by = dates, coords = coords,
                  projection = "+proj=longlat")

# Monthly mean
result <- period_apply_nc(ncdf_sxts, period = "months", FUN = "mean")

# Quarterly sum
result <- period_apply_nc(ncdf_sxts, period = "months",
                          period_multiplier = 3, FUN = "sum")
```
