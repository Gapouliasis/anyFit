# monthly_stats_nc

Computes calendar-month statistics on gridded data — either an sxts
object, a NetCDF file, or a raster object. For each calendar month
present in the dataset, the function calls
[`basic_stats_nc`](https://gapouliasis.github.io/anyFit/reference/basic_stats_nc.md),
returning a multi-layer raster of per-cell statistics: mean, min, max,
standard deviation, L-moments, probability dry, empirical quantiles, and
(unless `ignore_zeros = TRUE`) wet/dry transition statistics. The outer
loop over up to twelve calendar months is separable from the
already-vectorised per-cell computations inside `basic_stats_nc`, so
parallelisation is offered at the month level only. NetCDF input is
loaded once via
[`nc2xts`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md) and
then sliced by month index; raster input is converted to sxts with
automatic date parsing.

## Usage

``` r
monthly_stats_nc(
  data = NULL,
  filename = NA,
  varname = NA,
  ignore_zeros = FALSE,
  zero_threshold = 0.01,
  parallel = FALSE,
  ncores = 2,
  ...
)
```

## Arguments

- data:

  An sxts object, or a Raster\* object. Leave `NULL` when supplying
  `filename` and `varname`.

- filename:

  Optional; path to a NetCDF file to import when `data` is not provided.

- varname:

  Optional; name of the variable to extract from `filename`.

- ignore_zeros:

  Logical; if `TRUE`, zeros are ignored in the per-cell statistics.
  Default `FALSE`.

- zero_threshold:

  Numeric; threshold below which values are treated as zero. Default
  `0.01`.

- parallel:

  Logical; whether to compute per-month statistics in parallel. Default
  `FALSE`.

- ncores:

  Integer; number of cores for parallel computation. Default `2`.

- ...:

  Additional arguments passed to
  [`nc2xts`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md)
  when `filename` and `varname` are supplied.

## Value

A named list with one element per calendar month present in the data,
named after [`month.name`](https://rdrr.io/r/base/Constants.html). Each
element is the
[`basic_stats_nc`](https://gapouliasis.github.io/anyFit/reference/basic_stats_nc.md)
output raster for that month.

## Examples

``` r
# Synthetic sxts with 4 cells and 2 years of daily data
set.seed(123)
n <- 365 * 2
dates <- seq(as.POSIXct("2000-01-01", tz = "UTC"), by = "day", length.out = n)
vals <- matrix(pmax(0, rnorm(n * 4, mean = 3, sd = 5)), nrow = n, ncol = 4)
coords <- data.frame(x = c(0, 1, 0, 1), y = c(50, 50, 51, 51))
ts_sxts <- sxts(data = vals, order.by = dates, coords = coords,
                projection = "+proj=longlat +datum=WGS84")

ms_nc <- monthly_stats_nc(data = ts_sxts)
names(ms_nc)
#>  [1] "January"   "February"  "March"     "April"     "May"       "June"     
#>  [7] "July"      "August"    "September" "October"   "November"  "December" 
```
