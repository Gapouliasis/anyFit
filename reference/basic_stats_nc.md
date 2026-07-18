# Compute basic statistics for gridded (NetCDF or raster) data

Computes a comprehensive set of over 30 summary statistics for every
grid cell in a spatial time series. The function accepts an `sxts`
object, a `Raster*` object, or a NetCDF file path with variable name.
When `filename` and `varname` are supplied, the data are first imported
via [`nc2xts`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md)
and additional arguments are forwarded to that function. For raster
inputs the layer names are parsed for dates, and the raster values are
converted to an xts matrix. The returned raster has one layer per
statistic; layer names correspond to the rows of the `basic_stats`
output.

## Usage

``` r
basic_stats_nc(
  data = NULL,
  filename = NA,
  varname = NA,
  ignore_zeros = FALSE,
  zero_threshold = 0.01,
  ...
)
```

## Arguments

- data:

  An `sxts` object, or a `Raster*` object. Leave as `NULL` when
  supplying `filename` and `varname` instead.

- filename:

  A NetCDF file name to import. Ignored when `data` is supplied
  directly.

- varname:

  The name of the variable to extract from `filename`.

- ignore_zeros:

  Logical. If `TRUE`, values at or below `zero_threshold` are excluded
  from distributional statistics. Default is `FALSE`.

- zero_threshold:

  Numeric threshold below which values are treated as zero. Default is
  `0.01`.

- ...:

  Additional arguments passed to
  [`nc2xts`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md)
  when `filename` and `varname` are provided.

## Value

A `RasterBrick` with one layer per statistic. The layer names correspond
to the statistic labels (e.g. `"Mean"`, `"StDev"`, `"Skewness"`,
`"Pdr"`, etc.). The spatial reference (projection and extent) is
inherited from the input.

## Examples

``` r
if (FALSE) { # \dontrun{
# From an sxts object (e.g. output of nc2xts)
s <- nc2xts(filename = "precip.nc", varname = "pr")
r <- basic_stats_nc(data = s, ignore_zeros = TRUE)
raster::plot(r[["Mean"]])

# Direct from NetCDF (shortcut — nc2xts is called internally)
r <- basic_stats_nc(filename = "precip.nc", varname = "pr",
                     ignore_zeros = TRUE)
} # }
```
