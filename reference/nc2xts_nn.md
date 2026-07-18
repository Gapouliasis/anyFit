# nc2xts_nn

Extracts nearest-neighbour point time series from one or more NetCDF
files at the supplied coordinates, returning a plain `xts` object. For
each requested point, the value of the grid cell whose centre is closest
in Euclidean distance is read directly via ncdf4 hyperslab indexing — no
interpolation or extrapolation is performed. Points that fall outside
the grid extent (with a half-cell buffer) yield an all-`NA` column in
the output. Data are read one point at a time, each as a single-cell
indexed hyperslab of all time steps, keeping the I/O footprint small.
Multiple files (e.g. data downloaded in time chunks) are supported.

## Usage

``` r
nc2xts_nn(filename, varname, coords)
```

## Arguments

- filename:

  Character vector; path(s) to one or more NetCDF files. If multiple
  files are given, results are row-bound and ordered by date.

- varname:

  Character; name of the variable to extract (a single variable name
  used across all files).

- coords:

  A matrix or data.frame of coordinates; the first two columns are taken
  as x (longitude) and y (latitude). Each row is a separate extraction
  point.

## Value

An xts time series object with one column per supplied coordinate point.
Columns are named after the variable for single-point extraction or
`"Point1"`, `"Point2"`, etc. for multiple points.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single point from a NetCDF file
pt <- nc2xts_nn("era5_daily.nc", varname = "tp",
                coords = cbind(23.7, 37.9))

# Multiple points across multiple files
coords <- cbind(c(23.7, 24.0), c(37.9, 38.0))
pts <- nc2xts_nn(c("era5_2000.nc", "era5_2001.nc"),
                 varname = "tp", coords = coords)
} # }
```
