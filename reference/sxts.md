# Spatial xts (sxts) class

Spatial extension of the xts class that stores coordinates, projection,
and element count as object attributes rather than slots, keeping memory
overhead low and allowing standard xts operations to work transparently.
The constructor validates that the number of coordinate rows matches the
number of data columns and that `coords` has at least two columns (x and
y).

`str.sxts` Prints the structure of an sxts object, displaying spatial
extent, projection, number of elements, and the range of dates.

`summaruy.sxts` Summarises an sxts object by reporting element count,
date range, time step regularity, spatial projection and extent, and
value range.

`print.sxts` Prints a compact summary of an sxts object including
spatial extent, projection, and dimensions.

`attributes.sxts` Returns a named list of the three spatial attributes
(`elements`, `coords`, `projection`) stored in an sxts object.

`coords.sxts` Returns the spatial coordinates data frame attached to an
sxts object.

`projection.sxts` Returns the coordinate reference system (CRS) string
attached to an sxts object.

`elements.sxts` Returns the number of spatial locations (columns) in an
sxts object.

`[` Subset method for sxts objects that preserves spatial attributes.
When `j` is specified, coordinates are subset to match the selected
columns and the spatial element count is updated accordingly.

`is.sxts` Tests whether an object inherits from the `sxts` class.

`lag.sxts` Lag method that dispatches to `NextMethod("lag")` and
re-attaches spatial attributes via `restore_sxts`.

`diff.sxts` Diff method that dispatches to `NextMethod("diff")` and
re-attaches spatial attributes via `restore_sxts`.

`Ops.sxts` Group method for arithmetic and comparison operators.
Dispatches via `NextMethod("Ops")` and preserves spatial attributes from
whichever operand is an sxts object.

`rasterFromSxts.sxts` Converts an sxts object to a `RasterBrick`,
preserving coordinates, projection, and time-indexed layers.

`sxtsFromRaster` Converts a `Raster*` object to an sxts object,
auto-detecting date formats from layer names via
[`lubridate::parse_date_time`](https://lubridate.tidyverse.org/reference/parse_date_time.html).

## Usage

``` r
sxts(data, order.by, coords = NULL, projection = NULL)

# S3 method for class 'sxts'
str(object, ...)

# S3 method for class 'sxts'
summary(object, ...)

# S3 method for class 'sxts'
print(x, ...)

attributes.sxts(obj)

# S3 method for class 'sxts'
coords(x)

# S3 method for class 'sxts'
projection(x)

# S3 method for class 'sxts'
elements(x)

# S3 method for class 'sxts'
x[i, j, drop = FALSE, ...]

is.sxts(x)

# S3 method for class 'sxts'
lag(x, k = 1, ...)

# S3 method for class 'sxts'
diff(x, lag = 1, differences = 1, ...)

# S3 method for class 'sxts'
Ops(e1, e2)

# S3 method for class 'sxts'
rasterFromSxts(x, ...)

sxtsFromRaster(raster, ...)
```

## Arguments

- data:

  Matrix or data.frame of time series data.

- order.by:

  POSIXct vector for the time index.

- coords:

  A data.frame or matrix with coordinates (x, y columns).

- projection:

  Character string specifying the coordinate reference system.

- object, x, obj:

  An xts or sxts object.

- ...:

  Additional arguments passed to `xts()`.

- i, j:

  Row and column indices for subsetting.

- drop:

  Logical; whether to drop dimensions (passed to `NextMethod("[")`).

- k:

  Number of periods to shift (passed to `lag`).

- lag:

  Integer; number of periods to lag (passed to `diff`).

- differences:

  Integer; order of differencing (passed to `diff`).

- e1, e2:

  Operands for arithmetic and comparison operations.

- raster:

  A `Raster*` object to convert to an sxts object.

- xlim:

  Numeric vector of length 2; x-axis bounding box limits.

- ylim:

  Numeric vector of length 2; y-axis bounding box limits.

- shapefile_name:

  Character; path to a shapefile for spatial masking.

- mask:

  An `sf` or `Spatial*` object used as a spatial mask.

## Value

An sxts object.

## Examples

``` r
# Synthetic sxts
set.seed(123)
dates <- seq(as.POSIXct("2000-01-01"), as.POSIXct("2000-01-10"), by = "day")
coords <- data.frame(x = c(1, 2, 3, 4), y = c(1, 1, 1, 1))
vals <- matrix(rnorm(length(dates) * 4, 10, 5), nrow = length(dates), ncol = 4)
ncdf_sxts <- sxts(vals, order.by = dates, coords = coords,
                  projection = "+proj=longlat")
str(ncdf_sxts)
#> 'sxts' Spatial xts object
#>  Description:  This is a spatial xts (sxts) object 
#>  Projection:  +proj=longlat 
#>  Elements:  4  spatial locations
#> Spatial extent is: 1 , 4 , 1 , 1  (xmin, xmax, ymin, ymax) 
#>  Time series data:
#> Range of dates is from: 2000-01-01 to: 2000-01-10 
summary(ncdf_sxts)
#> Summary of sxts object:
#> Number of elements: 4 
#> Number of dates: 10 
#> Range of dates is from: 2000-01-01 to: 2000-01-10 
#> Time step is: Strict 
#> Time step is: 24 hours 
#> Number of dates: 10 
#> Spatial projection: +proj=longlat 
#> Spatial extent is: 1 , 4 , 1 , 1  (xmin, xmax, ymin, ymax) 
#> Min value is: 0.1669142 
#> Max value is: 18.93457 

# Subsetting preserves spatial attributes
ncdf_sxts[1:5, 1:2]
#> Spatial xts (sxts) object
#> Description: This is a spatial xts (sxts) object 
#> Projection: +proj=longlat 
#> Number of locations: 2 
#> Dimensions: 5 x 2 
#> 
#> Spatial extent:
#>   X range: 1 2 
#>   Y range: 1 1 
#> 
#> Time series data:
```
