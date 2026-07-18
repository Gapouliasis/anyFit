# nc2xts

Reads a variable from a NetCDF file directly via ncdf4 and converts it
into an sxts (spatial xts) object. The CRS is detected automatically
from the NetCDF metadata — checking the variable's `grid_mapping`
attribute and global `crs_wkt` / `spatial_ref` attributes — and
cross-checked against any user-supplied `projection` argument, raising
an error on mismatch. Axis coordinates are regularised to remove
floating-point rounding errors that would otherwise cause downstream
raster functions to reject the grid as irregular. Spatial subsetting is
handled either by a polygon mask (country, continent, or shapefile) or
by explicit `xlim`/`ylim` bounding-box coordinates. For polygon masks
the bounding box is computed in the grid's own CRS and used as a
hyperslab prefilter, avoiding a full-grid read; the prefilter is skipped
with a warning for rotated-pole grids where the projected bounding box
cannot be mapped to a contiguous index window.

## Usage

``` r
nc2xts(
  filename,
  varname,
  shapefile = NA,
  country = NA,
  continent = NA,
  xlim = NA,
  ylim = NA,
  projection = "+proj=longlat +datum=WGS84"
)
```

## Arguments

- filename:

  Character; path to the NetCDF file.

- varname:

  Character; name of the variable to extract.

- shapefile:

  Optional; path to a shapefile for spatial masking.

- country:

  Optional; country name for masking (must match `world_data$name`).

- continent:

  Optional; continent name for masking (must match
  `world_data$continent`).

- xlim:

  Optional; numeric vector of length 2 giving the x-axis (longitude)
  range in the grid's coordinate units. Requires `ylim`.

- ylim:

  Optional; numeric vector of length 2 giving the y-axis (latitude)
  range in the grid's coordinate units. Requires `xlim`.

- projection:

  Character; CRS string used as a fallback when the NetCDF file contains
  no projection metadata. If the file does contain CRS metadata, both
  are compared and an error is raised on mismatch. Default
  `"+proj=longlat +datum=WGS84"`.

## Value

An sxts object with rows = time steps and columns = spatial locations,
with spatial attributes (coordinates, projection, elements).

## Details

For large files, supplying a `country`, `continent`, `shapefile`, or
`xlim`/`ylim` restricts the NetCDF read to the bounding box of the
requested region, so only that hyperslab is loaded from disk rather than
the whole grid. For polygon masks the bounding box is computed in the
grid's own CRS, so the prefilter is correct for both geographic and
projected grids. The prefilter is skipped (full read, with a warning)
for rotated-pole grids.

## Examples

``` r
if (FALSE) { # \dontrun{
# Read precipitation from a NetCDF file
sxts_data <- nc2xts("era5_daily.nc", varname = "tp")

# Mask by country
sxts_fr <- nc2xts("era5_daily.nc", varname = "tp", country = "France")

# Mask by bounding box
sxts_bb <- nc2xts("era5_daily.nc", varname = "tp",
                  xlim = c(-10, 10), ylim = c(40, 60))
} # }
```
