# Spatial masking of sxts objects

Spatially subsets an sxts object by bounding box (`xlim`, `ylim`),
country name, continent name, or a shapefile. When using polygons,
point-in-polygon membership is determined via
[`sf::st_within`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
with CRS-aware coordinate transformation between the sxts CRS and the
mask CRS. Points falling outside the mask are dropped.

## Usage

``` r
mask.sxts(obj, xlim = NA, ylim = NA, shapefile_name = NA, mask = NA)
```
