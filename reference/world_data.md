# World country and continent boundary polygons

A global map of country and territory boundaries whose attributes back
anyFit's spatial helpers. The `name` and `continent` columns are the
ones matched when a `country` or `continent` is supplied to
[`mask.sxts()`](https://gapouliasis.github.io/anyFit/reference/mask.sxts.md)
and
[`zonal_stats()`](https://gapouliasis.github.io/anyFit/reference/zonal_stats.md),
so the exact spelling in these columns is what those arguments expect
(e.g. the United Kingdom is stored as
`"U.K. of Great Britain and Northern Ireland"`).

## Usage

``` r
world_data
```

## Format

A
[sp::SpatialPolygonsDataFrame](https://edzer.github.io/sp/reference/SpatialPolygons.html)
with 256 features and 8 attribute columns:

- iso3:

  Three-letter ISO 3166-1 alpha-3 country code.

- status:

  Sovereignty or status of the territory.

- color_code:

  Colour group used for cartographic display.

- name:

  Country or territory name (matched by `country =` arguments).

- continent:

  Continent name (matched by `continent =` arguments).

- region:

  Sub-regional grouping.

- iso_3166_1\_:

  ISO 3166-1 numeric code.

- french_shor:

  Short country name in French.

## Source

Bundled global country/territory boundary polygons.

## See also

[`mask.sxts()`](https://gapouliasis.github.io/anyFit/reference/mask.sxts.md),
[`zonal_stats()`](https://gapouliasis.github.io/anyFit/reference/zonal_stats.md)

## Examples

``` r
# Continents available for `continent =` masking
sort(unique(world_data$continent))
#> [1] "Africa"     "Americas"   "Antarctica" "Asia"       "Europe"    
#> [6] "Oceania"   

# Exact country name used by `country =` arguments
grep("Great Britain", world_data$name, value = TRUE)
#> [1] "U.K. of Great Britain and Northern Ireland"
```
