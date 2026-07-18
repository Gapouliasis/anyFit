# Aggregate sxts over shapefile polygons (zonal statistics)

For each polygon, finds all sxts spatial points that fall within it and
aggregates them row-wise (per time step) using `FUN`. Returns a plain
`xts` with polygon names as column names. Points that do not fall within
any polygon are silently dropped.

The polygon source can be supplied in three ways (mutually exclusive):

- `shapefile` — an `sf`/`Spatial*` object or a file path; `name_col`
  must also be given.

- `country` — a country name matching `world_data$name`; yields one
  output column named after the country.

- `continent` — a continent name matching `world_data$continent`; yields
  one output column per country within the continent.

Aggregates sxts cells within polygon zones (shapefile, country, or
continent), returning a plain xts with one column per zone. Points
outside all zones are dropped.

## Usage

``` r
zonal_stats(x, ...)

# S3 method for class 'sxts'
zonal_stats(
  x,
  shapefile = NULL,
  name_col = NULL,
  country = NA,
  continent = NA,
  FUN = mean,
  ...
)
```

## Arguments

- x:

  An `sxts` object.

- ...:

  Additional arguments forwarded to `FUN` (e.g. `na.rm = TRUE`).

- shapefile:

  An `sf` object, `Spatial*` object, or file path to a shapefile.

- name_col:

  Character; column in `shapefile` whose values become output column
  names.

- country:

  Character; country name (matches `world_data$name`).

- continent:

  Character; continent name (matches `world_data$continent`).

- FUN:

  Aggregation function applied row-wise within each polygon (default:
  `mean`).

## Value

A plain `xts` with one column per polygon that contains at least one
sxts point.
