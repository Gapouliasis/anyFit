# Changelog

## anyFit 2.0.0

- Introduced the **`sxts`** (spatial xts) S3 class, which stores
  coordinates and projection alongside the time series and preserves
  them through subsetting, arithmetic, lagging, differencing and
  aggregation.
- Added native spatial operations for `sxts`: bounding-box and shapefile
  masking
  ([`mask.sxts()`](https://gapouliasis.github.io/anyFit/reference/mask.sxts.md)),
  zonal statistics over polygons
  ([`zonal_stats()`](https://gapouliasis.github.io/anyFit/reference/zonal_stats.md)),
  and bidirectional raster conversion
  ([`rasterFromSxts()`](https://gapouliasis.github.io/anyFit/reference/rasterFromSxts.md),
  [`sxtsFromRaster()`](https://gapouliasis.github.io/anyFit/reference/sxts.md)).
- [`nc2xts()`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md)
  reads NetCDF directly and clips on load by country, continent,
  bounding box or shapefile, returning an `sxts`.
- Gridded counterparts of the point-based tools:
  [`basic_stats_nc()`](https://gapouliasis.github.io/anyFit/reference/basic_stats_nc.md),
  [`monthly_stats_nc()`](https://gapouliasis.github.io/anyFit/reference/monthly_stats_nc.md),
  [`period_apply_nc()`](https://gapouliasis.github.io/anyFit/reference/period_apply_nc.md),
  [`fitlm_nc()`](https://gapouliasis.github.io/anyFit/reference/fitlm_nc.md)
  and
  [`fitlm_monthly_nc()`](https://gapouliasis.github.io/anyFit/reference/fitlm_monthly_nc.md),
  with optional parallel fitting.
- Distribution fitting by the method of L-moments across sixteen
  families, including the flexible Burr type-XII, Dagum, Generalised
  Gamma and Exponentiated Weibull, with L-ratio diagnostic diagrams
  ([`LRatio_check()`](https://gapouliasis.github.io/anyFit/reference/LRatio_check.md)).
- A rebuilt documentation website with articles driven entirely by the
  bundled E-OBS example dataset.
