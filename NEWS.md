# anyFit 2.0.0

* Introduced the **`sxts`** (spatial xts) S3 class, which stores coordinates and
  projection alongside the time series and preserves them through subsetting,
  arithmetic, lagging, differencing and aggregation.
* Added native spatial operations for `sxts`: bounding-box and shapefile masking
  (`mask.sxts()`), zonal statistics over polygons (`zonal_stats()`), and
  bidirectional raster conversion (`rasterFromSxts()`, `sxtsFromRaster()`).
* `nc2xts()` reads NetCDF directly and clips on load by country, continent,
  bounding box or shapefile, returning an `sxts`.
* Gridded counterparts of the point-based tools: `basic_stats_nc()`,
  `monthly_stats_nc()`, `period_apply_nc()`, `fitlm_nc()` and
  `fitlm_monthly_nc()`, with optional parallel fitting.
* Distribution fitting by the method of L-moments across sixteen families,
  including the flexible Burr type-XII, Dagum, Generalised Gamma and
  Exponentiated Weibull, with L-ratio diagnostic diagrams (`LRatio_check()`).
* A rebuilt documentation website with articles driven entirely by the bundled
  E-OBS example dataset.
