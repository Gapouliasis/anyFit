# Changelog

## anyFit 2.1.0

- [`fitlm_burr()`](https://gapouliasis.github.io/anyFit/reference/fitlm_burr.md),
  [`fitlm_dagum()`](https://gapouliasis.github.io/anyFit/reference/fitlm_dagum.md),
  [`fitlm_gengamma()`](https://gapouliasis.github.io/anyFit/reference/fitlm_gengamma.md),
  [`fitlm_gengamma_loc()`](https://gapouliasis.github.io/anyFit/reference/fitlm_gengamma_loc.md)
  and
  [`fitlm_expweibull()`](https://gapouliasis.github.io/anyFit/reference/fitlm_expweibull.md)
  use a new shared L-moment solver. The two shapes are solved from the
  scale-free L-ratios by L-BFGS-B on log-shapes (Stacy (k, d)
  coordinates for the generalised gamma), and the scale is derived
  analytically. The fit is restarted once from the second-nearest seed
  if it stalls. The fits reproduce the sample L-ratios to about 1e-7
  instead of 1e-4, take 10–15 objective evaluations instead of 70–110,
  and no longer end on an artificial scale bound. This changes the
  generalised gamma and exponentiated Weibull fits substantially; Burr
  XII and Dagum fits change by about 1e-4.
- Wider shape bounds: Burr XII `shape1 <= 1e4`, Dagum `shape2 <= 0.95`,
  exponentiated Weibull `shape2 <= 1e4`. The generalised gamma requires
  `shape1 / shape2 >= 0.01`, where its quadrature is accurate. The
  `*_InitValues` seed tables are extended to cover the new bounds.
- The fitters return `Diag` (`status`, `resid`, `convergence`, `evals`,
  `restarts`). `status` is one of `"matched"`, `"matched_at_bound"`,
  `"outside_domain"` (the sample L-ratios lie outside the family’s
  attainable region; the fit is the closest the family can get) or
  `"not_converged"`.
- `order` now selects L-ratios: order 2 matches the L-CV, order r \>= 3
  matches tau_r, and lambda_1 always sets the scale. At least two ratios
  are required.
- The attainable-region polygons `BurrXII_LSpace`, `Dagum_LSpace`,
  `GGamma_LSpace` and `ExpWeibull_LSpace` used by
  [`LRatio_check()`](https://gapouliasis.github.io/anyFit/reference/LRatio_check.md)
  are regenerated over the new bounds (`data-raw/make_LSpace.R`). The
  exponentiated Weibull region was previously far too small.
- [`qburr()`](https://gapouliasis.github.io/anyFit/reference/BurXII.md),
  [`rburr()`](https://gapouliasis.github.io/anyFit/reference/BurXII.md)
  and
  [`dburr()`](https://gapouliasis.github.io/anyFit/reference/BurXII.md)
  are evaluated in log space and no longer overflow for large
  `shape1 * shape2`.

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
