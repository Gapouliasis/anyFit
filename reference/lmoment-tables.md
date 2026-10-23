# Precomputed L-moment lookup tables (internal)

Reference tables consumed internally by the distribution-fitting
routines. The `*_LSpace` objects delimit the attainable (L-CV,
L-skewness) region of a distribution over the shape bounds of its
fitter. They back
[`LRatio_check()`](https://gapouliasis.github.io/anyFit/reference/LRatio_check.md)
and decide the `"outside_domain"` status of
[`fitlm_burr()`](https://gapouliasis.github.io/anyFit/reference/fitlm_burr.md),
[`fitlm_dagum()`](https://gapouliasis.github.io/anyFit/reference/fitlm_dagum.md),
[`fitlm_gengamma()`](https://gapouliasis.github.io/anyFit/reference/fitlm_gengamma.md)
and
[`fitlm_expweibull()`](https://gapouliasis.github.io/anyFit/reference/fitlm_expweibull.md).
They are built by `data-raw/make_LSpace.R`, which maps a 120 x 120 log
grid over the shape box through the exact L-moment engine and unions the
grid-cell images. The `*_InitValues` tables hold unit-scale shape pairs
with their L-moments. They seed the numerical (L-BFGS-B) fit of these
four families by nearest-neighbour matching on the scale-free (L-CV,
\\\tau_3\\, \\\tau_4\\), and cover the fitter boxes
(`data-raw/make_*_InitValues.R`).

## Format

- `LMom_LSpace`: a data frame (17994 x 3) with columns `L-CV`,
  `L-Skewness`, `Dist`.

- `BurrXII_LSpace`, `Dagum_LSpace`:
  [sp::SpatialPolygonsDataFrame](https://edzer.github.io/sp/reference/SpatialPolygons.html)
  objects delimiting each distribution's attainable L-ratio region (the
  generalised gamma and exponentiated Weibull polygons are internal).

- `Burr_InitValues`, `Dagum_InitValues`, `ExpWeibull_InitValues`,
  `GG_InitValues`: data frames with 8 columns (`lambda_1`, `lambda_2`,
  `tau_3`, `tau_4`, `lcv`, `scale`, `shape1`, `shape2`), at `scale = 1`.
