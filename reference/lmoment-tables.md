# Precomputed L-moment lookup tables (internal)

Reference tables consumed internally by the distribution-fitting
routines. The `*_LSpace` objects delimit the feasible L-moment ratio
support region of a distribution and back
[`LRatio_check()`](https://gapouliasis.github.io/anyFit/reference/LRatio_check.md).
The `*_InitValues` tables map sampled L-moment ratios to good starting
parameters for the numerical (L-BFGS-B) fit of the distributions that
lack closed-form L-moments (Burr XII, Dagum, Exponentiated Weibull,
Generalised Gamma), as used by
[`fitlm_nc()`](https://gapouliasis.github.io/anyFit/reference/fitlm_nc.md)
and its relatives.

## Format

- `LMom_LSpace`: a data frame (17994 x 3) with columns `L-CV`,
  `L-Skewness`, `Dist`.

- `BurrXII_LSpace`, `Dagum_LSpace`:
  [sp::SpatialPolygonsDataFrame](https://edzer.github.io/sp/reference/SpatialPolygons.html)
  objects delimiting each distribution's L-ratio support polygon.

- `Burr_InitValues`, `Dagum_InitValues`, `ExpWeibull_InitValues`,
  `GG_InitValues`: data frames with 8 columns (`lambda_1`, `lambda_2`,
  `tau_3`, `tau_4`, `lcv`, `scale`, `shape1`, `shape2`).
