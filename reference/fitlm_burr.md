# fitlm_burr

Fits the three-parameter Burr Type XII distribution to an xts series by
the method of L-moments. The L-ratios (L-CV \\\lambda_2/\lambda_1\\ and
\\\tau_r\\) do not depend on the scale, so the two shapes are solved
from the sample L-ratios selected by `order`, and the scale follows
analytically as the sample \\\lambda_1\\ divided by the unit-scale
\\\lambda_1\\ of the fitted shapes. The shapes are found by L-BFGS-B in
log space, minimising the sum of squared differences between sample and
theoretical L-ratios; the theoretical L-moments are closed-form
Beta-function PWM expressions. The start is the row of the unit-scale
`Burr_InitValues` table nearest in min-max scaled (L-CV, \\\tau_3\\,
\\\tau_4\\); a fit that stalls is restarted once from the second-nearest
row. The shapes are bounded to `0.5 <= shape1 <= 1e4` and
`0.001 <= shape2 <= 1` (the mean exists only for \\\theta \< 1\\); the
scale is unbounded. The `shape2 -> 0` edge is the Weibull limit of the
family. Zero values below `zero_threshold` may be excluded via
`ignore_zeros`.

## Usage

``` r
fitlm_burr(x, ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))
```

## Arguments

- x:

  An xts object containing the time series data.

- ignore_zeros:

  Logical. If `TRUE`, values below `zero_threshold` are excluded.
  Default `FALSE`.

- zero_threshold:

  Numeric. Threshold below which values are treated as zero. Default
  `0.01`.

- order:

  Integer vector of L-moment orders to match. Order 2 matches the L-CV
  and order \\r \ge 3\\ matches \\\tau_r\\; \\\lambda_1\\ always sets
  the scale. At least two ratios are required, otherwise an error is
  raised: two give an exact match (e.g. `1:3`), more a least-squares
  fit. Default `1:5`.

## Value

A list with elements `Distribution`, `Param` (named list of fitted
`scale`, `shape1`, `shape2`), `TheorLMom` (theoretical L-moments),
`DataLMom` (sample L-moments), `GoF` (goodness-of-fit metrics) and
`Diag` (solver diagnostics): `status` is `"matched"` (fitted and sample
L-ratios agree to within 1e-3), `"matched_at_bound"` (the same, with a
shape on its bound), `"outside_domain"` (the sample (L-CV, \\\tau_3\\)
lies outside the attainable region of the family, `BurrXII_LSpace`, see
[`LRatio_check`](https://gapouliasis.github.io/anyFit/reference/LRatio_check.md);
the fit is the closest the family can get) or `"not_converged"`; `resid`
is the largest absolute difference between fitted and sample L-ratios
over the matched orders; `convergence` is the `optim` code; `evals` the
number of objective evaluations; `restarts` 0 or

1.  With more than two ratios the status follows the `optim` code
    instead of `resid`, since least squares cannot match every ratio.

## Examples

``` r
x <- xts::xts(rburr(365, scale = 1.5, shape1 = 2, shape2 = 0.3),
              order.by = as.Date("2020-01-01") + 0:364)
if (FALSE) { # \dontrun{
fit <- fitlm_burr(x)
fit$Param
} # }
```
