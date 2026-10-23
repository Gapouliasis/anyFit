# fitlm_dagum

Fits the three-parameter Dagum distribution to an xts series by the
method of L-moments. The L-ratios (L-CV \\\lambda_2/\lambda_1\\ and
\\\tau_r\\) do not depend on the scale, so the two shapes are solved
from the sample L-ratios selected by `order`, and the scale follows
analytically as the sample \\\lambda_1\\ divided by the unit-scale
\\\lambda_1\\ of the fitted shapes. The shapes are found by L-BFGS-B in
log space, minimising the sum of squared differences between sample and
theoretical L-ratios; the theoretical L-moments are closed-form
Beta-function PWM expressions. The start is the row of the unit-scale
`Dagum_InitValues` table nearest in min-max scaled (L-CV, \\\tau_3\\,
\\\tau_4\\); a fit that stalls is restarted once from the second-nearest
row. The shapes are bounded to `1e-4 <= shape1 <= 500` and
`1e-6 <= shape2 <= 0.95` (the mean exists only for \\\theta \< 1\\); the
scale is unbounded. Large `shape1` is an asymptote of the family. Zero
values below `zero_threshold` may be excluded via `ignore_zeros`.

## Usage

``` r
fitlm_dagum(x, ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))
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
`Diag` (solver diagnostics: `status`, `resid`, `convergence`, `evals`,
`restarts`; see
[`fitlm_burr`](https://gapouliasis.github.io/anyFit/reference/fitlm_burr.md)).
The attainable region that decides `"outside_domain"` is `Dagum_LSpace`.

## Examples

``` r
x <- xts::xts(rdagum(365, scale = 1.5, shape1 = 2, shape2 = 0.3),
              order.by = as.Date("2020-01-01") + 0:364)
if (FALSE) { # \dontrun{
fit <- fitlm_dagum(x)
fit$Param
} # }
```
