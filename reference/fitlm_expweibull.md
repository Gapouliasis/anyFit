# fitlm_expweibull

Fits the three-parameter Exponentiated Weibull distribution to an xts
series by the method of L-moments. The L-ratios (L-CV
\\\lambda_2/\lambda_1\\ and \\\tau_r\\) do not depend on the scale, so
the two shapes are solved from the sample L-ratios selected by `order`,
and the scale follows analytically as the sample \\\lambda_1\\ divided
by the unit-scale \\\lambda_1\\ of the fitted shapes. The theoretical
L-moments come from a tanh-sinh (double-exponential) quadrature engine,
well conditioned across the full positive parameter range. The shapes
are found by L-BFGS-B in log space, minimising the sum of squared
differences between sample and theoretical L-ratios, within
`0.05 <= shape1 <= 50` and `0.02 <= shape2 <= 1e4`; the scale is
unbounded. As `shape2` grows the family approaches the Gumbel skewness
(\\\tau_3 \to 0.17\\), so large `shape2` values are identified only
slowly. The start is the row of the unit-scale `ExpWeibull_InitValues`
table nearest in min-max scaled (L-CV, \\\tau_3\\, \\\tau_4\\); a fit
that stalls is restarted once from the second-nearest row. The
distribution is meant to be fitted only to non-negative data. Zero
values below `zero_threshold` may be excluded via `ignore_zeros`.

## Usage

``` r
fitlm_expweibull(
  x,
  ignore_zeros = FALSE,
  zero_threshold = 0.01,
  order = c(1:3)
)
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

  Integer vector of L-moment orders to match, mapped to L-ratios as in
  [`fitlm_burr`](https://gapouliasis.github.io/anyFit/reference/fitlm_burr.md).
  Default `1:3` (exact method-of-L-moments for the three parameters).

## Value

A list with elements `Distribution`, `Param` (named list of fitted
`scale`, `shape1`, `shape2`), `TheorLMom` (theoretical L-moments),
`DataLMom` (sample L-moments), `GoF` (goodness-of-fit metrics) and
`Diag` (solver diagnostics: `status`, `resid`, `convergence`, `evals`,
`restarts`; see
[`fitlm_burr`](https://gapouliasis.github.io/anyFit/reference/fitlm_burr.md)).
The attainable region that decides `"outside_domain"` is
`ExpWeibull_LSpace`.

## Examples

``` r
x <- xts::xts(rexpweibull(365, scale = 3, shape1 = 1.5, shape2 = 2),
              order.by = as.Date("2020-01-01") + 0:364)
if (FALSE) { # \dontrun{
fit <- fitlm_expweibull(x)
fit$Param
} # }
```
