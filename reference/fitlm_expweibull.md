# fitlm_expweibull

Fits the three-parameter Exponentiated Weibull distribution to an xts
series by numerical minimisation of L-moments. L-moments are computed at
each optimisation step via a fast tanh-sinh (double-exponential)
quadrature engine (`lmom_expweibull`), chosen for its well-conditioned
behaviour across the full positive parameter range. The L-BFGS-B
optimiser minimises the normalised root-sum-square error between the
sample and theoretical L-moments of orders specified by `order`. A
two-step seeding procedure initialises the shape parameters from the
`ExpWeibull_InitValues` lookup table and derives the scale analytically.
The distribution is meant to be fitted only to non-negative data. Zero
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

  Integer vector of L-moment orders matched by the optimiser. Default
  `1:3` (exact method-of-L-moments for the three parameters).

## Value

A list with elements `Distribution`, `Param` (named list of fitted
`scale`, `shape1`, `shape2`), `TheorLMom` (theoretical L-moments),
`DataLMom` (sample L-moments), and `GoF` (goodness-of-fit metrics).

## Examples

``` r
x <- xts::xts(rexpweibull(365, scale = 3, shape1 = 1.5, shape2 = 2),
              order.by = as.Date("2020-01-01") + 0:364)
if (FALSE) { # \dontrun{
fit <- fitlm_expweibull(x)
fit$Param
} # }
```
