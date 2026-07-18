# fitlm_dagum

Fits the three-parameter Dagum distribution to an xts series by
numerical minimisation of L-moments. Closed-form L-moment expressions
via Beta-function PWMs are evaluated at each optimisation step. The
L-BFGS-B optimiser minimises the normalised root-sum-square error
between the sample and theoretical L-moments of orders specified by
`order`. A two-step seeding procedure initialises the shape parameters
by min-max nearest-neighbour search over the scale-free L-ratios in the
`Dagum_InitValues` lookup table, and the scale is derived analytically
from the sample first L-moment and the matched shapes' unit-scale
L-moment. The mean exists only when \\\theta \< 1\\. Zero values below
`zero_threshold` may be excluded via `ignore_zeros`.

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

  Integer vector of L-moment orders matched by the optimiser. Default
  `1:5`.

## Value

A list with elements `Distribution`, `Param` (named list of fitted
`scale`, `shape1`, `shape2`), `TheorLMom` (theoretical L-moments),
`DataLMom` (sample L-moments), and `GoF` (goodness-of-fit metrics).

## Examples

``` r
x <- xts::xts(rdagum(365, scale = 1.5, shape1 = 2, shape2 = 0.3),
              order.by = as.Date("2020-01-01") + 0:364)
if (FALSE) { # \dontrun{
fit <- fitlm_dagum(x)
fit$Param
} # }
```
