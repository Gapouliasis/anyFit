# fitlm_gengamma

Fits the three-parameter Generalised Gamma (Stacy) distribution to an
xts series by numerical minimisation of L-moments. As no closed-form
L-moment expressions exist for this distribution, the L-moments are
computed at each optimisation step via Gauss-Legendre quadrature of the
PWMs. The L-BFGS-B optimiser minimises the normalised root-sum-square
error between the sample and theoretical L-moments of orders specified
by `order`. A two-step seeding procedure is employed: the shape
parameters are initialised by a min-max nearest-neighbour search over
the scale-free L-ratios in the `GG_InitValues` lookup table, and the
scale is then derived analytically from the sample first L-moment and
the matched shapes' unit-scale L-moment. Zero values below
`zero_threshold` may be excluded via `ignore_zeros`.

## Usage

``` r
fitlm_gengamma(x, ignore_zeros = FALSE, zero_threshold = 0.01, order = c(1:5))
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
x <- xts::xts(rgengamma(365, scale = 2, shape1 = 0.8, shape2 = 0.5),
              order.by = as.Date("2020-01-01") + 0:364)
if (FALSE) { # \dontrun{
fit <- fitlm_gengamma(x)
fit$Param
} # }
```
