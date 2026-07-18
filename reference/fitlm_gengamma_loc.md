# fitlm_gengamma_loc

Fits the four-parameter Generalised Gamma distribution with location to
an xts series. The location parameter is supplied by the user and is
subtracted from the data before fitting; the remaining three parameters
(scale, shape1, shape2) are then fitted by delegating to
[`fitlm_gengamma`](https://gapouliasis.github.io/anyFit/reference/fitlm_gengamma.md).
The theoretical first L-moment is adjusted by adding back the location
after the fit. Zero values below `zero_threshold` may be excluded via
`ignore_zeros`. Goodness-of-fit is assessed via `GOF_tests` evaluated on
the full four-parameter set.

## Usage

``` r
fitlm_gengamma_loc(
  x,
  location,
  ignore_zeros = FALSE,
  zero_threshold = 0.01,
  order = c(1:5)
)
```

## Arguments

- x:

  An xts object containing the time series data.

- location:

  Numeric. The location parameter of the distribution.

- ignore_zeros:

  Logical. If `TRUE`, values below `zero_threshold` are excluded.
  Default `FALSE`.

- zero_threshold:

  Numeric. Threshold below which values are treated as zero. Default
  `0.01`.

- order:

  Integer vector of L-moment orders passed to `fitlm_gengamma`. Default
  `1:5`.

## Value

A list with elements `Distribution`, `Param` (named list of fitted
`location`, `scale`, `shape1`, `shape2`), `TheorLMom` (theoretical
L-moments), `DataLMom` (sample L-moments), and `GoF` (goodness-of-fit
metrics).

## Examples

``` r
x <- xts::xts(rgengamma_loc(365, location = 1, scale = 2,
        shape1 = 0.8, shape2 = 0.5),
        order.by = as.Date("2020-01-01") + 0:364)
if (FALSE) { # \dontrun{
fit <- fitlm_gengamma_loc(x, location = 1)
fit$Param
} # }
```
