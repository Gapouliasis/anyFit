# fitlm_gev

Fits the three-parameter GEV distribution to an xts series by the method
of L-moments. Parameters are obtained in closed form via
[`pelgev`](https://rdrr.io/pkg/lmom/man/pel-functions.html), which
matches the sample L-moment ratios to the theoretical L-moment ratios of
the GEV distribution. Zero values below `zero_threshold` may be excluded
via `ignore_zeros`. Goodness-of-fit is assessed via `GOF_tests`, and
theoretical L-moments are computed from the fitted parameters for
comparison with the sample.

## Usage

``` r
fitlm_gev(x, ignore_zeros = FALSE, zero_threshold = 0.01)
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

## Value

A list with elements `Distribution`, `Param` (named list of fitted
parameters), `TheorLMom` (theoretical L-moments), `DataLMom` (sample
L-moments), and `GoF` (goodness-of-fit metrics).

## Examples

``` r
x <- xts::xts(rgev(365, location = 0, scale = 1, shape = -0.5),
              order.by = as.Date("2020-01-01") + 0:364)
fit <- fitlm_gev(x)
fit$Param
#> $location
#> [1] -0.02157823
#> 
#> $scale
#> [1] 0.9512176
#> 
#> $shape
#> [1] 0.4586823
#> 
```
