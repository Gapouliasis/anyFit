# fitlm_gamma3

Fits the three-parameter Gamma (Pearson type III) distribution to an xts
series by the method of L-moments. Parameters are obtained via
[`pelpe3`](https://rdrr.io/pkg/lmom/man/pel-functions.html) and then
transformed to the standard location, scale, and shape parameterisation
used throughout the package. Zero values below `zero_threshold` may be
excluded via `ignore_zeros`. Goodness-of-fit is assessed via
`GOF_tests`, and theoretical L-moments are computed from the fitted
parameters.

## Usage

``` r
fitlm_gamma3(x, ignore_zeros = FALSE, zero_threshold = 0.01)
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
x <- xts::xts(rgamma3(365, location = 5, scale = 0.5, shape = 3),
              order.by = as.Date("2020-01-01") + 0:364)
fit <- fitlm_gamma3(x)
fit$Param
#> $location
#> [1] 3.952122
#> 
#> $scale
#> [1] 0.6040903
#> 
#> $shape
#> [1] 4.269671
#> 
```
