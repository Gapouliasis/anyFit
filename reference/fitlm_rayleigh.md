# fitlm_rayleigh

Fits the two-parameter Rayleigh distribution to an xts series by the
method of L-moments. Closed-form expressions derived from gamma-function
integrals are used: the location is recovered from \\\lambda_1\\ and
\\\lambda_2\\, and the scale follows from the squared second L-moment
ratio. Zero values below `zero_threshold` may be excluded via
`ignore_zeros`. Goodness-of-fit is assessed via `GOF_tests`.

## Usage

``` r
fitlm_rayleigh(x, ignore_zeros = FALSE, zero_threshold = 0.01)
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
x <- xts::xts(rrayleigh(365, location = 0, scale = 2),
              order.by = as.Date("2020-01-01") + 0:364)
fit <- fitlm_rayleigh(x)
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
#> Warning: longer object length is not a multiple of shorter object length
fit$Param
#> $location
#>        l_1        l_2 
#> -6.1572535  0.1673074 
#> 
#> $scale
#>      l_1      l_2 
#> 6.947717 1.901447 
#> 
```
