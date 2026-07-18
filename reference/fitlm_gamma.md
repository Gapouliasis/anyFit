# fitlm_gamma

Fits the two-parameter Gamma distribution to an xts series by the method
of L-moments. Parameters are obtained in closed form via
[`pelgam`](https://rdrr.io/pkg/lmom/man/pel-functions.html), which
matches the sample L-moments to the theoretical L-moment ratios of the
Gamma distribution. Zero values below `zero_threshold` may be excluded
via `ignore_zeros`. Goodness-of-fit is assessed via `GOF_tests`, and
theoretical L-moments are computed from the fitted parameters for
comparison with the sample.

## Usage

``` r
fitlm_gamma(x, ignore_zeros = FALSE, zero_threshold = 0.01)
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
`scale` and `shape`), `TheorLMom` (theoretical L-moments), `DataLMom`
(sample L-moments), and `GoF` (goodness-of-fit metrics).

## Examples

``` r
x <- xts::xts(rgamma(365, scale = 3, shape = 0.5),
              order.by = as.Date("2020-01-01") + 0:364)
fit <- fitlm_gamma(x)
fit$Param
#> $scale
#> [1] 3.153804
#> 
#> $shape
#> [1] 0.4678799
#> 
```
