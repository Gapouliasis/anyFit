# fitlm_lognorm

Fits the three-parameter Log-Normal distribution to an xts series by the
method of L-moments. Parameters are obtained in closed form via
[`pelln3`](https://rdrr.io/pkg/lmom/man/pel-functions.html), which
matches the sample L-moment ratios to the theoretical L-moment ratios of
the Log-Normal distribution. If `bound` is supplied the location is
fixed; otherwise all three parameters are free. Zero values below
`zero_threshold` may be excluded via `ignore_zeros`. Goodness-of-fit is
assessed via `GOF_tests`.

## Usage

``` r
fitlm_lognorm(x, bound = NULL, ignore_zeros = FALSE, zero_threshold = 0.01)
```

## Arguments

- x:

  An xts object containing the time series data.

- bound:

  Numeric or `NULL`. Optional fixed lower bound (location). Default
  `NULL`.

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
x <- xts::xts(rlognorm(365, location = 0, scale = 1, shape = 0.5),
              order.by = as.Date("2020-01-01") + 0:364)
fit <- fitlm_lognorm(x)
fit$Param
#> $location
#> [1] -0.4250635
#> 
#> $scale
#> [1] 1.168126
#> 
#> $shape
#> [1] 0.4228812
#> 
```
