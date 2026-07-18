# fitlm_exp

Fits the two-parameter Exponential distribution to an xts series by the
method of L-moments. If `bound` is supplied, the location is fixed to
that value and the scale is derived from the second L-moment; otherwise
both parameters are obtained from
[`pelexp`](https://rdrr.io/pkg/lmom/man/pel-functions.html). Zero values
below `zero_threshold` may be excluded via `ignore_zeros`.
Goodness-of-fit is assessed via `GOF_tests`, and theoretical L-moments
are computed from the fitted parameters for comparison with the sample.

## Usage

``` r
fitlm_exp(x, bound = NULL, ignore_zeros = FALSE, zero_threshold = 0.01)
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
x <- xts::xts(rexp(365, location = 0, scale = 10), order.by = as.Date("2020-01-01") + 0:364)
fit <- fitlm_exp(x)
fit$Param
#> $location
#> [1] 0.1021415
#> 
#> $scale
#> [1] 9.977176
#> 
```
