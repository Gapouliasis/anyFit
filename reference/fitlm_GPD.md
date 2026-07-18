# fitlm_GPD

Fits the three-parameter Generalised Pareto distribution to an xts
series by the method of L-moments. Parameters are obtained in closed
form via [`pelgpa`](https://rdrr.io/pkg/lmom/man/pel-functions.html),
which matches the sample L-moment ratios to the theoretical L-moment
ratios of the GPD. If `bound` is supplied the location is fixed;
otherwise all three parameters are free. The shape parameter governs
tail behaviour: values greater than -1/2 are required for finite
variance, greater than -1/3 for finite skewness, and greater than -1/4
for finite kurtosis. Zero values below `zero_threshold` may be excluded
via `ignore_zeros`. Goodness-of-fit is assessed via `GOF_tests`.

## Usage

``` r
fitlm_GPD(x, bound = NULL, ignore_zeros = FALSE, zero_threshold = 0.01)
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
x <- xts::xts(rgpd(365, location = 0, scale = 1, shape = -0.2),
              order.by = as.Date("2020-01-01") + 0:364)
fit <- fitlm_GPD(x)
#> [1] "Shape parameter must be >-1/2 for finite variance"
#> [1] "Shape parameter must be >-1/3 for finite skewness"
#> [1] "Shape parameter must be >-1/4 for finite kurtosis"
fit$Param
#> $location
#> [1] -0.01926867
#> 
#> $scale
#> [1] 1.10862
#> 
#> $shape
#> [1] -0.1187803
#> 
```
