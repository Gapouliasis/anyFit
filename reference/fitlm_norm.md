# fitlm_norm

Fits the Normal distribution to an xts series by the method of
L-moments. Parameters are obtained in closed form via
[`pelnor`](https://rdrr.io/pkg/lmom/man/pel-functions.html), which maps
the first two sample L-moments directly to the mean and standard
deviation. Zero values below `zero_threshold` may be excluded via
`ignore_zeros`. Goodness-of-fit is assessed via `GOF_tests`, and
theoretical L-moments are computed from the fitted parameters for
comparison with the sample.

## Usage

``` r
fitlm_norm(x, ignore_zeros = FALSE, zero_threshold = 0.01)
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
`mean` and `sd`), `TheorLMom` (theoretical L-moments), `DataLMom`
(sample L-moments), and `GoF` (goodness-of-fit metrics).

## Examples

``` r
x <- xts::xts(rnorm(365, mean = 0, sd = 3),
              order.by = as.Date("2020-01-01") + 0:364)
fit <- fitlm_norm(x)
fit$Param
#> $mean
#> [1] 0.1994077
#> 
#> $sd
#> [1] 3.02824
#> 
```
