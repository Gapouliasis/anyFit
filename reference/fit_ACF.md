# fit_ACF

Fits one or more theoretical autocorrelation functions to the sample ACF
of an xts series. Three parametric models are available: the Cauchy-type
autocorrelation structure (CAS, a two-parameter model with algebraic
decay), the Hurst-Kolmogorov process (HK, a one-parameter fractional
Gaussian noise model), and the short-range dependence model (SRD, a
one-parameter exponential decay). Parameters are estimated by minimising
the root-mean-square error between the theoretical and sample ACFs up to
`lag_max`, using the L-BFGS-B optimiser with positivity constraints on
all parameters. The function returns fitted parameter vectors for each
model type, a data frame of theoretical ACF values at each lag, and a
`ggplot` overlay of the fitted curves on the empirical autocorrelation.

## Usage

``` r
fit_ACF(
  ts,
  lag_max,
  type = list("CAS", "HK", "SRD"),
  ignore_zeros = FALSE,
  zero_threshold = 0.01
)
```

## Arguments

- ts:

  An xts object containing the time series data.

- lag_max:

  Integer. Maximum lag to use in fitting.

- type:

  Character vector of ACF model types. Options are `"CAS"`, `"HK"`, and
  `"SRD"`. Default all three.

- ignore_zeros:

  Logical. If `TRUE`, values below `zero_threshold` are excluded.
  Default `FALSE`.

- zero_threshold:

  Numeric. Threshold below which values are treated as zero. Default
  `0.01`.

## Value

A list with elements `ACF_params` (named list of fitted parameter
vectors per model), `ACF_fitted` (data frame of theoretical ACF values
with a `lag` column), and `ACF_plot` (`ggplot` of sample ACF with fitted
curves overlaid).

## Examples

``` r
ts <- xts::xts(rnorm(365), order.by = as.Date("2020-01-01") + 0:364)
acfs <- fit_ACF(ts, lag_max = 10, type = c("CAS", "HK", "SRD"))
acfs$ACF_plot

```
