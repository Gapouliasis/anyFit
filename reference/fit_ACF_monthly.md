# fit_ACF_monthly

Wraps
[`fit_ACF`](https://gapouliasis.github.io/anyFit/reference/fit_ACF.md)
to fit ACF models separately to each calendar-month sub-series of an xts
object. The input series is split by calendar month, `fit_ACF` is
applied to each resulting sub-series, and the fitted parameters are
assembled into a parameter matrix with months as rows. A 12-panel
patchwork plot of fitted ACF curves is produced, with empty panels for
months absent from the data. The grid layout is controlled by `nrow` and
`ncol`. The three model types (CAS, HK, SRD) supported by `fit_ACF` are
available and may be specified singly or in combination.

## Usage

``` r
fit_ACF_monthly(ts, lag, type = list("CAS", "HK", "SRD"), nrow = 4, ncol = 3)
```

## Arguments

- ts:

  An xts object containing the time series data.

- lag:

  Integer. Maximum lag passed to `fit_ACF`.

- type:

  Character vector of ACF model types. Options are `"CAS"`, `"HK"`, and
  `"SRD"`. Default all three.

- nrow:

  Integer. Number of rows in the panel grid. Default `4`.

- ncol:

  Integer. Number of columns in the panel grid. Default `3`.

## Value

A list with elements `ACF_params_monthly` (matrix of fitted parameters
with months as rows) and `ACF_monthly_plot` (12-panel patchwork of
fitted ACF `ggplot` objects with a shared legend).

## Examples

``` r
ts <- xts::xts(rnorm(365), order.by = as.Date("2020-01-01") + 0:364)
monthly_acfs <- fit_ACF_monthly(ts, lag = 10)
monthly_acfs$ACF_monthly_plot

```
