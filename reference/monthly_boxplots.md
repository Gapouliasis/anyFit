# Monthly Boxplots

Produces a ggplot boxplot panel grouped by calendar month. Multi-column
xts input is faceted by variable; single-column input is coloured by
month.

## Usage

``` r
monthly_boxplots(
  ts,
  palette = "Set1",
  ignore_zeros = FALSE,
  zero_threshold = 0.01
)
```

## Arguments

- ts:

  An xts object containing the time series data.

- palette:

  An RColorBrewer palette name. Default `'Set1'`.

- ignore_zeros:

  Logical. If `TRUE`, zeros are excluded. Default `FALSE`.

- zero_threshold:

  Numeric. Values below this threshold are treated as zero. Default
  0.01.

## Value

A ggplot object with monthly boxplots.

## Examples

``` r
# Synthetic daily data
set.seed(42)
ts <- xts::xts(rgamma(3650, shape = 2, scale = 5),
          order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = 3650))
boxes <- monthly_boxplots(ts, palette = 'Set3', ignore_zeros = TRUE)
```
