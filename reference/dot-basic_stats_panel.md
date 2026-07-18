# Build the diagnostic panel for a single series

Constructs a four-panel `patchwork` plot (raw time series, PDF histogram
with density overlay, empirical CDF on log-probability scale, and sample
ACF) for a single-column xts series.

## Usage

``` r
.basic_stats_panel(
  ts,
  pstart,
  pend,
  show_label,
  label_prefix,
  show_table,
  xpos_label,
  ypos_label,
  xpos_table,
  ypos_table,
  nbins,
  ignore_zeros,
  zero_threshold
)
```

## Arguments

- ts:

  A single-column xts object.

- pstart:

  Plot start date (`NA` for series start).

- pend:

  Plot end date (`NA` for series end).

- show_label:

  Logical. Overlay series label.

- label_prefix:

  Character prefix for the label.

- show_table:

  Logical. Overlay statistics table.

- xpos_label:

  Horizontal label position (0–1 fraction).

- ypos_label:

  Vertical label position (0–1 fraction).

- xpos_table:

  Horizontal table position (0–1 fraction).

- ypos_table:

  Vertical table position (0–1 fraction).

- nbins:

  Number of PDF histogram bins.

- ignore_zeros:

  Logical. Exclude sub-threshold values.

- zero_threshold:

  Numeric zero threshold.

## Value

A `patchwork` object combining four ggplot panels.
