# Vectorised wet/dry transition statistics

Computes per-column transition statistics (mean and variance
after/before zero, mean and variance after non-zero, and transition
probabilities) using lag masks and a sentinel-based counting scheme.

## Usage

``` r
.basic_stats_transition(m, thr)
```

## Arguments

- m:

  A numeric matrix (columns = series, rows = observations).

- thr:

  Numeric threshold defining the zero/non-zero boundary.

## Value

A list with eight named numeric vectors, each of length `ncol(m)`.
