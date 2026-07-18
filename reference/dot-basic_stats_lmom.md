# Per-column L-moments

Computes L-mean, L-scale, L3, L4, L-CV, L-skewness, and L-kurtosis for a
single numeric vector via
[`samlmu`](https://rdrr.io/pkg/lmom/man/samlmu.html).

## Usage

``` r
.basic_stats_lmom(x)
```

## Arguments

- x:

  A numeric vector.

## Value

A numeric vector of length 7: LMean, LScale, L3, L4, L-CV, L-Skewness,
L-Kurtosis. Returns `NA` for all elements on error or if `x` is empty.
