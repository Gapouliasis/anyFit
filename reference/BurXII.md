# Burr Type XII Distribution

Three-parameter Burr Type XII distribution with scale \\s\\, first shape
\\\zeta\\ (`shape1`), and second shape \\\theta\\ (`shape2`, the tail
index). The distribution is supported on \\x \> 0\\ and is widely used
for modelling heavy-tailed positive variables. L-moments are computed
via closed-form expressions involving Beta-function PWMs; a two-step
seeding procedure initialises the shapes from a pre-computed lookup
table (`Burr_InitValues`) and derives the scale analytically. Fitting
uses L-BFGS-B minimisation of the normalised L-moment error.

The probability density function is: \$\$f(x) = \zeta s^{-\zeta}
x^{\zeta-1} \left(\zeta\theta\left(\frac{x}{s}\right)^\zeta +
1\right)^{-1/(\zeta\theta) - 1}, \quad x \> 0\$\$ where:

- \\s\\ — scale parameter (`scale`)

- \\\zeta\\ — first shape parameter (`shape1`), controls lower-tail
  behaviour

- \\\theta\\ — second shape parameter (`shape2`), tail index; the mean
  exists only when \\\theta \< 1\\

## Usage

``` r
dburr(x, scale, shape1, shape2, PW = 1)

pburr(q, scale, shape1, shape2, PW = 1)

qburr(p, scale, shape1, shape2, PW = 1)

rburr(n, scale, shape1, shape2, PW = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- scale:

  scale parameter.

- shape1:

  first shape parameter.

- shape2:

  second shape parameter.

- PW:

  probability weight (point mass at zero for zero-inflated variants).

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken as the
  number required.
