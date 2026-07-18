# Dagum Distribution

Three-parameter Dagum (Burr Type III) distribution with scale \\s\\,
first shape \\\alpha\\ (`shape1`), and second shape \\\theta\\
(`shape2`, the tail index). The distribution is supported on \\x \> 0\\
and, unlike the Burr XII, has a unimodal density with mode at the origin
for certain parameter ranges. L-moments are computed via closed-form
Beta-function PWMs; fitting uses L-BFGS-B minimisation seeded from the
`Dagum_InitValues` lookup table with scale-free L-ratio matching and
analytic scale derivation.

The probability density function is: \$\$f(x) = \frac{1}{s}
\left(\frac{x}{s}\right)^{-1/\theta - 1} \left(1 +
\frac{\theta}{\alpha}\left(\frac{x}{s}\right)^{-1/\theta}\right)^{-\alpha-1},
\quad x \> 0\$\$ where:

- \\s\\ — scale parameter (`scale`)

- \\\alpha\\ — first shape parameter (`shape1`), controls upper-tail
  heaviness

- \\\theta\\ — second shape parameter (`shape2`), tail index; the mean
  exists only when \\\theta \< 1\\

## Usage

``` r
pdagum(q, scale, shape1, shape2, PW = 1)

ddagum(x, scale, shape1, shape2, PW = 1)

qdagum(p, scale, shape1, shape2, PW = 1)

rdagum(n, scale, shape1, shape2, PW = 1)
```

## Arguments

- scale:

  scale parameter.

- shape1:

  first shape parameter.

- shape2:

  second shape parameter.

- PW:

  probability weight (point mass at zero for zero-inflated variants).

- x, q:

  vector of quantiles.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken as the
  number required.
