# Exponential Weibull Distribution

Three-parameter Exponentiated Weibull distribution with scale \\s\\,
first shape \\a\\ (`shape1`, the Weibull shape), and second shape \\b\\
(`shape2`, the exponentiation power). The distribution is supported on
\\x \> 0\\ and contains the ordinary Weibull (\\b = 1\\) and Exponential
(\\a = 1, b = 1\\) as special cases. L-moments are computed by tanh-sinh
(double-exponential) quadrature, which is well-conditioned across the
full positive parameter range. Fitting uses L-BFGS-B minimisation seeded
from the `ExpWeibull_InitValues` lookup table with scale-free L-ratio
matching and analytic scale derivation.

The probability density function is: \$\$f(x) = \frac{b a}{s}
\left(\frac{x}{s}\right)^{a-1}
\exp\left(-\left(\frac{x}{s}\right)^a\right) \left(1 -
\exp\left(-\left(\frac{x}{s}\right)^a\right)\right)^{b-1}, \quad x \>
0\$\$ where:

- \\s\\ — scale parameter (`scale`)

- \\a\\ — first shape parameter (`shape1`), the Weibull shape; a = 1
  gives the exponential baseline

- \\b\\ — second shape parameter (`shape2`), the exponentiation power; b
  = 1 recovers the ordinary Weibull

## Usage

``` r
pexpweibull(q, scale, shape1, shape2, log.p = FALSE)

dexpweibull(x, scale, shape1, shape2, log = FALSE)

qexpweibull(p, scale, shape1, shape2)

rexpweibull(n, scale, shape1, shape2)
```

## Arguments

- scale:

  scale parameter.

- shape1:

  first shape parameter.

- shape2:

  second shape parameter.

- x, q:

  vector of quantiles.

- log, log.p:

  logical; if `TRUE`, probabilities p are given as log(p).

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken as the
  number required.
