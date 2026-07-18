# 3-parameter Weibull Distribution

Three-parameter Weibull distribution with location \\\mu\\, scale
\\\sigma\\, and shape \\k\\. The distribution is lower-bounded at the
location parameter and generalises the Exponential (\\k = 1\\) and
Rayleigh (\\k = 2\\) distributions. Fitting is performed via closed-form
L-moments. An optional fixed lower bound is supported.

The probability density function is: \$\$f(x) = \frac{k}{\sigma}
\left(\frac{x - \mu}{\sigma}\right)^{k-1} \exp\left(-\left(\frac{x -
\mu}{\sigma}\right)^k\right), \quad x \ge \mu\$\$ where:

- \\\mu\\ — location parameter (`location`)

- \\\sigma\\ — scale parameter (`scale`)

- \\k\\ — shape parameter (`shape`)

## Usage

``` r
dweibull(x, location, scale, shape)

pweibull(q, location, scale, shape)

qweibull(p, location, scale, shape)

rweibull(n, location, scale, shape)
```

## Arguments

- x, q:

  vector of quantiles.

- location:

  location parameter.

- scale:

  scale parameter.

- shape:

  shape parameter.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken as the
  number required.
