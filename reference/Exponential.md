# Exponential Distribution

Two-parameter Exponential distribution with location \\\mu\\ and scale
\\\sigma\\. The distribution is lower-bounded at the location parameter.
Fitting is performed via closed-form L-moments. An optional fixed lower
bound `bound` is supported for cases where the minimum is known a
priori.

The probability density function is: \$\$f(x) = \frac{1}{\sigma}
\exp\left(-\frac{x-\mu}{\sigma}\right), \quad x \ge \mu\$\$ where:

- \\\mu\\ — location parameter (`location`)

- \\\sigma\\ — scale parameter (`scale`)

## Usage

``` r
dexp(x, location, scale)

pexp(q, location, scale)

qexp(p, location, scale)

rexp(n, location, scale)
```

## Arguments

- x, q:

  vector of quantiles.

- location:

  location parameter.

- scale:

  scale parameter.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken as the
  number required.
