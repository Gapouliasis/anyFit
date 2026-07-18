# Rayleigh Distribution

Two-parameter Rayleigh distribution with location \\\mu\\ and scale
\\\sigma\\. The distribution is lower-bounded at the location parameter.
Fitting uses closed-form L-moment expressions derived from
gamma-function relationships; the location is obtained from the first
two L-moments and the scale from the second L-moment ratio.

The probability density function is: \$\$f(x) = \frac{x - \mu}{\sigma^2}
\exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right), \quad x \ge \mu\$\$
where:

- \\\mu\\ — location parameter (`location`)

- \\\sigma\\ — scale parameter (`scale`)

## Usage

``` r
drayleigh(x, location, scale)

prayleigh(x, location, scale)

qrayleigh(p, location, scale)

rrayleigh(n, location, scale)
```

## Arguments

- x:

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
