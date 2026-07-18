# Gumbel Distribution

Two-parameter Gumbel (Extreme Value type I) distribution with location
\\\mu\\ and scale \\\sigma\\. The Gumbel distribution is the limiting
distribution of block maxima and is widely used in extreme-value
analysis of hydroclimatic variables. Density, distribution, quantile,
and random generation functions are provided via the FAdist package.
Fitting is performed via closed-form L-moments through
[`pelgum`](https://rdrr.io/pkg/lmom/man/pel-functions.html).

Define \\z = (x - \mu)/\sigma\\. The probability density function is:
\$\$f(x) = \frac{1}{\sigma} \exp(-z - e^{-z})\$\$ where:

- \\\mu\\ — location parameter (`location`)

- \\\sigma\\ — scale parameter (`scale`)

## Usage

``` r
dgumbel(x, location, scale)

pgumbel(q, location, scale)

qgumbel(p, location, scale)

rgumbel(n, location, scale)
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
