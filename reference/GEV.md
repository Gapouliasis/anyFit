# Generalized Extreme Value Distribution

Three-parameter Generalised Extreme Value (GEV) distribution with
location \\\mu\\, scale \\\sigma\\, and shape \\\xi\\. The GEV
encompasses the Gumbel (\\\xi = 0\\), Frechet (\\\xi \> 0\\), and
reversed Weibull (\\\xi \< 0\\) families. Density, distribution,
quantile, and random generation functions are implemented natively.
Fitting is performed via closed-form L-moments through
[`pelgev`](https://rdrr.io/pkg/lmom/man/pel-functions.html).

Define \\t(x) = (1 + \xi(x-\mu)/\sigma)^{-1/\xi}\\. The probability
density function is: \$\$f(x) = \frac{1}{\sigma} t(x)^{\xi+1}
\exp(-t(x)), \quad \xi \neq 0\$\$ \$\$f(x) = \frac{1}{\sigma}
\exp\left(-\frac{x-\mu}{\sigma} -
\exp\left(-\frac{x-\mu}{\sigma}\right)\right), \quad \xi = 0\$\$ with
\\1 + \xi(x-\mu)/\sigma \> 0\\. where:

- \\\mu\\ — location parameter (`location`)

- \\\sigma\\ — scale parameter (`scale`)

- \\\xi\\ — shape parameter (`shape`)

## Usage

``` r
dgev(x, location, scale, shape)

pgev(q, location, scale, shape)

qgev(p, shape, scale, location)

rgev(n, location, scale, shape)
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
