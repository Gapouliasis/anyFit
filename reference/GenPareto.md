# Generalized Pareto Distribution

Three-parameter Generalised Pareto distribution (GPD) with location
\\\mu\\, scale \\\sigma\\, and shape \\\xi\\. The GPD is the limiting
distribution of peaks over a threshold and is widely used in
extreme-value modelling. The distribution is lower-bounded at the
location parameter. Density, distribution, quantile, and random
generation functions are provided via the lmom package. Fitting is
performed via closed-form L-moments through
[`pelgpa`](https://rdrr.io/pkg/lmom/man/pel-functions.html); an optional
fixed lower bound is supported.

Define \\z = (x - \mu)/\sigma\\. The probability density function is:
\$\$f(x) = \frac{1}{\sigma} (1 - \xi z)^{1/\xi - 1}, \quad \xi \neq
0\$\$ \$\$f(x) = \frac{1}{\sigma} e^{-z}, \quad \xi = 0\$\$ with \\z \ge
0\\ and \\1 - \xi z \> 0\\. where:

- \\\mu\\ — location parameter (`location`)

- \\\sigma\\ — scale parameter (`scale`)

- \\\xi\\ — shape parameter (`shape`)

## Usage

``` r
dgpd(x, location, scale, shape)

pgpd(q, location, scale, shape)

qgpd(p, location, scale, shape)

rgpd(n, location, scale, shape)
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
