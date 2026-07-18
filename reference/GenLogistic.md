# Generalized Logistic Distribution

Three-parameter Generalised Logistic distribution with location \\\mu\\,
scale \\\sigma\\, and shape \\\xi\\. Density, distribution, and quantile
functions are provided via the lmomco package. The distribution
encompasses symmetric (\\\xi = 0\\) and asymmetric forms and is used in
regional frequency analysis of hydrological extremes. Fitting is
performed via closed-form L-moments through
[`pelglo`](https://rdrr.io/pkg/lmom/man/pel-functions.html).

Define \\z = (x - \mu)/\sigma\\. The probability density function is:
\$\$f(x) = \frac{1}{\sigma} \frac{e^{-z}}{(1+e^{-z})^2}, \quad \xi =
0\$\$ \$\$f(x) = \frac{1}{\sigma} \frac{(1-\xi z)^{1/\xi-1}}{(1+(1-\xi
z)^{1/\xi})^2}, \quad \xi \neq 0\$\$ where:

- \\\mu\\ — location parameter (`location`)

- \\\sigma\\ — scale parameter (`scale`)

- \\\xi\\ — shape parameter (`shape`)

## Usage

``` r
dgenlogi(x, location, scale, shape)

pgenlogi(q, location, scale, shape)

qgenlogi(p, location, scale, shape)
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
