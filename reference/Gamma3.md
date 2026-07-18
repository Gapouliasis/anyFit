# 3-parameter Gamma Distribution

Three-parameter Gamma (Pearson type III) distribution with location
\\\mu\\, scale \\\beta\\, and shape \\\alpha\\. The distribution is
lower-bounded at the location parameter and is widely employed in
flood-frequency analysis. Density, distribution, quantile, and random
generation functions are provided via the PearsonDS package. Fitting is
performed via closed-form L-moments through
[`pelpe3`](https://rdrr.io/pkg/lmom/man/pel-functions.html).

The probability density function is: \$\$f(x) =
\frac{\beta^\alpha}{\Gamma(\alpha)} (x - \mu)^{\alpha-1}
\exp\left(-\beta (x - \mu)\right), \quad x \ge \mu\$\$ where:

- \\\mu\\ — location parameter (`location`)

- \\\beta\\ — scale parameter (`scale`)

- \\\alpha\\ — shape parameter (`shape`)

## Usage

``` r
dgamma3(x, location, scale, shape)

pgamma3(q, location, scale, shape)

qgamma3(p, location, scale, shape)

rgamma3(n, location, scale, shape)
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
