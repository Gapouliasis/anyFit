# Gamma Distribution

Two-parameter Gamma distribution with scale \\\beta\\ and shape
\\\alpha\\. The distribution is supported on \\x \> 0\\ and is widely
used for modelling positively skewed variables such as precipitation
amounts. Density, distribution, quantile, and random generation
functions are provided via base R
[`dgamma`](https://rdrr.io/r/stats/GammaDist.html). Fitting is performed
via closed-form L-moments through
[`pelgam`](https://rdrr.io/pkg/lmom/man/pel-functions.html).

The probability density function is: \$\$f(x) = \frac{1}{\beta^\alpha
\Gamma(\alpha)} x^{\alpha-1} \exp\left(-\frac{x}{\beta}\right), \quad x
\> 0\$\$ where:

- \\\beta\\ — scale parameter (`scale`)

- \\\alpha\\ — shape parameter (`shape`)

## Usage

``` r
dgamma(x, scale, shape)

pgamma(q, scale, shape)

qgamma(p, scale, shape)

rgamma(n, scale, shape)
```

## Arguments

- x, q:

  vector of quantiles.

- scale:

  scale parameter.

- shape:

  shape parameter.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken as the
  number required.
