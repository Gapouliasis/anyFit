# Log-Normal Distribution

Three-parameter Log-Normal distribution with location \\\mu\\, scale
\\\beta\\, and shape \\\sigma\\. The distribution is lower-bounded at
\\\mu\\ and arises when the logarithm of \\x - \mu\\ follows a Normal
distribution. Density, distribution, quantile, and random generation
functions are provided via the FAdist package. Fitting is performed via
closed-form L-moments through
[`pelln3`](https://rdrr.io/pkg/lmom/man/pel-functions.html); an optional
fixed lower bound is supported.

The probability density function is: \$\$f(x) = \frac{1}{(x-\mu) \sigma
\sqrt{2\pi}}
\exp\left(-\frac{\left(\ln\left(\frac{x-\mu}{\beta}\right)\right)^2}{2\sigma^2}\right),
\quad x \> \mu\$\$ where:

- \\\mu\\ — location parameter (`location`)

- \\\beta\\ — scale parameter (`scale`)

- \\\sigma\\ — shape parameter (`shape`)

## Usage

``` r
dlognorm(x, location = 0, scale, shape)

plognorm(q, location = 0, scale, shape)

qlognorm(p, location = 0, scale, shape)

rlognorm(n, location = 0, scale, shape)
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
