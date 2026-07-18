# Generalized Gamma with Location Distribution

Four-parameter Generalised Gamma distribution with location \\\mu\\,
scale \\s\\, first shape \\\alpha\\, and second shape \\\beta\\. The
distribution is lower-bounded at the location parameter and extends the
three-parameter Stacy form by a shift. Fitting is delegated to
[`fitlm_gengamma`](https://gapouliasis.github.io/anyFit/reference/fitlm_gengamma.md)
after subtracting the location from the data; the location is supplied
by the user rather than fitted automatically. Density, distribution,
quantile, and random generation functions are provided via VGAM.

The probability density function is: \$\$f(x) =
\frac{\beta}{\Gamma(\alpha/\beta) \\ s^\alpha} (x-\mu)^{\alpha-1}
\exp\left(-\left(\frac{x-\mu}{s}\right)^\beta\right), \quad x \ge
\mu\$\$ where:

- \\\mu\\ — location parameter (`location`)

- \\s\\ — scale parameter (`scale`)

- \\\alpha\\ — first shape parameter (`shape1`)

- \\\beta\\ — second shape parameter (`shape2`)

## Usage

``` r
pgengamma_loc(q, location, scale, shape1, shape2)

dgengamma_loc(x, location, scale, shape1, shape2)

qgengamma_loc(p, location, scale, shape1, shape2)

rgengamma_loc(n, location, scale, shape1, shape2)
```

## Arguments

- location:

  location parameter.

- scale:

  scale parameter.

- shape1:

  first shape parameter.

- shape2:

  second shape parameter.

- x, q:

  vector of quantiles.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken as the
  number required.
