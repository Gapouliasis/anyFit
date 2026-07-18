# Generalized Gamma Distribution

Three-parameter Generalised Gamma (Stacy) distribution with scale \\s\\,
first shape \\\alpha\\ (`shape1`), and second shape \\\beta\\
(`shape2`). The distribution is supported on \\x \> 0\\ and contains the
ordinary Gamma (\\\beta = 1\\) and Weibull (\\\beta = \alpha\\) as
special cases. Density, distribution, quantile, and random generation
functions are provided via VGAM using the Stacy parameterisation.
Fitting is performed numerically via L-BFGS-B minimisation of a
normalised L-moment error function, seeded from a pre-computed lookup
table (`GG_InitValues`) with Gauss-Legendre quadrature evaluation of the
L-moments at each optimisation step.

The probability density function is: \$\$f(x) =
\frac{\beta}{\Gamma(\alpha/\beta) \\ s^\alpha} x^{\alpha-1}
\exp\left(-\left(\frac{x}{s}\right)^\beta\right), \quad x \> 0\$\$
where:

- \\s\\ — scale parameter (`scale`)

- \\\alpha\\ — first shape parameter (`shape1`); \\k = \alpha/\beta\\ is
  the Stacy shape

- \\\beta\\ — second shape parameter (`shape2`); the distribution
  reduces to Gamma(\\\alpha\\, \\s\\) when \\\beta=1\\, and to
  Weibull(\\s\\, \\\alpha\\) when \\\beta=\alpha\\

## Usage

``` r
dgengamma(x, scale, shape1, shape2)

pgengamma(q, scale, shape1, shape2)

qgengamma(p, scale, shape1, shape2)

rgengamma(n, scale, shape1, shape2)
```

## Arguments

- x, q:

  vector of quantiles.

- scale:

  scale parameter.

- shape1:

  first shape parameter.

- shape2:

  second shape parameter.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken as the
  number required.
