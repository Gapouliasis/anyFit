# Normal Distribution

Two-parameter Normal (Gaussian) distribution with mean \\\mu\\ and
standard deviation \\\sigma\\. Fitting is performed via closed-form
L-moments. The L-moment estimates coincide with the conventional
method-of-moments estimates (`mean` = \\\lambda_1\\, `sd` = \\\lambda_2
\sqrt{\pi}\\).

The probability density function is: \$\$f(x) =
\frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{1}{2}\left(\frac{x -
\mu}{\sigma}\right)^2\right)\$\$ where:

- \\\mu\\ — location/mean parameter (`mean`)

- \\\sigma\\ — scale/standard deviation parameter (`sd`)

## Usage

``` r
dnorm(x, mean = 0, sd = 1)

pnorm(q, mean = 0, sd = 1)

qnorm(p, mean = 0, sd = 1)

rnorm(n, mean = 0, sd = 1)
```

## Arguments

- x, q:

  vector of quantiles.

- mean:

  mean parameter.

- sd:

  standard deviation parameter.

- p:

  vector of probabilities.

- n:

  number of observations. If `length(n) > 1`, the length is taken as the
  number required.
