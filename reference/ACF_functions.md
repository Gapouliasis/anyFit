# Theoretical autocorrelation functions

Closed-form expressions for three commonly used theoretical
autocorrelation function (ACF) models: the Cauchy-type algebraic decay
(CAS), the Hurst–Kolmogorov ACF (HK), and the short-range dependence
exponential decay (SRD). Each function returns a data frame of ACF
values for lags 0 through `lag_max`.

Cauchy-type autocorrelation structure with algebraic decay, governed by
two parameters: `kappa` controls the scale and `beta` controls the decay
rate. The ACF decays hyperbolically, producing long-range dependence for
small \\\beta\\.

The ACF is: \$\$\rho(\tau) = (1 + \kappa\beta\tau)^{-1/\beta}\$\$ where:

- \\\kappa\\ — scale parameter (`kappa`)

- \\\beta\\ — shape parameter (`beta`), controls the algebraic decay
  rate

- \\\tau\\ — lag

Hurst–Kolmogorov autocorrelation function, governed by a single Hurst
exponent \\H\\. Values \\H \> 0.5\\ indicate long-range dependence, \\H
= 0.5\\ corresponds to white noise, and \\H \< 0.5\\ indicates
anti-persistence.

The ACF is: \$\$\rho(\tau) = \frac{1}{2}\left(\|\tau-1\|^{2H} -
2\|\tau\|^{2H} + \|\tau+1\|^{2H}\right)\$\$ where:

- \\H\\ — Hurst exponent (`H`), \\0 \< H \< 1\\; \\H \> 0.5\\ indicates
  long-range dependence, \\H = 0.5\\ is white noise, \\H \< 0.5\\
  indicates anti-persistence

- \\\tau\\ — lag

Short-range dependence autocorrelation function with exponential decay,
corresponding to a first-order Markov (AR1) process. Governed by a
single decay parameter \\\kappa\\, the lag-1 autocorrelation is \\\rho_1
= \exp(-\kappa)\\.

The ACF is: \$\$\rho(\tau) = \exp(-\kappa\tau)\$\$ where:

- \\\kappa\\ — decay rate (`kappa`); the lag-1 autocorrelation is
  \\\rho_1 = \exp(-\kappa)\\

- \\\tau\\ — lag

## Usage

``` r
CAS_ACF(kappa, beta, lag_max = 10)

HK_ACF(H, lag_max = 10)

SRD_ACF(kappa, lag_max = 10)
```

## Arguments

- kappa:

  Decay rate. Positive numeric scalar.

- beta:

  Shape parameter. Controls the algebraic decay rate. Positive numeric
  scalar.

- lag_max:

  Maximum lag for which ACF values are computed. Default is 10.

- H:

  Hurst exponent. Numeric scalar in \\(0, 1)\\.

## Value

A data frame with columns `lag` (integer, 0 to `lag_max`) and `ACF`
(numeric, the theoretical autocorrelation at each lag).

A data frame with columns `lag` (integer, 0 to `lag_max`) and `ACF`
(numeric, the theoretical autocorrelation at each lag).

A data frame with columns `lag` (integer, 0 to `lag_max`) and `ACF`
(numeric, the theoretical autocorrelation at each lag).

## Examples

``` r
# Cauchy-type ACF with slow decay
CAS_ACF(kappa = 0.5, beta = 0.3)
#>    lag        ACF
#> 1    0 1.00000000
#> 2    1 0.62758689
#> 3    2 0.41705067
#> 4    3 0.28980552
#> 5    4 0.20873730
#> 6    5 0.15483644
#> 7    6 0.11771216
#> 8    7 0.09137354
#> 9    8 0.07220896
#> 10   9 0.05795716
#> 11  10 0.04715560

# Faster decay with larger beta
CAS_ACF(kappa = 0.5, beta = 1.5, lag_max = 20)
#>    lag       ACF
#> 1    0 1.0000000
#> 2    1 0.6886121
#> 3    2 0.5428835
#> 4    3 0.4557686
#> 5    4 0.3968503
#> 6    5 0.3538921
#> 7    6 0.3209408
#> 8    7 0.2947225
#> 9    8 0.2732759
#> 10   9 0.2553478
#> 11  10 0.2400974
#> 12  11 0.2269371
#> 13  12 0.2154435
#> 14  13 0.2053026
#> 15  14 0.1962764
#> 16  15 0.1881811
#> 17  16 0.1808719
#> 18  17 0.1742334
#> 19  18 0.1681724
#> 20  19 0.1626123
#> 21  20 0.1574901

# Long-range dependence (H > 0.5)
HK_ACF(H = 0.8)
#>    lag       ACF
#> 1    0 1.0000000
#> 2    1 0.5157166
#> 3    2 0.3683399
#> 4    3 0.3109639
#> 5    4 0.2765057
#> 6    5 0.2526226
#> 7    6 0.2347187
#> 8    7 0.2206062
#> 9    8 0.2090851
#> 10   9 0.1994322
#> 11  10 0.1911809

# Anti-persistence (H < 0.5)
HK_ACF(H = 0.3, lag_max = 20)
#>    lag          ACF
#> 1    0  1.000000000
#> 2    1 -0.242141717
#> 3    2 -0.049125544
#> 4    3 -0.026625407
#> 5    4 -0.017541785
#> 6    5 -0.012751424
#> 7    6 -0.009844225
#> 8    7 -0.007916697
#> 9    8 -0.006557919
#> 10   9 -0.005555839
#> 11  10 -0.004790730
#> 12  11 -0.004190245
#> 13  12 -0.003708293
#> 14  13 -0.003314222
#> 15  14 -0.002986920
#> 16  15 -0.002711408
#> 17  16 -0.002476789
#> 18  17 -0.002274961
#> 19  18 -0.002099788
#> 20  19 -0.001946540
#> 21  20 -0.001811522

# Moderate short-range dependence
SRD_ACF(kappa = 0.5)
#>    lag         ACF
#> 1    0 1.000000000
#> 2    1 0.606530660
#> 3    2 0.367879441
#> 4    3 0.223130160
#> 5    4 0.135335283
#> 6    5 0.082084999
#> 7    6 0.049787068
#> 8    7 0.030197383
#> 9    8 0.018315639
#> 10   9 0.011108997
#> 11  10 0.006737947

# Very short memory
SRD_ACF(kappa = 2, lag_max = 15)
#>    lag          ACF
#> 1    0 1.000000e+00
#> 2    1 1.353353e-01
#> 3    2 1.831564e-02
#> 4    3 2.478752e-03
#> 5    4 3.354626e-04
#> 6    5 4.539993e-05
#> 7    6 6.144212e-06
#> 8    7 8.315287e-07
#> 9    8 1.125352e-07
#> 10   9 1.522998e-08
#> 11  10 2.061154e-09
#> 12  11 2.789468e-10
#> 13  12 3.775135e-11
#> 14  13 5.109089e-12
#> 15  14 6.914400e-13
#> 16  15 9.357623e-14
```
