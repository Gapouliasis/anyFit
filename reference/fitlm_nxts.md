# Fit Distributions to Multiple Time Series

Fits a set of candidate probability distributions to every column of an
xts object using the L-moment method. The function wraps
[`fitlm_multi`](https://gapouliasis.github.io/anyFit/reference/fitlm_multi.md)
for each series and collects the estimated parameters, goodness-of-fit
metrics, and diagnostic plots (Q-Q, P-P) into a list. When
`parallel = TRUE`, the work is distributed across available cores via
the future framework, with either file-backed shared memory (bigmemory
memory-mapped matrices for single-machine parallelism) or serialised
column chunks (for multi-node clusters). File-backed shared memory can
be enabled by `shared_memory = TRUE`. This is reccomended for single
machine usage and especially windows for efficiency and reduced RAM
consumption. Panel layouts of the per-series Q-Q and P-P plots are
assembled via patchwork in pages of `nrow` by `ncol` sub-plots. This
function is designed for large ensembles of grid cells or stations where
fitting individually would be impractical.

## Usage

``` r
fitlm_nxts(
  ts,
  candidates,
  nrow = 5,
  ncol = 4,
  ignore_zeros = FALSE,
  zero_threshold = 0.01,
  parallel = FALSE,
  ncores = 2,
  shared_memory = TRUE,
  diagnostic_plots = TRUE,
  order = NULL
)
```

## Arguments

- ts:

  An xts object containing the time series data (multiple columns).

- candidates:

  A list of distribution names to fit (e.g.
  `list('exp','gamma3','gengamma')`).

- nrow:

  Number of rows per panel page in the diagnostic-plot grid.

- ncol:

  Number of columns per panel page in the diagnostic-plot grid.

- ignore_zeros:

  Logical. If `TRUE`, zeros are excluded before fitting. Default
  `FALSE`.

- zero_threshold:

  Numeric. Values below this threshold are treated as zero. Default
  0.01.

- parallel:

  Logical. If `TRUE`, use parallel processing via future. Default
  `FALSE`.

- ncores:

  Number of worker processes when `parallel = TRUE` and no pre-existing
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  is set. Default 2.

- shared_memory:

  Logical. When `parallel = TRUE`, share the data grid with workers via
  a file-backed bigmatrix (memory-mapped, single machine only). Set
  `FALSE` for multi-node `plan(cluster)` setups. Default `TRUE`.

- diagnostic_plots:

  Logical. If `TRUE`, produce Q-Q and P-P diagnostic plots. Default
  `TRUE`.

- order:

  Optional named list mapping candidate distribution names to the vector
  of L-moment orders used by their optimiser, e.g.
  `list(gengamma = 1:5, expweibull = 1:3)`. Passed through to
  [`fitlm_multi`](https://gapouliasis.github.io/anyFit/reference/fitlm_multi.md).
  Default `NULL`.

## Value

A list with elements `params` (per-column parameter tables and GoF
summaries), and, if `diagnostic_plots = TRUE`, `diagnostic_plots`,
`QQ_plots`, `PP_plots`, `QQ_panels`, and `PP_panels`.

## Examples

``` r
# Synthetic daily data: 3 stations, 10 years
set.seed(123)
n <- 3650
ts <- xts::xts(cbind(station1 = rgamma(n, shape = 2, scale = 5),
                station2 = rgamma(n, shape = 3, scale = 7),
                station3 = rgamma(n, shape = 2.5, scale = 6)),
          order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = n))

fits <- fitlm_nxts(ts, candidates = list('exp','gamma3'),
                   nrow = 2, ncol = 2, ignore_zeros = FALSE)
fits$params[[1]]
#> $params
#> $params$exp
#> $params$exp$Distribution
#> $params$exp$Distribution$FXs
#> [1] "qexp"
#> 
#> $params$exp$Distribution$bound
#> NULL
#> 
#> 
#> $params$exp$Param
#> $params$exp$Param$location
#> [1] 2.618052
#> 
#> $params$exp$Param$scale
#> [1] 7.176231
#> 
#> 
#> $params$exp$TheorLMom
#>  lambda_1  lambda_2     tau_3     tau_4 
#> 9.7942835 3.5881157 0.3333333 0.1666667 
#> 
#> $params$exp$DataLMom
#>  lambda_1  lambda_2     tau_3     tau_4 
#> 9.7942835 3.5881157 0.2257197 0.1434131 
#> 
#> $params$exp$GoF
#> $params$exp$GoF$MLE
#> [1] 14975.7
#> 
#> $params$exp$GoF$CM
#> [1] 2.537503
#> 
#> $params$exp$GoF$KS
#> [1] 0.09507189
#> 
#> $params$exp$GoF$MSEquant
#> [1] 0.9175677
#> 
#> $params$exp$GoF$DiffOfMax
#> [1] 16.60772
#> 
#> $params$exp$GoF$MeanDiffOf10Max
#> [1] 6.04083
#> 
#> 
#> 
#> $params$gamma3
#> $params$gamma3$Distribution
#> $params$gamma3$Distribution$FXs
#> [1] "qgamma3"
#> 
#> 
#> $params$gamma3$Param
#> $params$gamma3$Param$location
#> [1] -0.09710559
#> 
#> $params$gamma3$Param$scale
#> [1] 0.2179958
#> 
#> $params$gamma3$Param$shape
#> [1] 2.156281
#> 
#> 
#> $params$gamma3$TheorLMom
#>  lambda_1  lambda_2     tau_3     tau_4 
#> 9.7942835 3.5881158 0.2257190 0.1399896 
#> 
#> $params$gamma3$DataLMom
#>       l_1       l_2       t_3       t_4 
#> 9.7942835 3.5881157 0.2257197 0.1434131 
#> 
#> $params$gamma3$GoF
#> $params$gamma3$GoF$MLE
#> [1] 11513.06
#> 
#> $params$gamma3$GoF$CM
#> [1] 0.008314168
#> 
#> $params$gamma3$GoF$KS
#> [1] 0.00736687
#> 
#> $params$gamma3$GoF$MSEquant
#> [1] 0.02230082
#> 
#> $params$gamma3$GoF$DiffOfMax
#> [1] -4.455013
#> 
#> $params$gamma3$GoF$MeanDiffOf10Max
#> [1] 1.97211
#> 
#> 
#> 
#> 
#> $GoF_summary
#>                          exp        gamma3
#> MLE             1.497570e+04  1.151306e+04
#> CM              2.537503e+00  8.314168e-03
#> KS              9.507189e-02  7.366870e-03
#> MSEquant        9.175677e-01  2.230082e-02
#> DiffOfMax       1.660772e+01 -4.455013e+00
#> MeanDiffOf10Max 6.040830e+00  1.972110e+00
#> 
```
