# Package index

## Loading & conversion

Read point and gridded data into xts / sxts.

- [`nc2xts()`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md)
  : nc2xts
- [`nc2xts_nn()`](https://gapouliasis.github.io/anyFit/reference/nc2xts_nn.md)
  : nc2xts_nn
- [`delim2xts()`](https://gapouliasis.github.io/anyFit/reference/delim2xts.md)
  : delim2xts

## The sxts class & spatial operations

The spatial-xts class, its accessors, and native spatial operations.

- [`sxts()`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`str(`*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`summary(`*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`print(`*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`attributes.sxts()`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`coords(`*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`projection(`*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`elements(`*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`` `[`( ``*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`is.sxts()`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`lag(`*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`diff(`*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`Ops(`*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`rasterFromSxts(`*`<sxts>`*`)`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  [`sxtsFromRaster()`](https://gapouliasis.github.io/anyFit/reference/sxts.md)
  : Spatial xts (sxts) class
- [`coords()`](https://gapouliasis.github.io/anyFit/reference/coords.md)
  : Spatial coordinates of an sxts object
- [`projection()`](https://gapouliasis.github.io/anyFit/reference/projection.md)
  : Coordinate reference system of an sxts object
- [`elements()`](https://gapouliasis.github.io/anyFit/reference/elements.md)
  : Number of spatial locations in an sxts object
- [`mask.sxts()`](https://gapouliasis.github.io/anyFit/reference/mask.sxts.md)
  : Spatial masking of sxts objects
- [`zonal_stats()`](https://gapouliasis.github.io/anyFit/reference/zonal_stats.md)
  : Aggregate sxts over shapefile polygons (zonal statistics)
- [`rasterFromSxts()`](https://gapouliasis.github.io/anyFit/reference/rasterFromSxts.md)
  : Convert sxts to RasterBrick
- [`rbindlist.sxts()`](https://gapouliasis.github.io/anyFit/reference/rbindlist.sxts.md)
  : Fast row-bind of a list of sxts objects

## Summary statistics

Per-series and per-cell summary statistics, including L-moments.

- [`basic_stats()`](https://gapouliasis.github.io/anyFit/reference/basic_stats.md)
  : Basic statistics and diagnostic panel for time series data
- [`basic_stats_nc()`](https://gapouliasis.github.io/anyFit/reference/basic_stats_nc.md)
  : Compute basic statistics for gridded (NetCDF or raster) data
- [`lmom_stats()`](https://gapouliasis.github.io/anyFit/reference/lmom_stats.md)
  : L-Moment Statistics
- [`period_stats()`](https://gapouliasis.github.io/anyFit/reference/period_stats.md)
  : Period-based summary statistics
- [`monthly_stats()`](https://gapouliasis.github.io/anyFit/reference/monthly_stats.md)
  : monthly_stats
- [`monthly_stats_nc()`](https://gapouliasis.github.io/anyFit/reference/monthly_stats_nc.md)
  : monthly_stats_nc
- [`check_missing()`](https://gapouliasis.github.io/anyFit/reference/check_missing.md)
  : Check for missing values in time series data by period

## Aggregation & transforms

Temporal aggregation and value transformations.

- [`aggregate_xts()`](https://gapouliasis.github.io/anyFit/reference/aggregate_xts.md)
  : Aggregate an xts time series to coarser temporal scales
- [`period_apply_nc()`](https://gapouliasis.github.io/anyFit/reference/period_apply_nc.md)
  : Temporal aggregation of gridded data
- [`normalise_xts()`](https://gapouliasis.github.io/anyFit/reference/normalise_xts.md)
  : normalise_xts

## Distribution fitting

L-moment fitting, goodness of fit, and L-ratio diagnostics.

- [`fitlm_multi()`](https://gapouliasis.github.io/anyFit/reference/fitlm_multi.md)
  : Multi-Distribution Fitting
- [`fitlm_nc()`](https://gapouliasis.github.io/anyFit/reference/fitlm_nc.md)
  : Distribution Fitting on Gridded Data
- [`fitlm_nxts()`](https://gapouliasis.github.io/anyFit/reference/fitlm_nxts.md)
  : Fit Distributions to Multiple Time Series
- [`fitlm_monthly()`](https://gapouliasis.github.io/anyFit/reference/fitlm_monthly.md)
  : Monthly Distribution Fitting
- [`fitlm_monthly_nc()`](https://gapouliasis.github.io/anyFit/reference/fitlm_monthly_nc.md)
  : Monthly Distribution Fitting on Gridded Data
- [`fitlm_exp()`](https://gapouliasis.github.io/anyFit/reference/fitlm_exp.md)
  : fitlm_exp
- [`fitlm_rayleigh()`](https://gapouliasis.github.io/anyFit/reference/fitlm_rayleigh.md)
  : fitlm_rayleigh
- [`fitlm_gamma()`](https://gapouliasis.github.io/anyFit/reference/fitlm_gamma.md)
  : fitlm_gamma
- [`fitlm_gamma3()`](https://gapouliasis.github.io/anyFit/reference/fitlm_gamma3.md)
  : fitlm_gamma3
- [`fitlm_genlogi()`](https://gapouliasis.github.io/anyFit/reference/fitlm_genlogi.md)
  : fitlm_genlogi
- [`fitlm_norm()`](https://gapouliasis.github.io/anyFit/reference/fitlm_norm.md)
  : fitlm_norm
- [`fitlm_weibull()`](https://gapouliasis.github.io/anyFit/reference/fitlm_weibull.md)
  : fitlm_weibull
- [`fitlm_gumbel()`](https://gapouliasis.github.io/anyFit/reference/fitlm_gumbel.md)
  : fitlm_gumbel
- [`fitlm_lognorm()`](https://gapouliasis.github.io/anyFit/reference/fitlm_lognorm.md)
  : fitlm_lognorm
- [`fitlm_gev()`](https://gapouliasis.github.io/anyFit/reference/fitlm_gev.md)
  : fitlm_gev
- [`fitlm_GPD()`](https://gapouliasis.github.io/anyFit/reference/fitlm_GPD.md)
  : fitlm_GPD
- [`fitlm_gengamma()`](https://gapouliasis.github.io/anyFit/reference/fitlm_gengamma.md)
  : fitlm_gengamma
- [`fitlm_gengamma_loc()`](https://gapouliasis.github.io/anyFit/reference/fitlm_gengamma_loc.md)
  : fitlm_gengamma_loc
- [`fitlm_burr()`](https://gapouliasis.github.io/anyFit/reference/fitlm_burr.md)
  : fitlm_burr
- [`fitlm_dagum()`](https://gapouliasis.github.io/anyFit/reference/fitlm_dagum.md)
  : fitlm_dagum
- [`fitlm_expweibull()`](https://gapouliasis.github.io/anyFit/reference/fitlm_expweibull.md)
  : fitlm_expweibull
- [`LRatio_check()`](https://gapouliasis.github.io/anyFit/reference/LRatio_check.md)
  : Check Sample L-Ratios Against Distribution Support
- [`fit_diagnostics()`](https://gapouliasis.github.io/anyFit/reference/fit_diagnostics.md)
  : Diagnostics for a Fitted Distribution

## Autocorrelation

Sample and theoretical autocorrelation structures (CAS, HK, SRD).

- [`fit_ACF()`](https://gapouliasis.github.io/anyFit/reference/fit_ACF.md)
  : fit_ACF
- [`fit_ACF_monthly()`](https://gapouliasis.github.io/anyFit/reference/fit_ACF_monthly.md)
  : fit_ACF_monthly
- [`CAS_ACF()`](https://gapouliasis.github.io/anyFit/reference/ACF_functions.md)
  [`HK_ACF()`](https://gapouliasis.github.io/anyFit/reference/ACF_functions.md)
  [`SRD_ACF()`](https://gapouliasis.github.io/anyFit/reference/ACF_functions.md)
  : Theoretical autocorrelation functions

## Visualization

ggplot2-based maps, diagnostics and seasonal plots.

- [`nc_ggplot()`](https://gapouliasis.github.io/anyFit/reference/nc_ggplot.md)
  : nc_ggplot
- [`correl_plots()`](https://gapouliasis.github.io/anyFit/reference/correl_plots.md)
  : correl_plots
- [`monthly_boxplots()`](https://gapouliasis.github.io/anyFit/reference/monthly_boxplots.md)
  : Monthly Boxplots
- [`monthly_violins()`](https://gapouliasis.github.io/anyFit/reference/monthly_violins.md)
  : monthly_violins
- [`monthly_ecdf()`](https://gapouliasis.github.io/anyFit/reference/monthly_ecdf.md)
  : Monthly Empirical CDF Grid
- [`monthly_fdiagnostics()`](https://gapouliasis.github.io/anyFit/reference/monthly_fdiagnostics.md)
  : Monthly Distribution Diagnostics
- [`ridge_plots()`](https://gapouliasis.github.io/anyFit/reference/ridge_plots.md)
  : Ridge density plots for time series
- [`plot_missing()`](https://gapouliasis.github.io/anyFit/reference/plot_missing.md)
  : Plot missing value positions

## Distribution functions (d / p / q / r)

Density, distribution, quantile and random-generation functions.

- [`dexp()`](https://gapouliasis.github.io/anyFit/reference/Exponential.md)
  [`pexp()`](https://gapouliasis.github.io/anyFit/reference/Exponential.md)
  [`qexp()`](https://gapouliasis.github.io/anyFit/reference/Exponential.md)
  [`rexp()`](https://gapouliasis.github.io/anyFit/reference/Exponential.md)
  : Exponential Distribution
- [`drayleigh()`](https://gapouliasis.github.io/anyFit/reference/Rayleigh.md)
  [`prayleigh()`](https://gapouliasis.github.io/anyFit/reference/Rayleigh.md)
  [`qrayleigh()`](https://gapouliasis.github.io/anyFit/reference/Rayleigh.md)
  [`rrayleigh()`](https://gapouliasis.github.io/anyFit/reference/Rayleigh.md)
  : Rayleigh Distribution
- [`dgamma()`](https://gapouliasis.github.io/anyFit/reference/Gamma.md)
  [`pgamma()`](https://gapouliasis.github.io/anyFit/reference/Gamma.md)
  [`qgamma()`](https://gapouliasis.github.io/anyFit/reference/Gamma.md)
  [`rgamma()`](https://gapouliasis.github.io/anyFit/reference/Gamma.md)
  : Gamma Distribution
- [`dnorm()`](https://gapouliasis.github.io/anyFit/reference/Normal.md)
  [`pnorm()`](https://gapouliasis.github.io/anyFit/reference/Normal.md)
  [`qnorm()`](https://gapouliasis.github.io/anyFit/reference/Normal.md)
  [`rnorm()`](https://gapouliasis.github.io/anyFit/reference/Normal.md)
  : Normal Distribution
- [`dlognorm()`](https://gapouliasis.github.io/anyFit/reference/Log-Normal.md)
  [`plognorm()`](https://gapouliasis.github.io/anyFit/reference/Log-Normal.md)
  [`qlognorm()`](https://gapouliasis.github.io/anyFit/reference/Log-Normal.md)
  [`rlognorm()`](https://gapouliasis.github.io/anyFit/reference/Log-Normal.md)
  : Log-Normal Distribution
- [`dgenlogi()`](https://gapouliasis.github.io/anyFit/reference/GenLogistic.md)
  [`pgenlogi()`](https://gapouliasis.github.io/anyFit/reference/GenLogistic.md)
  [`qgenlogi()`](https://gapouliasis.github.io/anyFit/reference/GenLogistic.md)
  : Generalized Logistic Distribution
- [`dweibull()`](https://gapouliasis.github.io/anyFit/reference/Weibull3.md)
  [`pweibull()`](https://gapouliasis.github.io/anyFit/reference/Weibull3.md)
  [`qweibull()`](https://gapouliasis.github.io/anyFit/reference/Weibull3.md)
  [`rweibull()`](https://gapouliasis.github.io/anyFit/reference/Weibull3.md)
  : 3-parameter Weibull Distribution
- [`dgumbel()`](https://gapouliasis.github.io/anyFit/reference/Gumbel.md)
  [`pgumbel()`](https://gapouliasis.github.io/anyFit/reference/Gumbel.md)
  [`qgumbel()`](https://gapouliasis.github.io/anyFit/reference/Gumbel.md)
  [`rgumbel()`](https://gapouliasis.github.io/anyFit/reference/Gumbel.md)
  : Gumbel Distribution
- [`dgamma3()`](https://gapouliasis.github.io/anyFit/reference/Gamma3.md)
  [`pgamma3()`](https://gapouliasis.github.io/anyFit/reference/Gamma3.md)
  [`qgamma3()`](https://gapouliasis.github.io/anyFit/reference/Gamma3.md)
  [`rgamma3()`](https://gapouliasis.github.io/anyFit/reference/Gamma3.md)
  : 3-parameter Gamma Distribution
- [`dgev()`](https://gapouliasis.github.io/anyFit/reference/GEV.md)
  [`pgev()`](https://gapouliasis.github.io/anyFit/reference/GEV.md)
  [`qgev()`](https://gapouliasis.github.io/anyFit/reference/GEV.md)
  [`rgev()`](https://gapouliasis.github.io/anyFit/reference/GEV.md) :
  Generalized Extreme Value Distribution
- [`dgpd()`](https://gapouliasis.github.io/anyFit/reference/GenPareto.md)
  [`pgpd()`](https://gapouliasis.github.io/anyFit/reference/GenPareto.md)
  [`qgpd()`](https://gapouliasis.github.io/anyFit/reference/GenPareto.md)
  [`rgpd()`](https://gapouliasis.github.io/anyFit/reference/GenPareto.md)
  : Generalized Pareto Distribution
- [`dgengamma()`](https://gapouliasis.github.io/anyFit/reference/GenGamma.md)
  [`pgengamma()`](https://gapouliasis.github.io/anyFit/reference/GenGamma.md)
  [`qgengamma()`](https://gapouliasis.github.io/anyFit/reference/GenGamma.md)
  [`rgengamma()`](https://gapouliasis.github.io/anyFit/reference/GenGamma.md)
  : Generalized Gamma Distribution
- [`pgengamma_loc()`](https://gapouliasis.github.io/anyFit/reference/GenGamma-Location.md)
  [`dgengamma_loc()`](https://gapouliasis.github.io/anyFit/reference/GenGamma-Location.md)
  [`qgengamma_loc()`](https://gapouliasis.github.io/anyFit/reference/GenGamma-Location.md)
  [`rgengamma_loc()`](https://gapouliasis.github.io/anyFit/reference/GenGamma-Location.md)
  : Generalized Gamma with Location Distribution
- [`dburr()`](https://gapouliasis.github.io/anyFit/reference/BurXII.md)
  [`pburr()`](https://gapouliasis.github.io/anyFit/reference/BurXII.md)
  [`qburr()`](https://gapouliasis.github.io/anyFit/reference/BurXII.md)
  [`rburr()`](https://gapouliasis.github.io/anyFit/reference/BurXII.md)
  : Burr Type XII Distribution
- [`pdagum()`](https://gapouliasis.github.io/anyFit/reference/Dagum.md)
  [`ddagum()`](https://gapouliasis.github.io/anyFit/reference/Dagum.md)
  [`qdagum()`](https://gapouliasis.github.io/anyFit/reference/Dagum.md)
  [`rdagum()`](https://gapouliasis.github.io/anyFit/reference/Dagum.md)
  : Dagum Distribution
- [`pexpweibull()`](https://gapouliasis.github.io/anyFit/reference/ExpWeibull.md)
  [`dexpweibull()`](https://gapouliasis.github.io/anyFit/reference/ExpWeibull.md)
  [`qexpweibull()`](https://gapouliasis.github.io/anyFit/reference/ExpWeibull.md)
  [`rexpweibull()`](https://gapouliasis.github.io/anyFit/reference/ExpWeibull.md)
  : Exponential Weibull Distribution

## Datasets

Data bundled with the package.

- [`world_data`](https://gapouliasis.github.io/anyFit/reference/world_data.md)
  : World country and continent boundary polygons
