# anyFit

Timeseries analysis of hydro-climatic variables is crucial for the study of relevant physical processes as well as the design and forcing of a plethora of physical and statistical models. However, physical processes exhibit a wide range of distinct features that complicate their study. For example:

* Seasonality 
* Over-annual trends - long-term persistence
* Intermittency 
* Highly skewed distributions 
* Different forms of autocorrelation 

Importantly, many of these physical characteristics act on different scales, ranging from very fine (e.g. hourly) to very large (e.g. annual or decadal). Additionally, advances in different technology sectors, with a prominent position of the remote sensing industry, have made publicly accessible, large databases of ground (e.g. NASA MODIS, KNMI) and remote observations,  hindcast and reanalysis models. All of the above imply that tools developed for timeseries analysis, or more commonly known exploratory data analysis (EDA) should:

* Be data agnostic
* Be able to handle multiple temporal scales
* Be easily extensible 
* Be able to handle large sets of data from multiple stations/locations efficiently
* Be able to work with NetCDF files

AnyFit aims to address the above demands in a robust manner. It is based on the eXtensible Time Series (xts) data format, which is a powerful package that provides an extensible time series class, enabling uniform handling of many R time series classes by extending zoo. Moreover, all visualization functions are implemented in ggplot2 package. AnyFit is developed with practitioners and researchers in mind. In this manner it streamlines the EDA with the stochastic modelling of hydrocliamtic variables using the anySim package, thus setting the foundations for a timeseries analysis and modeling ecosystem. Notable examples of anyFit functionality include:

* Read directly from delimited files to xts format
* Identify gaps in data
* Produce summary statistics at various temporal scales (e.g. monthly, annual) with ease
* Perform distribution fitting at various temporal scales (e.g. monthly, annual) with ease
* Provides support for some "exotic" distributions, such as the Dagum and BurXII
* Extent summary statistics and distribution fitting functionality directly to spatial data and produce NetCDF files with results

In detail the distributions supported by anyFit are the following:
* Exponential
* Rayleigh
* Gamma
* Normal
* Log-Normal
* Generalized Logistic 
* 3-parameter Weibull 
* Gumbel
* 3-parameter Gamma (PearsonIII)
* Generalized Extreme Value Distribution (GEV) 
* Generalized Pareto Distribution (GPD) 
* Generalized Gamma 
* Generalized Gamma with location 
* BurXII 
* Dagum 
* Exponentiated Weibull

This package is still in alpha testing. Use at your own risk.

Installation: devtools::install_github("Gapouliasis/anyFit")

Installation with vignette: devtools::install_github("Gapouliasis/anyFit", build_vignette = TRUE)

To view the vignette run browseVignettes("anyFit")
