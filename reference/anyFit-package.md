# anyFit: anyFit: Exploratory Data Analysis on xts timeseries objects and spatial netCDF objects

Package for exploratory data analysis, timeseries analysis and
visualisation of point and spatially distributed environmental
variables.

## Details

anyFit streamlines exploratory data analysis of point and gridded
hydroclimatic time series: loading (NetCDF and delimited files), summary
statistics, temporal aggregation, and distribution fitting by the method
of L-moments, all built on the xts time-series format and ggplot2
visualisation. Gridded data are carried in the
[sxts](https://gapouliasis.github.io/anyFit/reference/sxts.md) (spatial
xts) class, which stores coordinates and projection alongside the
series.

The package ships a real example dataset used throughout the vignettes
and examples: an E-OBS ensemble-mean daily rainfall grid over Europe
(0.25 degrees, 2011-2022), installed at
`system.file("extdata", "rr_ens_mean_0.25deg_reg_2011-2022_v27.0e.nc", package = "anyFit")`.
Load it with
[`nc2xts()`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md),
typically clipped to a region, e.g.
`nc2xts(f, "rr", country = "U.K. of Great Britain and Northern Ireland")`.
A daily station series in
`system.file("extdata", "KNMI_Daily.csv", package = "anyFit")`
illustrates
[`delim2xts()`](https://gapouliasis.github.io/anyFit/reference/delim2xts.md).

## See also

Useful links:

- <https://gapouliasis.github.io/anyFit/>

- <https://github.com/Gapouliasis/anyFit>

- Report bugs at <https://github.com/Gapouliasis/anyFit/issues>

## Author

**Maintainer**: George Pouliasis <gapouliasis@gmail.com>

Authors:

- George Pouliasis <gapouliasis@gmail.com>
