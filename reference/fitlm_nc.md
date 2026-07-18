# Distribution Fitting on Gridded Data

Fits a list of candidate distributions to every grid cell of a NetCDF
raster or an `sxts` object via the L-moments method. The input is
normalised to an `sxts` object, accepting `sxts`, `Raster*`, or NetCDF
filename and variable name. The results are returned as per-candidate
`RasterBrick` objects for the parameters, theoretical L-moments, data
L-moments, and GoF metrics. A density comparison plot of the six GoF
metrics across all grid cells is assembled with `ggplot2`, faceted by
metric and coloured by candidate distribution. When `parallel = TRUE`,
the work is distributed across available cores via the future framework,
with either file-backed shared memory (bigmemory memory-mapped matrices
for single-machine parallelism) or serialised column chunks (for
multi-node clusters). File-backed shared memory can be enabled by
`shared_memory = TRUE`. This is reccomended for single machine usage and
especially windows for efficiency and reduced RAM consumption.

## Usage

``` r
fitlm_nc(
  data = NULL,
  filename = NA,
  varname = NA,
  candidates = "norm",
  ignore_zeros = FALSE,
  zero_threshold = 0.01,
  parallel = FALSE,
  ncores = 2,
  shared_memory = TRUE,
  order = NULL,
  ...
)
```

## Arguments

- data:

  An `sxts` object, or a `Raster*` object. Leave `NULL` when supplying
  `filename` and `varname`.

- filename:

  A NetCDF file name to import if `data` is not provided.

- varname:

  The name of the variable to extract from `filename`.

- candidates:

  A character vector of distribution names to fit.

- ignore_zeros:

  A logical value, if `TRUE` zeros will be ignored. Default is `FALSE`.

- zero_threshold:

  The threshold below which values are considered zero. Default is 0.01.

- parallel:

  Logical, whether to use parallel processing.

- ncores:

  Number of cores to use for parallel computations. Default is 2.

- shared_memory:

  Logical, when parallel, share the grid with workers via a file-backed
  `big.matrix` (mmap, single machine) instead of serialising column
  chunks. Set `FALSE` for multi-node `plan(cluster)` setups. Default
  `TRUE`.

- order:

  Optional named list mapping a candidate name to the vector of L-moment
  orders matched by its optimiser, e.g.
  `list(gengamma = 1:5, expweibull = 1:3)`. Only the numerically-fitted
  distributions accept it; passed through to
  [`fitlm_nxts`](https://gapouliasis.github.io/anyFit/reference/fitlm_nxts.md).
  Default `NULL`.

- ...:

  Additional arguments passed to
  [`nc2xts`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md)
  when `filename` and `varname` are provided.

## Value

A list with components `fit_results` (a named list, one element per
candidate, each containing `raster_params`, `raster_TheorLMom`,
`raster_DataLMom`, and `raster_GoF` `RasterBrick` objects) and
`gof_plots` (a `ggplot` density comparison of the six GoF metrics across
all grid cells, faceted by metric and coloured by candidate).

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulated 3-cell grid
n <- 365
dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = n)
vals <- cbind(cell1 = rgamma(n, shape = 0.8, scale = 3),
              cell2 = rgamma(n, shape = 1.2, scale = 2),
              cell3 = rgamma(n, shape = 0.6, scale = 4))
coords <- data.frame(lon = c(10, 11, 12), lat = c(45, 46, 47))
grid <- sxts(vals, order.by = dates, coords = coords,
             projection = "+proj=longlat +datum=WGS84")

fits <- fitlm_nc(grid, candidates = c("exp", "gamma3"), ignore_zeros = TRUE)
names(fits$fit_results)
fits$gof_plots
} # }
```
