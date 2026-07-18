# Monthly Distribution Fitting on Gridded Data

Gridded counterpart of
[`fitlm_monthly`](https://gapouliasis.github.io/anyFit/reference/fitlm_monthly.md)
that fits a list of candidate distributions to every calendar month of a
NetCDF raster or an `sxts` object via the L-moments method. The input is
first normalised to an `sxts` object (accepting `sxts`, `Raster*`, or
NetCDF filename and variable name), and each calendar month present in
the record is extracted as a sub-grid. Per-month fitting is dispatched
to
[`fitlm_nc`](https://gapouliasis.github.io/anyFit/reference/fitlm_nc.md),
which delegates to
[`fitlm_nxts`](https://gapouliasis.github.io/anyFit/reference/fitlm_nxts.md)
for per-cell L-moment fitting and returns parameter rasters, theoretical
and sample L-moment rasters, GoF rasters, and a density GoF comparison
plot. Parallelism is controlled by `parallel_by`: `"cells"` parallelises
the cell-level fits within each month, whereas `"months"` parallelises
across months but is serial across cells. The `shared_memory` flag
determines whether the grid is shared with workers via a file-backed
`big.matrix` or serialised in column chunks (suitable for multi-node
plans). This is reccomended for single machine usage and especially
windows for efficiency and reduced RAM consumption.The returned list is
named by month and contains the full `fitlm_nc` output for each.

## Usage

``` r
fitlm_monthly_nc(
  data = NULL,
  filename = NA,
  varname = NA,
  candidates = "norm",
  ignore_zeros = FALSE,
  zero_threshold = 0.01,
  parallel = FALSE,
  ncores = 2,
  shared_memory = TRUE,
  parallel_by = c("cells", "months"),
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

  Logical, when parallel on the cell axis, share the grid with workers
  via a file-backed `big.matrix` (mmap, single machine) instead of
  serialising column chunks. Set `FALSE` for multi-node `plan(cluster)`
  setups. Default `TRUE`.

- parallel_by:

  Character, the axis to parallelise over when `parallel = TRUE`:
  `"cells"` (default) parallelises the grid-cell fits within each month
  via `fitlm_nxts`; `"months"` parallelises the per-month fits, each run
  serially across cells.

- order:

  Optional named list mapping a candidate name to the vector of L-moment
  orders matched by its optimiser, e.g.
  `list(gengamma = 1:5, expweibull = 1:3)`. Only the numerically-fitted
  distributions accept it; passed through to
  [`fitlm_nc`](https://gapouliasis.github.io/anyFit/reference/fitlm_nc.md).
  Default `NULL`.

- ...:

  Additional arguments passed to
  [`nc2xts`](https://gapouliasis.github.io/anyFit/reference/nc2xts.md)
  when `filename` and `varname` are provided.

## Value

A named list with one element per calendar month present in the data
(named after `month.name`). Each element is the standard
[`fitlm_nc`](https://gapouliasis.github.io/anyFit/reference/fitlm_nc.md)
output: a list with `fit_results` (per-candidate rasters) and
`gof_plots`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulated 3-cell grid over two years
n <- 730
dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = n)
vals <- cbind(cell1 = rgamma(n, shape = 0.8, scale = 3),
              cell2 = rgamma(n, shape = 1.2, scale = 2),
              cell3 = rgamma(n, shape = 0.6, scale = 4))
coords <- data.frame(lon = c(10, 11, 12), lat = c(45, 46, 47))
grid <- sxts(vals, order.by = dates, coords = coords,
             projection = "+proj=longlat +datum=WGS84")

monthly_fits <- fitlm_monthly_nc(grid, candidates = c("exp", "gamma3"),
                                  ignore_zeros = TRUE)
names(monthly_fits)
} # }
```
