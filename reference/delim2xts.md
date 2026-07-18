# delim2xts

Reads a delimited text file and returns an xts object with auto-detected
timestamps. The date column is parsed via
[`parse_date_time`](https://lubridate.tidyverse.org/reference/parse_date_time.html)
using a comprehensive set of order formats (ymd, dmy, mdy, and variants
with hours, minutes, and seconds). When `strict_step = TRUE`, the
function detects the regular time step from the first difference of
parsed dates, constructs a complete index spanning the full date range,
and fills missing positions with `NA` entries. Leap days (29 February)
may be excluded from leap years via `exc_leaps`. A sentinel no-data
value may be specified via `no_value` and replaced with `NA`. The
resulting xts can optionally be saved to a tab-delimited text file.

## Usage

``` r
delim2xts(
  file_path,
  time_zone,
  date_index = 1,
  strict_step = TRUE,
  delim = "\t",
  skip_rows = 0,
  col_names = TRUE,
  no_value = NA,
  exc_leaps = TRUE,
  save_Xts = FALSE,
  filename = NA
)
```

## Arguments

- file_path:

  Character. Path to the delimited data file.

- time_zone:

  Character. Time zone for date parsing (e.g. `"UTC"`).

- date_index:

  Integer. Column position of the date/time variable. Default `1`.

- strict_step:

  Logical. If `TRUE` (default), gaps are filled with `NA` at the
  detected regular interval.

- delim:

  Character. Field delimiter. Default `"\t"`.

- skip_rows:

  Integer. Number of header rows to skip. Default `0`.

- col_names:

  Logical. If `TRUE` (default), the first non-skipped row contains
  column names.

- no_value:

  Value to be treated as missing (`NA`). Default `NA`.

- exc_leaps:

  Logical. If `TRUE` (default), 29 February is excluded from leap years.

- save_Xts:

  Logical. If `TRUE`, the xts is saved as a tab-delimited text file.
  Default `FALSE`.

- filename:

  Character. Output filename (without extension) used when
  `save_Xts = TRUE`. Default `NA`.

## Value

An [`xts`](https://rdrr.io/pkg/xts/man/xts.html) object with parsed
timestamps as the index and the remaining columns as the data matrix.

## Examples

``` r
# Create a temporary delimited file
tf <- tempfile(fileext = ".csv")
writeLines(c("date,value", "2020-01-01,10.5", "2020-01-02,12.3",
             "2020-01-03,11.0", "2020-01-04,9.8"), tf)
x <- delim2xts(tf, time_zone = "UTC", delim = ",", strict_step = TRUE)
#> Time difference of 1 days
head(x)
#>            value
#> 2020-01-01  10.5
#> 2020-01-02  12.3
#> 2020-01-03  11.0
#> 2020-01-04   9.8
```
