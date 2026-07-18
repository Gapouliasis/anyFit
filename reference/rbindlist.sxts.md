# Fast row-bind of a list of sxts objects

Drop-in replacement for `do.call(rbind.xts, list)` that binds a list of
sxts objects by rows in a single pass, avoiding the O(n^2) complexity of
pairwise merging. Coordinates and projection are taken from the first
element.

## Usage

``` r
rbindlist.sxts(sxts_list)
```

## Arguments

- sxts_list:

  A list of `sxts` objects sharing the same columns (spatial locations).
  Coordinates and projection are taken from the first element.
