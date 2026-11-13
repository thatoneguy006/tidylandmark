# Print a landmark dataset

Prints a summary of the landmark configuration and per-landmark outcome
counts (at risk, censored, and events by cause), followed by the usual
tibble printout. Automatically updates after
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
or `[` subsetting.

## Usage

``` r
# S3 method for class 'landmark_tbl'
print(x, ...)
```

## Arguments

- x:

  A `landmark_tbl`, as returned by
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md).

- ...:

  Passed on to the tibble print method

## Value

`x`, invisibly
