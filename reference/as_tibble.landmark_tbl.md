# Convert a landmark dataset to a plain tibble

Drops the `landmark_tbl` class and its stored configuration, returning
an ordinary tibble with the same rows and columns.

## Usage

``` r
# S3 method for class 'landmark_tbl'
as_tibble(x, ...)
```

## Arguments

- x:

  A `landmark_tbl`, as returned by
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md).

- ...:

  Ignored.

## Value

A tibble
