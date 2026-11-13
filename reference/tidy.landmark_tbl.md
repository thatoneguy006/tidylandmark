# Tidy a landmark dataset for additional summarization

Returns a tibble with the summary of the landmark configuration and
per-landmark outcome counts (at risk, censored, and events by cause).
Automatically updates after
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
or `[` subsetting.

## Usage

``` r
# S3 method for class 'landmark_tbl'
tidy(x, status_labels = FALSE, show_proportion = FALSE, ...)
```

## Arguments

- x:

  A `landmark_tbl`, as returned by
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md).

- status_labels:

  Should status labels be shown in tibble if `status` was provided as
  factor? Default FALSE.

- show_proportion:

  Indicate proportion of full cohort at risk? Default FALSE.

- ...:

  Ignored.

## Value

A tibble
