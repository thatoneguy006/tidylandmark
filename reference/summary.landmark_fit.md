# Summarise a landmark model fit

A shortcut for `tidy(object, exponentiate = TRUE, conf.int = TRUE)`:
hazard ratios (or subdistribution hazard ratios) with confidence
intervals for every landmark and term.

## Usage

``` r
# S3 method for class 'landmark_fit'
summary(object, conf.level = 0.95, ...)
```

## Arguments

- object:

  A `landmark_fit`, as returned by
  [`fit_landmark()`](https://thatoneguy006.github.io/tidylandmark/reference/fit_landmark.md).

- conf.level:

  Confidence level. Default `0.95`.

- ...:

  Must be empty.

## Value

A tibble; see
[`tidy.landmark_fit()`](https://thatoneguy006.github.io/tidylandmark/reference/tidy.landmark_fit.md).
