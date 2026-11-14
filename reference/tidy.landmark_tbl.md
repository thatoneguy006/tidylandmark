# Tidy a landmark dataset

Returns the per-landmark outcome counts that
[`print.landmark_tbl()`](https://thatoneguy006.github.io/tidylandmark/reference/print.landmark_tbl.md)
shows, as a tibble with one row per landmark: the number at risk, the
number censored within the horizon, and the number of events of each
type. Counts are recomputed from the current rows, so they follow
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
and `[` subsetting.

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

  If `TRUE` and `status` was supplied to
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md)
  as a factor, name the event columns after the factor levels (for
  example `n_relapse`, `n_death`) instead of the integer codes
  (`n_event_1`, `n_event_2`). Levels are made syntactic with
  [`make.names()`](https://rdrr.io/r/base/make.names.html). When
  `status` was numeric there are no labels to use, so the columns keep
  their code names and a warning is given. Default `FALSE`.

- show_proportion:

  If `TRUE`, add `prop_at_risk`: the number at risk divided by the
  number of subjects in the source data. Default `FALSE`.

- ...:

  Must be empty. Present for compatibility with the generic.

## Value

A tibble with columns `landmark`, `prop_at_risk` (only when
`show_proportion = TRUE`), `n_at_risk`, `n_censored`, and one
event-count column per event type.

## Examples

``` r
subjects <- data.frame(
  patient_id    = 1:10,
  followup_days = c(20, 30, 48, 65, 85, 110, 145, 190, 250, 270),
  event         = c(1, 2, 1, 2, 0, 1, 2, 1, 0, 1)
)
lm_dat <- landmark_data(
  subjects, patient_id, followup_days, event,
  landmarks = c(30, 60, 90), horizon = 180
)
tidy(lm_dat)
#> # A tibble: 3 × 5
#>   landmark n_at_risk n_censored n_event_1 n_event_2
#>      <dbl>     <int>      <int>     <int>     <int>
#> 1       30         8          3         3         2
#> 2       60         7          3         2         2
#> 3       90         5          1         3         1
tidy(lm_dat, show_proportion = TRUE)
#> # A tibble: 3 × 6
#>   landmark prop_at_risk n_at_risk n_censored n_event_1 n_event_2
#>      <dbl>        <dbl>     <int>      <int>     <int>     <int>
#> 1       30          0.8         8          3         3         2
#> 2       60          0.7         7          3         2         2
#> 3       90          0.5         5          1         3         1

# With a factor status, event columns can be named after its levels
subjects$event <- factor(
  subjects$event, levels = 0:2, labels = c("censored", "relapse", "death")
)
lm_fac <- landmark_data(
  subjects, patient_id, followup_days, event,
  landmarks = c(30, 60, 90), horizon = 180
)
tidy(lm_fac, status_labels = TRUE)
#> # A tibble: 3 × 5
#>   landmark n_at_risk n_censored n_relapse n_death
#>      <dbl>     <int>      <int>     <int>   <int>
#> 1       30         8          3         3       2
#> 2       60         7          3         2       2
#> 3       90         5          1         3       1
```
