# Construct landmark risk sets

`landmark_data()` builds the outcome side of a landmark analysis. For
each landmark time \\s\\ it keeps every subject who is still under
observation and event-free at \\s\\. It moves the time origin to \\s\\
and, optionally, administratively censors follow-up at a fixed
post-landmark `horizon`. The result stacks all landmarks into one long
table with one row per eligible subject-landmark pair.

Longitudinal covariates are not handled here. Attach them afterwards
with `add_landmark_covariates()`.

## Usage

``` r
landmark_data(data, id, time, status, landmarks, horizon = NULL)
```

## Arguments

- data:

  A data frame with **one row per subject**, holding follow-up time,
  event status and any baseline predictors.

- id:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  The subject identifier column. It must be unique and non-missing. Its
  type is preserved in `.id`.

- time:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric, non-negative, finite time from the common origin (time zero)
  to the first event or censoring.

- status:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Event indicator. It is either:

  - **numeric**: `0` for censoring and positive integers `1, 2, ...` for
    event types. Codes other than the target cause are treated as
    competing events by `fit_landmark()`.

  - **factor**: the first level is censoring and later levels are event
    types. They are recoded to `0, 1, 2, ...` in level order. The
    mapping is stored in the object's configuration.

- landmarks:

  A non-empty numeric vector of unique, non-negative landmark times on
  the original time scale. Landmarks are processed in ascending order
  regardless of input order.

- horizon:

  `NULL` (default) for follow-up to each subject's observed time, or a
  single positive number giving the follow-up window *relative to each
  landmark*. For example, `horizon = 180` at landmark `90` ends
  follow-up at original time `270`.

## Value

A `landmark_tbl`, which is a tibble subclass. Its leading columns are:

- `.id`: subject identifier, with the original type.

- `.landmark`: landmark time on the original scale.

- `.time`: follow-up since the landmark. It is always strictly positive.

- `.status`: integer, `0` for censored or `1, 2, ...` for the event
  type.

- `.horizon`: the configured horizon, or `NA` when follow-up is
  unbounded.

The remaining baseline columns follow. Rows are ordered by landmark,
then by the original row order of `data`. Use
[`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
to drop the class and metadata.

## Details

**Eligibility.** A subject is in the risk set at landmark \\s\\ exactly
when \\T_i \> s\\. A subject whose event or censoring is recorded *at*
\\s\\ is not eligible.

**Time and status.** For eligible subjects, `.time` is \\T_i - s\\ when
`horizon` is `NULL`, and \\\min(T_i - s, w)\\ otherwise. `.status` keeps
the original code when \\T_i - s \le w\\, so an event exactly at the
horizon is retained. Otherwise `.status` becomes `0`.

**Columns.** The original `id`, `time` and `status` columns are replaced
by the reserved columns below. All other columns are carried over
unchanged, and subjects are never dropped for missing baseline
predictors. If `data` already contains a reserved column name, the
function throws an error.

A landmark with no eligible subjects triggers a warning. It contributes
zero rows but still appears in the configuration and in the printed
counts.

## Examples

``` r
subjects <- data.frame(
  patient_id    = 1:10,
  followup_days = c(20, 30, 48, 65, 85, 110, 145, 190, 250, 270),
  event         = c(1, 2, 1, 2, 0, 1, 2, 1, 0, 1),
  age           = c(12, 15, 11, 17, 9, 14, 13, 16, 10, 15)
)

lm_dat <- landmark_data(
  subjects,
  id = patient_id,
  time = followup_days,
  status = event,
  landmarks = c(30, 60, 90),
  horizon = 180
)
lm_dat
#> <landmark_tbl> 3 landmark(s): 30, 60, 90
#> Horizon: 180 (relative to each landmark)
#> Subjects: 8 at risk in current rows (10 in source data)
#> Source columns: id = patient_id, time = followup_days, status = event
#>  landmark n_at_risk n_censored n_event_1 n_event_2
#>        30         8          3         3         2
#>        60         7          3         2         2
#>        90         5          1         3         1
#> 
#> # A tibble: 20 × 6
#>      .id .landmark .time .status .horizon   age
#>    <int>     <dbl> <dbl>   <int>    <dbl> <dbl>
#>  1     3        30    18       1      180    11
#>  2     4        30    35       2      180    17
#>  3     5        30    55       0      180     9
#>  4     6        30    80       1      180    14
#>  5     7        30   115       2      180    13
#>  6     8        30   160       1      180    16
#>  7     9        30   180       0      180    10
#>  8    10        30   180       0      180    15
#>  9     4        60     5       2      180    17
#> 10     5        60    25       0      180     9
#> 11     6        60    50       1      180    14
#> 12     7        60    85       2      180    13
#> 13     8        60   130       1      180    16
#> 14     9        60   180       0      180    10
#> 15    10        60   180       0      180    15
#> 16     6        90    20       1      180    14
#> 17     7        90    55       2      180    13
#> 18     8        90   100       1      180    16
#> 19     9        90   160       0      180    10
#> 20    10        90   180       1      180    15

# Without a horizon, follow-up runs to each subject's observed time
landmark_data(subjects, patient_id, followup_days, event, landmarks = 60)
#> <landmark_tbl> 1 landmark(s): 60
#> Horizon: none (follow-up to observed time)
#> Subjects: 7 at risk in current rows (10 in source data)
#> Source columns: id = patient_id, time = followup_days, status = event
#>  landmark n_at_risk n_censored n_event_1 n_event_2
#>        60         7          2         3         2
#> 
#> # A tibble: 7 × 6
#>     .id .landmark .time .status .horizon   age
#>   <int>     <dbl> <dbl>   <int>    <dbl> <dbl>
#> 1     4        60     5       2       NA    17
#> 2     5        60    25       0       NA     9
#> 3     6        60    50       1       NA    14
#> 4     7        60    85       2       NA    13
#> 5     8        60   130       1       NA    16
#> 6     9        60   190       0       NA    10
#> 7    10        60   210       1       NA    15
```
