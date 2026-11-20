# Attach pre-landmark covariates to a landmark dataset

`add_landmark_covariates()` adds one column per selected longitudinal
variable to a `landmark_tbl`. For every subject-landmark row, the value
is derived only from that subject's measurements taken at or before the
landmark, so no information from after the landmark leaks into the
covariates. Each variable is derived in one of two ways:

- **`locf`**: last observation carried forward, optionally limited to a
  `lookback` window.

- **`lme`**: the subject-specific value at the landmark predicted from a
  linear mixed-effects model fitted separately for every variable and
  landmark.

Rows are never added or removed, and `.time`/`.status` are untouched.
Subjects with no usable measurement keep their row with a missing value.

## Usage

``` r
add_landmark_covariates(
  landmark_data,
  longitudinal_data,
  id,
  time,
  locf = NULL,
  lme = NULL,
  include_at_landmark = TRUE,
  lookback = NULL,
  lme_args = list(),
  names = NULL,
  progress = getOption("tidylandmark.progress", interactive())
)
```

## Arguments

- landmark_data:

  A `landmark_tbl`, as returned by
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md)
  (or by a previous call to `add_landmark_covariates()`).

- longitudinal_data:

  A data frame with **one row per measurement occasion** (multiple rows
  per subject). Rows whose `id` does not appear in `landmark_data` are
  ignored, and their number is reported in a message.

- id:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  The subject identifier column in `longitudinal_data`. It must be
  non-missing and have the same type as the `id` used in
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md);
  types are never coerced implicitly.

- time:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Numeric, finite, non-missing measurement time. It must be on the
  **same time scale and origin** as the `time` given to
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md);
  this cannot be checked automatically. Measurements before time zero
  (e.g. screening visits) are allowed.

- locf, lme:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Columns of `longitudinal_data` to derive by last observation carried
  forward and by linear mixed-effects model, e.g.
  `locf = c(treatment, stage)` and `lme = c(bili, albumin)`. Either may
  be `NULL` (the default), but not both. Rename inside the selection to
  name a derived column, e.g. `lme = c(bili_lme = bili)`. `lme`
  variables must be continuous numeric.

- include_at_landmark:

  If `TRUE` (default), a measurement taken exactly at the landmark (\\t
  = s\\) may be used. If `FALSE`, only measurements strictly before the
  landmark (\\t \< s\\) are used. Applies to both methods.

- lookback:

  `locf` only. `NULL` (default) carries values forward without a limit.
  A single non-negative number only uses measurements within that
  window: \\s - t \le\\ `lookback`. Older values give `NA`.

- lme_args:

  `lme` only. A named list of model options:

  - `random_slope`: `FALSE` (default) fits a random intercept; `TRUE`
    adds a random slope for time.

  - `adjust`: `NULL` (default) or a one-sided formula of fixed-effect
    adjustment terms, e.g. `~ age + sex`. Variables come from
    `landmark_data`, so they are baseline values (or covariates added by
    an earlier call). Use `.t` for the time since the landmark, e.g.
    `~ sex * .t` for sex-specific slopes.

  - `REML`: passed to
    [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html). Default
    `TRUE`.

  - `control`: an
    [`lme4::lmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html)
    object. Default
    [`lme4::lmerControl()`](https://rdrr.io/pkg/lme4/man/lmerControl.html).

  - `min_subjects`: the fewest subjects with measurements needed to fit
    a model (default `5`). With `random_slope = TRUE`, this many
    subjects must have at least two measurements.

- names:

  `NULL` (default) or a template for derived column names, in which
  `"{.col}"` is replaced by the source variable and `"{.method}"` by
  `"locf"` or `"lme"`, e.g. `"{.col}_{.method}"`. Columns renamed in the
  selection keep their new name. A derived name that already exists in
  `landmark_data` is an error; existing columns are never overwritten.

- progress:

  Show a progress bar while fitting mixed models? Defaults to the
  `tidylandmark.progress` option, or
  [`interactive()`](https://rdrr.io/r/base/interactive.html) if it is
  unset.

## Value

`landmark_data` with one new column per selected variable (`locf`
columns first, then `lme`), appended after the existing columns. The
class and configuration are kept, and the configuration records how each
derived column was made: method, source variable, options, the LOCF
audit trail, and for LME the fitted models and their diagnostics. The
printed header lists the derived columns and their method.

## Details

### LOCF

For landmark \\s\\, subject \\i\\ and variable \\v\\, the value is the
last *non-missing* measurement of \\v\\ at time \\t \le s\\ (or \\t \<
s\\ when `include_at_landmark = FALSE`), provided it is within
`lookback`. Each variable is handled separately, so a visit where one
variable is missing still contributes the others. Selection depends only
on `id` and `time`, never on row order. Any column type works, and the
type (including factor levels) is preserved. The time since each
carried-forward measurement is stored in the metadata for auditing.

If a subject has two *different* non-missing values of an `locf`
variable at the same time, the function errors and lists the affected
subjects and times rather than choosing one arbitrarily. Exact
duplicates are collapsed.

### LME

For every `lme` variable \\v\\ and landmark \\s\\, an independent
Gaussian model is fitted with
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) to the
non-missing measurements, taken at or before \\s\\, of the subjects in
the risk set at \\s\\:

\$\$Y\_{ij} = \beta_0 + \beta_1 (t\_{ij} - s) + \gamma^\top x_i +
b\_{0i} \[+ b\_{1i}(t\_{ij} - s)\] + \epsilon\_{ij}\$\$

where \\x_i\\ are the optional `adjust` terms. Time is centred at the
landmark and rescaled internally for numerical stability, which does not
change the predictions. The derived value is the conditional prediction
at \\t = s\\ including the subject's random effects
(`predict(fit, re.form = NULL)`), i.e. \\\hat\beta_0 + \hat\gamma^\top
x_i + \hat b\_{0i}\\. A subject with a single measurement is kept; their
prediction is shrunk towards the population trajectory. A subject with
no measurement before the landmark, or with missing `adjust` values,
gets `NA` rather than a population-only prediction. Different values
recorded at the same time are treated as replicate measurements.

Each fit gets a status: `"ok"`; `"warning"` when
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) warned (e.g.
about convergence) or the fit is singular
([`lme4::isSingular()`](https://rdrr.io/pkg/lme4/man/isSingular.html));
or `"failed"` when there were too few subjects or
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) gave an error.
Failed fits give `NA` for that variable at that landmark without
affecting other fits, and a warning lists every fit that was not `"ok"`.
There is no automatic fallback from a random slope to a random
intercept.

Standard errors in a second-stage survival model treat LME-derived
values as known, so they understate uncertainty.

## See also

[`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md)
to build the risk sets first.

## Examples

``` r
subjects <- data.frame(
  patient_id    = 1:10,
  followup_days = c(20, 30, 48, 65, 85, 110, 145, 190, 250, 270),
  event         = c(1, 2, 1, 2, 0, 1, 2, 1, 0, 1)
)
visits <- data.frame(
  patient_id = c(6L, 6L, 6L, 7L, 7L, 7L, 10L, 10L),
  visit_day  = c(0, 45, 90, 0, 50, 120, 60, 140),
  treatment  = c(0, 1, 1, 0, 0, 1, 1, 1),
  biomarker  = c(25, 20, 17, 23, NA, 18, 24, 20)
)

lm_dat <- subjects |>
  landmark_data(
    patient_id, followup_days, event,
    landmarks = c(30, 60, 90), horizon = 180
  ) |>
  add_landmark_covariates(
    visits, id = patient_id, time = visit_day,
    locf = c(treatment, biomarker)
  )
lm_dat
#> <landmark_tbl> 3 landmark(s): 30, 60, 90
#> Horizon: 180 (relative to each landmark)
#> Subjects: 8 at risk in current rows (10 in source data)
#> Source columns: id = patient_id, time = followup_days, status = event
#> Derived covariates: treatment (LOCF), biomarker (LOCF)
#>  landmark n_at_risk n_censored n_event_1 n_event_2
#>        30         8          3         3         2
#>        60         7          3         2         2
#>        90         5          1         3         1
#> 
#> # A tibble: 20 × 7
#>      .id .landmark .time .status .horizon treatment biomarker
#>    <int>     <dbl> <dbl>   <int>    <dbl>     <dbl>     <dbl>
#>  1     3        30    18       1      180        NA        NA
#>  2     4        30    35       2      180        NA        NA
#>  3     5        30    55       0      180        NA        NA
#>  4     6        30    80       1      180         0        25
#>  5     7        30   115       2      180         0        23
#>  6     8        30   160       1      180        NA        NA
#>  7     9        30   180       0      180        NA        NA
#>  8    10        30   180       0      180        NA        NA
#>  9     4        60     5       2      180        NA        NA
#> 10     5        60    25       0      180        NA        NA
#> 11     6        60    50       1      180         1        20
#> 12     7        60    85       2      180         0        23
#> 13     8        60   130       1      180        NA        NA
#> 14     9        60   180       0      180        NA        NA
#> 15    10        60   180       0      180         1        24
#> 16     6        90    20       1      180         1        17
#> 17     7        90    55       2      180         0        23
#> 18     8        90   100       1      180        NA        NA
#> 19     9        90   160       0      180        NA        NA
#> 20    10        90   180       1      180         1        24

# Only carry values forward for up to 30 days, under new names
landmark_data(subjects, patient_id, followup_days, event, landmarks = 60) |>
  add_landmark_covariates(
    visits, patient_id, visit_day,
    locf = biomarker, lookback = 30, names = "{.col}_recent"
  )
#> <landmark_tbl> 1 landmark(s): 60
#> Horizon: none (follow-up to observed time)
#> Subjects: 7 at risk in current rows (10 in source data)
#> Source columns: id = patient_id, time = followup_days, status = event
#> Derived covariates: biomarker_recent (LOCF of biomarker, lookback 30)
#>  landmark n_at_risk n_censored n_event_1 n_event_2
#>        60         7          2         3         2
#> 
#> # A tibble: 7 × 6
#>     .id .landmark .time .status .horizon biomarker_recent
#>   <int>     <dbl> <dbl>   <int>    <dbl>            <dbl>
#> 1     4        60     5       2       NA               NA
#> 2     5        60    25       0       NA               NA
#> 3     6        60    50       1       NA               20
#> 4     7        60    85       2       NA               NA
#> 5     8        60   130       1       NA               NA
#> 6     9        60   190       0       NA               NA
#> 7    10        60   210       1       NA               24

# Mixed models with the PBC data: log bilirubin by LME, adjusted
# for age, alongside LOCF for albumin
visits_pbc <- transform(pbcseq_visits, log_bili = log(bili))
pbc_lm <- pbcseq_baseline |>
  landmark_data(id, futime, status, landmarks = c(365, 730), horizon = 1825) |>
  add_landmark_covariates(
    visits_pbc, id, day,
    locf = albumin,
    lme = log_bili,
    lme_args = list(adjust = ~ age)
  )
#> ℹ Ignored 32 rows of `longitudinal_data` whose `id` is not in `landmark_data`.
pbc_lm
#> <landmark_tbl> 2 landmark(s): 365, 730
#> Horizon: 1825 (relative to each landmark)
#> Subjects: 290 at risk in current rows (312 in source data)
#> Source columns: id = id, time = futime, status = status
#> Derived covariates: albumin (LOCF), log_bili (LME, adjusted ~age)
#>  landmark n_at_risk n_censored n_event_1 n_event_2
#>       365       290        194        20        76
#>       730       278        176        26        76
#> 
#> # A tibble: 568 × 10
#>      .id .landmark .time .status .horizon   trt   age sex   albumin log_bili
#>    <int>     <dbl> <dbl>   <int>    <dbl> <int> <dbl> <fct>   <dbl>    <dbl>
#>  1     1       365    35       2     1825     1  58.8 f        2.94   2.68  
#>  2     2       365  1825       0     1825     1  56.4 f        3.55  -0.0310
#>  3     3       365   647       2     1825     1  70.1 m        3.57   0.272 
#>  4     4       365  1560       2     1825     1  54.7 f        2.88   0.507 
#>  5     5       365  1140       1     1825     0  38.1 f        3.28   0.889 
#>  6     6       365  1825       0     1825     0  66.3 f        3.98  -0.166 
#>  7     7       365  1825       0     1825     0  55.5 f        4.09   0.0356
#>  8     8       365  1825       0     1825     0  53.1 f        3.86  -0.867 
#>  9     9       365  1825       0     1825     1  42.5 f        3.1    1.46  
#> 10    11       365  1825       0     1825     0  53.7 f        3.71   0.256 
#> # ℹ 558 more rows
```
