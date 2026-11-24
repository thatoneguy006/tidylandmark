# Fit a survival model at every landmark

`fit_landmark()` fits one independent regression model per landmark of a
`landmark_tbl`: a Cox proportional hazards model
([`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)) or
a Fine-Gray subdistribution hazards model
([`survival::finegray()`](https://rdrr.io/pkg/survival/man/finegray.html)
followed by a weighted
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)). The
outcome is built internally from `.time` and `.status`, so `formula`
only lists the covariates.

## Usage

``` r
fit_landmark(
  data,
  formula,
  model = c("cox", "finegray"),
  cause = 1,
  progress = getOption("tidylandmark.progress", interactive())
)
```

## Arguments

- data:

  A `landmark_tbl`, as returned by
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md)
  or
  [`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md).

- formula:

  A one-sided formula of covariates, e.g. `~ age + sex + bili`. Any
  right-hand side that
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
  accepts can be used, including factors, interactions, transformations
  (`log(bili)`), splines (`splines::ns(age, 3)`) and `strata()`. Leave
  out `Surv()`. The formula can't use `.`, reserved columns (`.time`,
  `.status`, ...), `tt()` or `cluster()`.

- model:

  `"cox"` (default) or `"finegray"`.

- cause:

  The event of interest: one of the event codes of the `status` given to
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md)
  (e.g. `1`), or, if `status` was a factor, one of its event levels
  (e.g. `"relapse"`). Default `1`.

- progress:

  Show a progress bar while fitting? Defaults to the
  `tidylandmark.progress` option, or
  [`interactive()`](https://rdrr.io/r/base/interactive.html) if it is
  unset.

## Value

A `landmark_fit` object: a list with

- `data`: the `landmark_tbl` that was fitted.

- `models`: one fitted
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
  object per landmark, named by landmark (`NULL` where no model was
  fitted). For Fine-Gray models this is the weighted fit to the expanded
  data.

- `model_info`: a tibble with one row per landmark: `.landmark`,
  `.model`, `.cause`, `n_eligible` (rows at risk), `n_model` (complete
  cases used), `n_event` and `n_competing` (target and competing events
  in the model sample), `.fit_status` (`"ok"`, `"warning"`,
  `"unestimable"` or `"failed"`) and `.message` (`NA` when `"ok"`).

- `coefficients`: a tibble with one row per landmark and model term:
  `.landmark`, `term`, `estimate` and `std.error` (log scale) and
  `origin` (`"baseline"`, `"locf"`, `"lme"`, or e.g. `"baseline:lme"`
  for an interaction).

- `warnings`: a tibble of the warnings captured from each fit.

- `config`: the formula, model, cause and related settings.

Use
[tidy()](https://thatoneguy006.github.io/tidylandmark/reference/tidy.landmark_fit.md)
for a results table and [`print()`](https://rdrr.io/r/base/print.html)
or [`summary()`](https://rdrr.io/r/base/summary.html) for an overview.

## Details

### Models

With `model = "cox"`, each landmark gets
`coxph(Surv(.time, .status == cause) ~ <formula>)` with Efron's method
for ties. If the status coding has a single event type, this is an
ordinary Cox model. If it has other event types, they are treated as
censored at their event time, so the model estimates **cause-specific**
hazard ratios; the `.model` column then reads `"cause-specific cox"`.

With `model = "finegray"`,
[`survival::finegray()`](https://rdrr.io/pkg/survival/man/finegray.html)
expands each landmark's data with time-dependent censoring weights for
`cause`, and a weighted
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) is
fitted to the expanded data with a robust variance clustered on subject,
as the
[`survival::finegray()`](https://rdrr.io/pkg/survival/man/finegray.html)
documentation recommends. Coefficients are log subdistribution hazard
ratios. The censoring distribution is estimated on each landmark's
complete-case model sample and is not stratified; `strata()` terms in
`formula` stratify the weighted Cox model only.

### Missing covariates and failed fits

Each landmark's model uses the rows with no missing values in the
variables of `formula`; `n_eligible` and `n_model` in the results show
how many rows that leaves. A landmark is not fitted (status
`"unestimable"`) when it has no subjects at risk, no complete cases, no
target events in the model sample, or a covariate with a single observed
value. A fit that warns (for example about convergence, or a coefficient
that may be infinite) gets status `"warning"`, and a fit that errors
gets `"failed"`. Neither affects the other landmarks, and one warning at
the end lists every landmark that did not fit cleanly. A Fine-Gray model
at a landmark with no competing events is still fitted (it then equals
the Cox model) with status `"warning"`. If a covariate made by an `lme`
model in
[`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md)
did not fit cleanly at a landmark, that landmark's status is at least
`"warning"` and its message names the covariate.

### Interpretation

Each model is conditional on being event-free and under follow-up at its
landmark, so different landmarks describe different populations. Fits at
different landmarks share subjects and are not independent, so
differences between landmarks are not tested. Standard errors treat
covariates derived by
[`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md)
as known values; for `lme` covariates they are therefore too small,
because the mixed models were estimated from the same subjects.

## See also

[`tidy.landmark_fit()`](https://thatoneguy006.github.io/tidylandmark/reference/tidy.landmark_fit.md),
[`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md),
[`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md)

## Examples

``` r
# PBC: status 0 = censored, 1 = transplant, 2 = death
pbc_lm <- pbcseq_baseline |>
  landmark_data(id, futime, status, landmarks = c(365, 730), horizon = 1825) |>
  add_landmark_covariates(pbcseq_visits, id, day, locf = c(bili, albumin))
#> ℹ Ignored 32 rows of `longitudinal_data` whose `id` is not in `landmark_data`.

# Cause-specific Cox model for death (transplant is a competing event)
fit <- fit_landmark(pbc_lm, ~ age + log(bili) + albumin, cause = 2)
fit
#> <landmark_fit> Cause-specific Cox model for cause 2 at 2 landmark(s)
#> Formula: ~age + log(bili) + albumin
#> Horizon: 1825 (relative to each landmark)
#>  landmark n_eligible n_model n_event n_competing fit_status
#>       365        290     290      76          20         ok
#>       730        278     278      76          26         ok
#> 
#> Use tidy() for coefficients and summary() for hazard ratios.
tidy(fit, exponentiate = TRUE, conf.int = TRUE)
#> # A tibble: 6 × 17
#>   .landmark term      estimate std.error statistic  p.value conf.low conf.high
#>       <dbl> <chr>        <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1       365 age          1.05     0.0108      4.49 6.97e- 6    1.03      1.07 
#> 2       365 log(bili)    3.39     0.127       9.60 8.36e-22    2.64      4.35 
#> 3       365 albumin      0.572    0.280      -2.00 4.57e- 2    0.330     0.989
#> 4       730 age          1.05     0.0108      4.81 1.55e- 6    1.03      1.08 
#> 5       730 log(bili)    3.32     0.129       9.30 1.45e-20    2.58      4.27 
#> 6       730 albumin      0.610    0.271      -1.82 6.82e- 2    0.358     1.04 
#> # ℹ 9 more variables: origin <chr>, n_eligible <int>, n_model <int>,
#> #   n_event <int>, n_competing <int>, .model <chr>, .cause <int>,
#> #   .fit_status <chr>, .message <chr>

# Fine-Gray model for death
fit_fg <- fit_landmark(
  pbc_lm, ~ age + log(bili) + albumin, model = "finegray", cause = 2
)
summary(fit_fg)
#> # A tibble: 6 × 17
#>   .landmark term      estimate std.error statistic  p.value conf.low conf.high
#>       <dbl> <chr>        <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1       365 age          1.06     0.0124     4.35  1.33e- 5    1.03       1.08
#> 2       365 log(bili)    3.26     0.113     10.5   8.93e-26    2.62       4.07
#> 3       365 albumin      0.624    0.300     -1.57  1.17e- 1    0.347      1.12
#> 4       730 age          1.06     0.0132     4.71  2.52e- 6    1.04       1.09
#> 5       730 log(bili)    3.03     0.119      9.30  1.45e-20    2.40       3.82
#> 6       730 albumin      0.740    0.308     -0.976 3.29e- 1    0.404      1.35
#> # ℹ 9 more variables: origin <chr>, n_eligible <int>, n_model <int>,
#> #   n_event <int>, n_competing <int>, .model <chr>, .cause <int>,
#> #   .fit_status <chr>, .message <chr>

# The native coxph() fit for one landmark
fit$models[["365"]]
#> Call:
#> survival::coxph(formula = Surv(.time, .status == 2L) ~ age + 
#>     log(bili) + albumin, data = d, model = TRUE)
#> 
#>               coef exp(coef) se(coef)      z        p
#> age        0.04848   1.04968  0.01079  4.495 6.97e-06
#> log(bili)  1.22076   3.38977  0.12722  9.595  < 2e-16
#> albumin   -0.55909   0.57173  0.27987 -1.998   0.0457
#> 
#> Likelihood ratio test=138.2  on 3 df, p=< 2.2e-16
#> n= 290, number of events= 76 
```
