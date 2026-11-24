# Tidy a landmark model fit

Returns one row per landmark and model term, combining the coefficients
of every landmark model with that landmark's counts and fit status.
Landmarks that were not fitted (`"unestimable"` or `"failed"`) keep one
row, with `term` and the estimates `NA`, so they don't silently
disappear from a results table.

## Usage

``` r
# S3 method for class 'landmark_fit'
tidy(x, exponentiate = FALSE, conf.int = FALSE, conf.level = 0.95, ...)
```

## Arguments

- x:

  A `landmark_fit`, as returned by
  [`fit_landmark()`](https://thatoneguy006.github.io/tidylandmark/reference/fit_landmark.md).

- exponentiate:

  If `TRUE`, report hazard ratios (subdistribution hazard ratios for
  Fine-Gray models) instead of log hazard ratios. `std.error` stays on
  the log scale. Default `FALSE`.

- conf.int:

  If `TRUE`, add `conf.low` and `conf.high`. Default `FALSE`.

- conf.level:

  Confidence level for the interval. Default `0.95`.

- ...:

  Must be empty. Present for compatibility with the generic.

## Value

A tibble with columns

- `.landmark`, `term`

- `estimate`, `std.error`, `statistic`, `p.value`, and `conf.low`,
  `conf.high` when `conf.int = TRUE`

- `origin`: how the term's variables were obtained: `"baseline"`,
  `"locf"`, `"lme"`, or a combination such as `"baseline:lme"` for an
  interaction

- `n_eligible`, `n_model`, `n_event`, `n_competing`

- `.model`, `.cause`, `.fit_status`, `.message`

## Details

Inference is Wald-based, as in
[`summary()`](https://rdrr.io/r/base/summary.html) of a
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) fit:
`statistic` is `estimate / std.error`, `p.value` is two-sided and
confidence limits are `estimate` \\\pm z\_{1-\alpha/2}\\ `std.error` on
the log scale. Fine-Gray models use the robust standard error. Standard
errors treat covariates derived by
[`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md)
as known values, so they are too small for `lme` covariates.

## See also

[`fit_landmark()`](https://thatoneguy006.github.io/tidylandmark/reference/fit_landmark.md)

## Examples

``` r
pbc_lm <- pbcseq_baseline |>
  landmark_data(id, futime, status, landmarks = c(365, 730), horizon = 1825) |>
  add_landmark_covariates(pbcseq_visits, id, day, locf = bili)
#> ℹ Ignored 32 rows of `longitudinal_data` whose `id` is not in `landmark_data`.
fit <- fit_landmark(pbc_lm, ~ age + log(bili), cause = 2)

tidy(fit)
#> # A tibble: 4 × 15
#>   .landmark term      estimate std.error statistic  p.value origin   n_eligible
#>       <dbl> <chr>        <dbl>     <dbl>     <dbl>    <dbl> <chr>         <int>
#> 1       365 age         0.0546    0.0103      5.31 1.09e- 7 baseline        290
#> 2       365 log(bili)   1.29      0.121      10.7  1.61e-26 locf            290
#> 3       730 age         0.0576    0.0103      5.62 1.94e- 8 baseline        278
#> 4       730 log(bili)   1.29      0.119      10.8  3.94e-27 locf            278
#> # ℹ 7 more variables: n_model <int>, n_event <int>, n_competing <int>,
#> #   .model <chr>, .cause <int>, .fit_status <chr>, .message <chr>
tidy(fit, exponentiate = TRUE, conf.int = TRUE)
#> # A tibble: 4 × 17
#>   .landmark term      estimate std.error statistic  p.value conf.low conf.high
#>       <dbl> <chr>        <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1       365 age           1.06    0.0103      5.31 1.09e- 7     1.04      1.08
#> 2       365 log(bili)     3.62    0.121      10.7  1.61e-26     2.85      4.58
#> 3       730 age           1.06    0.0103      5.62 1.94e- 8     1.04      1.08
#> 4       730 log(bili)     3.62    0.119      10.8  3.94e-27     2.87      4.58
#> # ℹ 9 more variables: origin <chr>, n_eligible <int>, n_model <int>,
#> #   n_event <int>, n_competing <int>, .model <chr>, .cause <int>,
#> #   .fit_status <chr>, .message <chr>
```
