# Two-stage bootstrap for propagating uncertainty from LME predictions

`bootstrap_landmark()` bootstraps the complete landmark workflow. In
each replicate, the risk set at every landmark is resampled with
replacement, the user's `pipeline` is re-run on the resampled data
(refitting any longitudinal mixed-effects models, regenerating their
predictions at each landmark and refitting the landmark survival
models), and the coefficients are collected.

The bootstrap distribution quantifies uncertainty from estimating the
longitudinal trajectories and its propagation into the landmark survival
models, which the standard errors from
[`fit_landmark()`](https://thatoneguy006.github.io/tidylandmark/reference/fit_landmark.md)
ignore. Results are reported as percentile confidence intervals,
percentile bootstrap p-values and the median of the bootstrap
distribution, next to the full-data fit.

## Usage

``` r
bootstrap_landmark(
  baseline,
  pipeline,
  n_boot = 200,
  exponentiate = FALSE,
  conf.level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  progress = getOption("tidylandmark.progress", interactive())
)
```

## Arguments

- baseline:

  A `landmark_tbl`, as returned by
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md).

- pipeline:

  A function of one argument, a `landmark_tbl`, that runs the rest of
  the workflow (e.g.
  [`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md)
  then
  [`fit_landmark()`](https://thatoneguy006.github.io/tidylandmark/reference/fit_landmark.md))
  and returns a `landmark_fit`.

- n_boot:

  Number of bootstrap replicates. Default `200`.

- exponentiate:

  If `TRUE`, report hazard ratios (subdistribution hazard ratios for
  Fine-Gray models) instead of log hazard ratios.

- conf.level:

  Confidence level for the intervals. Default `0.95`.

- alternative:

  Alternative hypothesis for the bootstrap p-values: one of
  `"two.sided"` (default), `"greater"` or `"less"`.

- progress:

  Show a progress bar over replicates? Defaults to the
  `tidylandmark.progress` option, or
  [`interactive()`](https://rdrr.io/r/base/interactive.html) if it is
  unset.

## Value

A tibble with one row per landmark and term of the full-data fit:

- `.landmark`, `.model`, `.cause`, `origin`, `term`: as in
  [`tidy.landmark_fit()`](https://thatoneguy006.github.io/tidylandmark/reference/tidy.landmark_fit.md).

- `estimate__full`, `conf.low__full`, `conf.high__full`,
  `p.value__full`: the full-data fit with model-based (Wald) inference.

- `estimate__boot`, `conf.low__boot`, `conf.high__boot`,
  `p.value__boot`: median, percentile interval and percentile p-value of
  the bootstrap distribution.

- `n_boot_ok`: replicates with an estimate for this term.

- `.fit_status__full`, `.message__full`: fit status of the full-data
  fit.

- `.estimates__boot`: one estimate per replicate, in draw order.

- `.fit_status__boot`, `.message__boot`: list-columns with the fit
  status and message for the same replicate, landmark and term.

## Details

### Resampling

Rows are resampled independently within each landmark, which matches
[`fit_landmark()`](https://thatoneguy006.github.io/tidylandmark/reference/fit_landmark.md)
fitting one independent model per landmark and
[`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md)
fitting one mixed model per landmark from that landmark's risk set. The
bootstrap therefore supports per-landmark inference, but not comparisons
between landmarks. A subject drawn more than once is treated as distinct
subjects, each with its own random effects in the mixed models.

Risk sets are read from the rows of `baseline` as they are now, so a
`landmark_tbl` filtered after
[`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md)
is resampled as filtered. Create any LME predictors used by the survival
model inside `pipeline` so their mixed models are refitted on each draw.
For spline coefficients, set the interior and boundary knots before
resampling and pass them to the spline function. The bootstrap checks
the fitted basis at each landmark.

### Failed replicates

Each replicate runs quietly: progress bars (when `pipeline` uses the
default `progress`), messages and warnings are suppressed, and a
replicate whose `pipeline` errors is recorded as `"failed"` at every
landmark instead of stopping the bootstrap. `n_boot_ok` gives the number
of replicates with an estimate for each term; summaries use only those.
Check `.fit_status__boot` when `n_boot_ok` is well below `n_boot`. An
unavailable coefficient (for example, a comparison to an absent factor
level) has an `NA` estimate in that draw. Intervals and p-values require
at least two usable estimates.

### Reproducibility and p-values

Results depend on the random number generator; call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible results. A zero p-value means that no usable draw fell in
the relevant tail. Its Monte Carlo precision depends on `n_boot_ok`.

## See also

[`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md),
[`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md),
[`fit_landmark()`](https://thatoneguy006.github.io/tidylandmark/reference/fit_landmark.md),
[`tidy.landmark_fit()`](https://thatoneguy006.github.io/tidylandmark/reference/tidy.landmark_fit.md)

## Examples

``` r
# \donttest{
# PBC: status 0 = censored, 1 = transplant, 2 = death
set.seed(2026)
pbcseq_baseline |>
  landmark_data(id, futime, status, landmarks = c(365, 730), horizon = 1825) |>
  bootstrap_landmark(
    pipeline = \(data) {
      data |>
        add_landmark_covariates(pbcseq_visits, id, day, lme = c(bili, albumin)) |>
        fit_landmark(~ age + log(bili) + albumin, cause = 2)
    },
    n_boot = 50
  )
#> ℹ Ignored 32 rows of `longitudinal_data` whose `id` is not in `landmark_data`.
#> ℹ Standard errors treat LME-derived covariates as known, so they understate
#>   uncertainty.
#>   Use `bootstrap_landmark()` for intervals that account for this.
#> This message is displayed once per session.
#> # A tibble: 6 × 19
#>   .landmark .model             .cause origin term  estimate__full conf.low__full
#>       <dbl> <chr>               <int> <chr>  <chr>          <dbl>          <dbl>
#> 1       365 cause-specific cox      2 basel… age           0.0365         0.0157
#> 2       365 cause-specific cox      2 lme    log(…         1.31           1.02  
#> 3       365 cause-specific cox      2 lme    albu…        -2.28          -3.66  
#> 4       730 cause-specific cox      2 basel… age           0.0387         0.0173
#> 5       730 cause-specific cox      2 lme    log(…         1.96           1.57  
#> 6       730 cause-specific cox      2 lme    albu…        -1.72          -2.75  
#> # ℹ 12 more variables: conf.high__full <dbl>, p.value__full <dbl>,
#> #   estimate__boot <dbl>, conf.low__boot <dbl>, conf.high__boot <dbl>,
#> #   p.value__boot <dbl>, n_boot_ok <int>, .fit_status__full <chr>,
#> #   .message__full <chr>, .estimates__boot <list>, .fit_status__boot <list>,
#> #   .message__boot <list>
# }
```
