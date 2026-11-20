# Inspect how derived covariates were made

`landmark_covariates()` returns one row per covariate added by
[`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md),
with the method and options used and how complete it is in the current
rows.

`landmark_lme_fits()` returns one row per mixed model (each `lme`
variable at each landmark), with fit status, sample sizes, variance
components and the fitted
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) object in the
`model` column.

Only covariates still present in `x` are included, so columns removed
with
[`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
drop out.

## Usage

``` r
landmark_covariates(x)

landmark_lme_fits(x)
```

## Arguments

- x:

  A `landmark_tbl` returned by
  [`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md).

## Value

`landmark_covariates()`: a tibble with columns `column`, `source`,
`method`, `include_at_landmark`, `lookback`, `random_slope`, `adjust`,
`n_rows`, `n_missing`, and, for LOCF, `median_lag` and `max_lag` (time
between each carried-forward measurement and its landmark). Counts and
lags are computed on the current rows, so they follow
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).

`landmark_lme_fits()`: a tibble with columns `column`, `source`,
`.landmark`, `status` (`"ok"`, `"warning"` or `"failed"`), `message`,
`n_eligible`, `n_subjects`, `n_obs`, `median_obs`, `time_min`,
`time_max`, `time_scale`, `slope` (per original time unit),
`var_intercept`, `var_slope`, `var_residual`, `singular` and `model`
(`NULL` for failed fits). In each model, time (`.t`) is centred at the
landmark and divided by `time_scale`.

Both return a zero-row tibble with the same columns when there is
nothing to report.

## See also

[`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md)

## Examples

``` r
visits_pbc <- transform(pbcseq_visits, log_bili = log(bili))
pbc_lm <- pbcseq_baseline |>
  landmark_data(id, futime, status, landmarks = c(365, 730), horizon = 1825) |>
  add_landmark_covariates(
    visits_pbc, id, day,
    locf = albumin, lme = log_bili, lme_args = list(adjust = ~ age)
  )
#> ℹ Ignored 32 rows of `longitudinal_data` whose `id` is not in `landmark_data`.

landmark_covariates(pbc_lm)
#> # A tibble: 2 × 11
#>   column   source method include_at_landmark lookback random_slope adjust n_rows
#>   <chr>    <chr>  <chr>  <lgl>                  <dbl> <lgl>        <chr>   <int>
#> 1 albumin  album… locf   TRUE                      NA NA           NA        568
#> 2 log_bili log_b… lme    TRUE                      NA FALSE        ~age      568
#> # ℹ 3 more variables: n_missing <int>, median_lag <dbl>, max_lag <dbl>

fits <- landmark_lme_fits(pbc_lm)
fits
#> # A tibble: 2 × 18
#>   column  source .landmark status message n_eligible n_subjects n_obs median_obs
#>   <chr>   <chr>      <dbl> <chr>  <chr>        <int>      <int> <int>      <dbl>
#> 1 log_bi… log_b…       365 ok     NA             290        290   642          2
#> 2 log_bi… log_b…       730 ok     NA             278        278   845          3
#> # ℹ 9 more variables: time_min <dbl>, time_max <dbl>, time_scale <dbl>,
#> #   slope <dbl>, var_intercept <dbl>, var_slope <dbl>, var_residual <dbl>,
#> #   singular <lgl>, model <list>
summary(fits$model[[1]])
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: .value ~ .t + age + (1 | .id)
#>    Data: data
#> Control: control
#> 
#> REML criterion at convergence: 1304.3
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.8287 -0.4163  0.0110  0.3826  3.5283 
#> 
#> Random effects:
#>  Groups   Name        Variance Std.Dev.
#>  .id      (Intercept) 0.8505   0.9222  
#>  Residual             0.1277   0.3574  
#> Number of obs: 642, groups:  .id, 290
#> 
#> Fixed effects:
#>              Estimate Std. Error t value
#> (Intercept)  0.702170   0.272507   2.577
#> .t          -0.007086   0.015214  -0.466
#> age         -0.005363   0.005360  -1.000
#> 
#> Correlation of Fixed Effects:
#>     (Intr) .t    
#> .t   0.096       
#> age -0.973  0.010
```
