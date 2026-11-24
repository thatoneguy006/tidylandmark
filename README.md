
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidylandmark

<!-- badges: start -->

[![R-CMD-check](https://github.com/thatoneguy006/tidylandmark/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/thatoneguy006/tidylandmark/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/thatoneguy006/tidylandmark/graph/badge.svg)](https://app.codecov.io/gh/thatoneguy006/tidylandmark)

<!-- badges: end -->

Provides a tidy framework for landmark survival analysis, facilitating
the construction of landmark-specific risk sets, incorporation of
longitudinal covariates, and estimation of associations with
time-to-event outcomes. Supports covariate estimation using last
observation carried forward (LOCF) or linear mixed-effects models (LME),
alongside Cox proportional hazards and Fine-Gray competing risks
regression. Offers a consistent and reproducible workflow for conducting
analyses across multiple landmark times.

## Installation

### CRAN

Coming soon…

### Development

You can install the development version of tidylandmark from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("thatoneguy006/tidylandmark")
```

## Use

We illustrate the package’s functionality with a basic example using
data from the Mayo Clinic primary biliary cirrhosis trial. This data was
originally bundled with the `survival` package as
[`survival::pbcseq`](https://rdrr.io/cran/survival/man/pbcseq.html). We
have included a copy of this data split across two datasets,
`pbcseq_baseline` and `pbcseq_visits` for baseline/status and follow-up
visits, respectively.

The data files include information on **312** randomized patients.
Patients are followed after enrollment and have statuses of `0` =
*Censored*, `1` = *Received Transplant*, or `2` = *Died*. A sample of
baseline patient characteristics and final status is shown below:

``` r
library(tidylandmark)
head(pbcseq_baseline)
#>   id trt      age sex futime status
#> 1  1   1 58.76523   f    400      2
#> 2  2   1 56.44627   f   5169      0
#> 3  3   1 70.07255   m   1012      2
#> 4  4   1 54.74059   f   1925      2
#> 5  5   0 38.10541   f   1505      1
#> 6  6   0 66.25873   f   2503      2
```

Laboratory values and clinical statuses monitored longitudinally are
included in `pbcseq_visits`. A total of **1,945** observations for the
**312** randomized patients are included, with a sample shown below:

``` r
head(pbcseq_visits)
#>   id day ascites hepato spiders edema bili chol albumin alk.phos   ast platelet
#> 1  1   0       1      1       1     1 14.5  261    2.60     1718 138.0      190
#> 2  1 192       1      1       1     1 21.3   NA    2.94     1612   6.2      183
#> 3  2   0       0      1       1     0  1.1  302    4.14     7395 113.5      221
#> 4  2 182       0      1       1     0  0.8   NA    3.60     2107 139.5      188
#> 5  2 365       0      1       1     0  1.0   NA    3.55     1711 144.2      161
#> 6  2 768       0      1       1     0  1.9   NA    3.92     1365 144.2      122
#>   protime stage
#> 1    12.2     4
#> 2    11.2     4
#> 3    10.6     3
#> 4    11.0     3
#> 5    11.6     3
#> 6    10.6     3
```

### Construct Landmark Risk Sets

Data preparation for landmarking is first performed with the
`landmark_data` function, which constructs the various landmark risk
sets and included baseline covariates.

Here, we specify a landmark analysis setup where we aim to estimate from
landmarks at approximately one through five years post-enrollment. In
this example, we are interested in the six-month period after each
landmark time and thus specify `horizon = 180`.

A summary of the landmark analysis set-up is printed along with the
output, indicating the landmark times, number at risk at the landmark
time, and the number censored or having events before the horizon, if
specified.

``` r
# Construct landmark risk sets
lm_dat <- pbcseq_baseline |>
  landmark_data(
    id = id,
    time = futime,
    status = status,
    landmarks = 365 * 1:5,
    horizon = 180
  )
```

``` r
# Print summary of risk set specifications
lm_dat
#> <landmark_tbl> 5 landmark(s): 365, 730, 1095, 1460, 1825
#> Horizon: 180 (relative to each landmark)
#> Subjects: 290 at risk in current rows (312 in source data)
#> Source columns: id = id, time = futime, status = status
#>  landmark n_at_risk n_censored n_event_1 n_event_2
#>       365       290        285         1         4
#>       730       278        260         4        14
#>      1095       245        237         0         8
#>      1460       225        218         2         5
#>      1825       202        194         3         5
#> 
#> # A tibble: 1,240 × 8
#>      .id .landmark .time .status .horizon   trt   age sex  
#>    <int>     <dbl> <dbl>   <int>    <dbl> <int> <dbl> <fct>
#>  1     1       365    35       2      180     1  58.8 f    
#>  2     2       365   180       0      180     1  56.4 f    
#>  3     3       365   180       0      180     1  70.1 m    
#>  4     4       365   180       0      180     1  54.7 f    
#>  5     5       365   180       0      180     0  38.1 f    
#>  6     6       365   180       0      180     0  66.3 f    
#>  7     7       365   180       0      180     0  55.5 f    
#>  8     8       365   180       0      180     0  53.1 f    
#>  9     9       365   180       0      180     1  42.5 f    
#> 10    11       365   180       0      180     0  53.7 f    
#> # ℹ 1,230 more rows
```

A `tidy` method is also available yielding a tibble containing landmark
risk set summaries.

``` r
# Pull risk set summaries as a tibble
lm_dat |>
  tidy(show_proportion = TRUE)
#> # A tibble: 5 × 6
#>   landmark prop_at_risk n_at_risk n_censored n_event_1 n_event_2
#>      <dbl>        <dbl>     <int>      <int>     <int>     <int>
#> 1      365        0.929       290        285         1         4
#> 2      730        0.891       278        260         4        14
#> 3     1095        0.785       245        237         0         8
#> 4     1460        0.721       225        218         2         5
#> 5     1825        0.647       202        194         3         5
```

### Add Longitudinal Covariates

`add_landmark_covariates()` attaches each patient’s pre-landmark
measurements from `pbcseq_visits` to every landmark risk set. Variables
listed in `locf` use the last observation carried forward. Variables
listed in `lme` use a subject-specific prediction from a linear
mixed-effects model (`lme4::lmer()`), fitted separately at each
landmark.

``` r
# Add longitudinal covariates to the landmark risk sets
lm_cov <- lm_dat |>
  add_landmark_covariates(
    longitudinal_data = dplyr::mutate(pbcseq_visits, log_bili = log(bili)),
    id = id,
    time = day,
    locf = c(ascites, edema),
    lme = c(log_bili, albumin),
    lme_args = list(adjust = ~ age)
  )
#> ℹ Ignored 32 rows of `longitudinal_data` whose `id` is not in `landmark_data`.
```

``` r
# Inspect the risk sets
lm_cov
#> <landmark_tbl> 5 landmark(s): 365, 730, 1095, 1460, 1825
#> Horizon: 180 (relative to each landmark)
#> Subjects: 290 at risk in current rows (312 in source data)
#> Source columns: id = id, time = futime, status = status
#> Derived covariates: ascites (LOCF), edema (LOCF), log_bili (LME, adjusted ~age), albumin (LME, adjusted ~age)
#>  landmark n_at_risk n_censored n_event_1 n_event_2
#>       365       290        285         1         4
#>       730       278        260         4        14
#>      1095       245        237         0         8
#>      1460       225        218         2         5
#>      1825       202        194         3         5
#> 
#> # A tibble: 1,240 × 12
#>      .id .landmark .time .status .horizon   trt   age sex   ascites edema
#>    <int>     <dbl> <dbl>   <int>    <dbl> <int> <dbl> <fct>   <int> <dbl>
#>  1     1       365    35       2      180     1  58.8 f           1   1  
#>  2     2       365   180       0      180     1  56.4 f           0   0  
#>  3     3       365   180       0      180     1  70.1 m           0   0.5
#>  4     4       365   180       0      180     1  54.7 f           0   0.5
#>  5     5       365   180       0      180     0  38.1 f           0   0  
#>  6     6       365   180       0      180     0  66.3 f           0   0  
#>  7     7       365   180       0      180     0  55.5 f           0   0  
#>  8     8       365   180       0      180     0  53.1 f           0   0  
#>  9     9       365   180       0      180     1  42.5 f           0   0  
#> 10    11       365   180       0      180     0  53.7 f           0   0  
#> # ℹ 1,230 more rows
#> # ℹ 2 more variables: log_bili <dbl>, albumin <dbl>
```

### Fit the Landmark Models

`fit_landmark()` fits a separate survival model at each landmark. Supply
a one-sided survival formula; the outcome is built internally from each
risk set’s follow-up time and status. Here we model death (`cause = 2`).
Because transplant is a competing event, the default Cox model estimates
cause-specific hazard ratios.

``` r
# Fit a cause-specific Cox model at every landmark
fit <- lm_cov |>
  fit_landmark(
    formula = ~ age + log_bili + albumin + edema,
    cause = 2
  )

fit
#> <landmark_fit> Cause-specific Cox model for cause 2 at 5 landmark(s)
#> Formula: ~age + log_bili + albumin + edema
#> Horizon: 180 (relative to each landmark)
#>  landmark n_eligible n_model n_event n_competing fit_status
#>       365        290     290       4           1         ok
#>       730        278     278      14           4         ok
#>      1095        245     245       8           0         ok
#>      1460        225     225       5           2         ok
#>      1825        202     202       5           3         ok
#> 
#> Standard errors treat LME-derived covariates (log_bili, albumin) as known.
#> Use tidy() for coefficients and summary() for hazard ratios.
```

`tidy()` returns one row per landmark and term, with the number of
events behind each estimate and the `origin` of each term (baseline,
LOCF or LME). Landmarks that could not be estimated keep a row with
their status, so they don’t silently drop out of the results.

``` r
# Hazard ratios with 95% confidence intervals
fit |>
  tidy(exponentiate = TRUE, conf.int = TRUE) |>
  dplyr::select(.landmark, term, estimate, 
                conf.low, conf.high, origin, n_event) |> 
  head()
#> ℹ Standard errors treat LME-derived covariates as known, so they understate
#>   uncertainty.
#>   Bootstrap intervals that account for this are planned for a future release.
#> This message is displayed once per session.
#> # A tibble: 6 × 7
#>   .landmark term     estimate conf.low conf.high origin   n_event
#>       <dbl> <chr>       <dbl>    <dbl>     <dbl> <chr>      <int>
#> 1       365 age         1.07  0.971         1.18 baseline       4
#> 2       365 log_bili    1.63  0.557         4.79 lme            4
#> 3       365 albumin     0.128 0.000182     90.0  lme            4
#> 4       365 edema      12.4   0.731       209.   locf           4
#> 5       730 age         1.06  0.999         1.11 baseline      14
#> 6       730 log_bili    3.67  1.98          6.80 lme           14
```

Fine-Gray subdistribution hazard models use `model = "finegray"`, and
`summary()` is a shortcut for the exponentiated `tidy()` output. The
native `survival::coxph()` fit for each landmark is kept in
`fit$models`.

``` r
# Fine-Gray model for death, with transplant as a competing event
lm_cov |>
  fit_landmark(~ age + log_bili + albumin + edema, 
               model = "finegray", cause = 2) |>
  summary() |>
  dplyr::select(.landmark, term, estimate, conf.low, conf.high, .fit_status) |> 
  head()
#> Warning: 1 of 5 landmark models did not fit cleanly.
#> • Landmark 1095 (warning): no competing events: estimates equal the Cox model
#> ℹ Unestimable and failed landmarks have no estimates. See `$model_info` for
#>   details.
#> # A tibble: 6 × 6
#>   .landmark term     estimate conf.low conf.high .fit_status
#>       <dbl> <chr>       <dbl>    <dbl>     <dbl> <chr>      
#> 1       365 age         1.07  0.974         1.17 ok         
#> 2       365 log_bili    1.63  0.910         2.93 ok         
#> 3       365 albumin     0.128 0.000370     44.3  ok         
#> 4       365 edema      12.4   0.514       297.   ok         
#> 5       730 age         1.06  0.995         1.12 ok         
#> 6       730 log_bili    3.67  2.08          6.48 ok
```
