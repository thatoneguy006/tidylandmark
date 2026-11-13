
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
*Censored*, `1` = *Received Transplant*, or`2` = *Died*. A sample of
baseline patient characteristics and final status status is shown below:

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

Data preparation for landmarking is fist performed with the
`landmark_data` function, which constructs the various landmark risk
sets and included baseline covariates.

Here, we specify a landmark analysis setup where we aim to estimate from
landmarks at approximately one through five years post-enrollment. In
this example, we are interested in the six-month period after each
landmark time and thus specify `horizon = 180`.

A summary of the landmark analysis set-up is printed along the output,
indicating the landmark times, number at risk at the landmark time, and
the number censored or having events before the horizon, if specified.

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
