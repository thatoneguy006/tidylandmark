# Mayo Clinic Primary Biliary Cirrhosis Baseline and Status Data for Randomized Patients

Adapted from
[survival::pbcseq](https://rdrr.io/pkg/survival/man/pbcseq.html), this
dataset contains baseline and outcome status data among randomized
patients in the Mayo Clinic primary biliary cirrhosis trial. Accompanies
[pbcseq_visits](https://thatoneguy006.github.io/tidylandmark/reference/pbcseq_visits.md)
which contains follow-up visit data preceding `status`.

## Usage

``` r
pbcseq_baseline
```

## Format

A data frame with 312 observations of six variables:

- id:

  case number

- age:

  in years

- sex:

  m/f

- trt:

  1/2 for D-penicillmain, placebo

- futime:

  Time at `status`

- status:

  status at endpoint, 0/1/2 for censored, transplant, dead

## Source

Adapted from
[survival::pbcseq](https://rdrr.io/pkg/survival/man/pbcseq.html).

T Therneau and P Grambsch, "Modeling Survival Data: Extending the Cox
Model", Springer-Verlag, New York, 2000. ISBN: 0-387-98784-3.
