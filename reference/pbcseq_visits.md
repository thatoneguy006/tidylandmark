# Mayo Clinic Primary Biliary Cirrhosis Follow-up Visit Data for Randomized Patients

Adapted from
[survival::pbcseq](https://rdrr.io/pkg/survival/man/pbcseq.html), this
dataset contains follow-up visit data preceding final study status among
randomized patients in the Mayo Clinic primary biliary cirrhosis trial.
Accompanies
[pbcseq_baseline](https://thatoneguy006.github.io/tidylandmark/reference/pbcseq_baseline.md)
which contains baseline and outcome status data.

## Usage

``` r
pbcseq_visits
```

## Format

A data frame with 1,945 observations of 14 variables:

- id:

  case number

- day:

  number of days between enrollment and this visit date all measurements
  below refer to this date

- ascites:

  presence of ascites

- hepato:

  presence of hepatomegaly or enlarged liver

- spiders:

  blood vessel malformations in the skin

- edema:

  0 no edema, 0.5 untreated or successfully treated edema despite
  diuretic therapy

- bili:

  serum bilirunbin (mg/dl)

- chol:

  serum cholesterol (mg/dl)

- albumin:

  serum albumin (mg/dl)

- alk.phos:

  alkaline phosphotase (U/liter)

- ast:

  aspartate aminotransferase, once called SGOT (U/ml)

- platelet:

  platelet count

- protime:

  standardised blood clotting time

- stage:

  histologic stage of disease (needs biopsy)

## Source

Adapted from
[survival::pbcseq](https://rdrr.io/pkg/survival/man/pbcseq.html).

T Therneau and P Grambsch, "Modeling Survival Data: Extending the Cox
Model", Springer-Verlag, New York, 2000. ISBN: 0-387-98784-3.
