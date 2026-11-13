# tidylandmark: Tidy Workflows for Landmark Survival Analysis

Provides a tidy framework for landmark survival analysis, facilitating
the construction of landmark-specific risk sets, incorporation of
longitudinal covariates, and estimation of associations with
time-to-event outcomes. Supports covariate estimation using last
observation carried forward (LOCF) or linear mixed-effects models (LME),
alongside Cox proportional hazards and Fine-Gray competing risks
regression. Offers a consistent and reproducible workflow for conducting
analyses across multiple landmark times.

## Package options

- `tidylandmark.progress`: `TRUE` or `FALSE`. The default for the
  `progress` argument of long-running functions. When unset, progress
  bars are shown only in interactive sessions, so nothing prints inside
  knitr, Quarto or `R CMD check`. Set
  `options(tidylandmark.progress = FALSE)` to turn them off everywhere.

## See also

Useful links:

- <https://thatoneguy006.github.io/tidylandmark/>

- <https://github.com/thatoneguy006/tidylandmark>

- Report bugs at <https://github.com/thatoneguy006/tidylandmark/issues>

## Author

**Maintainer**: Carson Richardson <carson.richardson@outlook.com>
([ORCID](https://orcid.org/0009-0004-8710-1097)) \[copyright holder\]

Authors:

- Carson Richardson <carson.richardson@outlook.com>
  ([ORCID](https://orcid.org/0009-0004-8710-1097)) \[copyright holder\]

- Chance Alvarado ([ORCID](https://orcid.org/0000-0001-7266-6344))
  \[copyright holder\]
