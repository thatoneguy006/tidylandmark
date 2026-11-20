# Changelog

## tidylandmark 0.0.0.9000 (development version)

### Major Changes

- Implemented
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md)
  which constructs landmark risk sets with an optional user-supplied
  time horizon.
- Implemented
  [`add_landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/add_landmark_covariates.md),
  which attaches pre-landmark values of longitudinal variables to each
  subject-landmark row. Variables selected with `locf` use the last
  non-missing value (optionally within a `lookback` window); variables
  selected with `lme` use subject-specific predictions from a linear
  mixed-effects model fitted separately for each variable and landmark,
  with options in `lme_args` (random slope, adjustment terms,
  `lmerControl()`).

### Minor Changes

- Added [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) methods to
  examine output from
  [`landmark_data()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_data.md).
- Implemented
  [`landmark_covariates()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_covariates.md)
  and
  [`landmark_lme_fits()`](https://thatoneguy006.github.io/tidylandmark/reference/landmark_covariates.md)
  to inspect and audit both LOCF and LME fits.
