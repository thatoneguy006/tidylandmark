# tidylandmark 0.0.0.9000 (development version)

## Major Changes

- Implemented `landmark_data()`, which constructs landmark risk sets with an 
optional user-supplied time horizon.
- Implemented `add_landmark_covariates()`, which attaches pre-landmark values
of longitudinal variables to each subject-landmark row. Variables selected with
`locf` use the last non-missing value (optionally within a `lookback` window);
variables selected with `lme` use subject-specific predictions from a linear
mixed-effects model fitted separately for each variable and landmark, with
options in `lme_args` (random slope, adjustment terms, `lmerControl()`).
- Implemented `fit_landmark()`, which fits an independent Cox, cause-specific,
or Fine-Gray model at each landmark. 
- Implemented `bootstrap_landmark()`, which resamples each landmark's risk set,
re-runs a user-supplied pipeline (mixed models, predictions and survival
models) on every replicate, and reports percentile intervals and p-values that
account for uncertainty in LME-derived covariates.

## Minor Changes

- Added `print()`, `summary()`, and `tidy()` methods to examine output from 
`landmark_data()`.
- Added `landmark_covariates()` and `landmark_lme_fits()` to inspect and 
audit both LOCF and LME fits.
- Added `print()`, `summary()`, and `tidy()` methods for `fit_landmark()`.
