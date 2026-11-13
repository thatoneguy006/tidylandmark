# Shared fixtures. testthat sources every helper-*.R file before running tests.

# Subject-level example.
subjects <- tibble::tibble(
  patient_id    = 1:10,
  followup_days = c(20, 30, 48, 65, 85, 110, 145, 190, 250, 270),
  event         = c(1, 2, 1, 2, 0, 1, 2, 1, 0, 1),
  age           = c(12, 15, 11, 17, 9, 14, 13, 16, 10, 15),
  sex           = rep(c("F", "M"), 5)
)

# Convenience wrapper so tests only state what they vary.
lm_spec <- function(data = subjects, landmarks = c(30, 60, 90), horizon = 180) {
  landmark_data(
    data,
    id = patient_id, time = followup_days, status = event,
    landmarks = landmarks, horizon = horizon
  )
}
