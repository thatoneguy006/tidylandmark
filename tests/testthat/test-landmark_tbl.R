# tidy.landmark_tbl() ---------------------------------------------------------
# Expected counts use landmarks 30/60/90 and horizon 180 (see
# test-landmark_data.R for the hand-computed risk sets).

# Factor version of the fixture's status: 0 = censored, 1 = relapse, 2 = death.
subjects_fct <- function(levels = c("censored", "relapse", "death")) {
  d <- subjects
  d$event <- factor(levels[subjects$event + 1], levels = levels)
  d
}

test_that("tidy() returns a plain tibble of per-landmark counts", {
  out <- tidy(lm_spec())

  expect_s3_class(out, c("tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_null(attr(out, "landmark_config"))
  expect_named(out, c("landmark", "n_at_risk", "n_censored", "n_event_1", "n_event_2"))
  expect_type(out$landmark, "double")
  expect_true(all(purrr::map_lgl(out[-1], is.integer)))

  expect_equal(out$landmark, c(30, 60, 90))
  expect_equal(out$n_at_risk, c(8L, 7L, 5L))
  expect_equal(out$n_censored, c(3L, 3L, 1L))
  expect_equal(out$n_event_1, c(3L, 2L, 3L))
  expect_equal(out$n_event_2, c(2L, 2L, 1L))
})

test_that("counts add up: at risk = censored + all events", {
  out <- tidy(lm_spec())
  expect_equal(out$n_at_risk, out$n_censored + out$n_event_1 + out$n_event_2)
})

test_that("tidy() matches the counts that print() shows", {
  lm <- lm_spec()
  expect_equal(tidy(lm), lm_counts(lm))
})

test_that("tidy() follows filter() and `[` subsetting", {
  lm <- lm_spec()

  only_30 <- tidy(dplyr::filter(lm, .landmark == 30))
  expect_equal(only_30$landmark, 30)
  expect_equal(only_30$n_at_risk, 8L)

  expect_equal(tidy(lm[lm$.landmark == 30, ]), only_30)

  # Filtering on a baseline column keeps every landmark but changes the counts
  older <- tidy(dplyr::filter(lm, age > 12))
  expect_equal(older$landmark, c(30, 60, 90))
  expect_equal(
    older$n_at_risk,
    as.integer(table(lm$.landmark[lm$age > 12])[c("30", "60", "90")])
  )
})

test_that("event columns are kept, as zeros, when filtering removes a cause", {
  out <- tidy(dplyr::filter(lm_spec(), .status != 2L))
  expect_named(out, c("landmark", "n_at_risk", "n_censored", "n_event_1", "n_event_2"))
  expect_equal(out$n_event_2, c(0L, 0L, 0L))
})

test_that("event columns follow the observed codes", {
  # A single cause gives a single event column
  d <- subjects
  d$event <- as.integer(d$event > 0)
  expect_named(tidy(lm_spec(d)), c("landmark", "n_at_risk", "n_censored", "n_event_1"))

  # Codes need not be consecutive
  d$event <- c(1, 3, 1, 3, 0, 1, 3, 1, 0, 1)
  out <- tidy(lm_spec(d))
  expect_named(out, c("landmark", "n_at_risk", "n_censored", "n_event_1", "n_event_3"))
  expect_equal(out$n_event_3, c(2L, 2L, 1L))
})

test_that("an empty landmark appears as a row of zeros", {
  lm <- suppressWarnings(lm_spec(landmarks = c(30, 500)))
  out <- tidy(lm, show_proportion = TRUE)

  expect_equal(out$landmark, c(30, 500))
  expect_equal(unlist(out[2, -1], use.names = FALSE), rep(0, 5))
})

test_that("show_proportion adds prop_at_risk after landmark, against the source n", {
  lm <- lm_spec()
  out <- tidy(lm, show_proportion = TRUE)

  expect_named(
    out,
    c("landmark", "prop_at_risk", "n_at_risk", "n_censored", "n_event_1", "n_event_2")
  )
  expect_equal(out$prop_at_risk, c(0.8, 0.7, 0.5))

  # The denominator is the source data (10 subjects), not the filtered rows.
  only_30 <- tidy(dplyr::filter(lm, .landmark == 30), show_proportion = TRUE)
  expect_equal(only_30$prop_at_risk, 0.8)
})

test_that("show_proportion does not raise tidyselect deprecation warnings", {
  # `.data$col` inside relocate() is deprecated in tidyselect >= 1.2.0.
  rlang::local_options(lifecycle_verbosity = "warning")
  expect_no_warning(tidy(lm_spec(), show_proportion = TRUE))
})

test_that("status_labels names event columns after the factor levels", {
  lm_fac <- lm_spec(subjects_fct())

  labelled <- tidy(lm_fac, status_labels = TRUE)
  expect_named(labelled, c("landmark", "n_at_risk", "n_censored", "n_relapse", "n_death"))

  # Only the names change
  coded <- tidy(lm_fac)
  expect_equal(unname(as.list(labelled)), unname(as.list(coded)))
  expect_equal(coded, tidy(lm_spec()))
})

test_that("status_labels and show_proportion combine", {
  out <- tidy(lm_spec(subjects_fct()), status_labels = TRUE, show_proportion = TRUE)
  expect_named(
    out,
    c("landmark", "prop_at_risk", "n_at_risk", "n_censored", "n_relapse", "n_death")
  )
})

test_that("status_labels makes level names syntactic", {
  lm_fac <- lm_spec(subjects_fct(c("censored", "liver transplant", "death")))
  expect_named(
    tidy(lm_fac, status_labels = TRUE),
    c("landmark", "n_at_risk", "n_censored", "n_liver.transplant", "n_death")
  )
})

test_that("status_labels warns and keeps code names when status was numeric", {
  lm <- lm_spec()
  expect_warning(out <- tidy(lm, status_labels = TRUE), "no effect")
  expect_equal(out, tidy(lm))
})

test_that("status_labels errors instead of creating duplicate column names", {
  # An event level called "censored" would collide with n_censored
  lm_clash <- lm_spec(subjects_fct(c("alive", "censored", "death")))
  expect_error(tidy(lm_clash, status_labels = TRUE), "not be unique")
  expect_no_error(tidy(lm_clash))

  # Two levels that make.names() maps to the same name
  lm_dup <- lm_spec(subjects_fct(c("censored", "a b", "a.b")))
  expect_error(tidy(lm_dup, status_labels = TRUE), "n_a.b")
})

test_that("tidy() validates its arguments", {
  lm <- lm_spec()
  expect_error(tidy(lm, status_labels = NA), "TRUE.*FALSE")
  expect_error(tidy(lm, status_labels = "yes"), "TRUE.*FALSE")
  expect_error(tidy(lm, show_proportion = c(TRUE, TRUE)), "TRUE.*FALSE")
  expect_error(tidy(lm, show_proportion = 1), "TRUE.*FALSE")

  # A misspelled argument is caught rather than silently ignored
  expect_error(tidy(lm, show_proportions = TRUE), "must be empty")
})

test_that("tidy() errors informatively on a broken landmark_tbl", {
  lm <- lm_spec()

  # Dropping a reserved column used to return `stats::df` instead of a table
  expect_error(tidy(dplyr::select(lm, -.status)), "\\.status")

  no_cfg <- lm
  attr(no_cfg, "landmark_config") <- NULL
  expect_error(tidy(no_cfg), "landmark_config")
})
