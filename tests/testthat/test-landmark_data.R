# Expected values below use landmarks 30/60/90 and horizon 180.

test_that("risk sets, shifted times and statuses match hand-computed values", {
  lm <- lm_spec()
  at <- function(s) lm[lm$.landmark == s, ]

  # s = 30: id 2 (T = 30) is excluded; id 10 (T - s = 240) is censored at 180.
  expect_equal(at(30)$.id, 3:10)
  expect_equal(at(30)$.time, c(18, 35, 55, 80, 115, 160, 180, 180))
  expect_equal(at(30)$.status, c(1, 2, 0, 1, 2, 1, 0, 0))

  # s = 60: ids 9 and 10 exceed the horizon.
  expect_equal(at(60)$.id, 4:10)
  expect_equal(at(60)$.time, c(5, 25, 50, 85, 130, 180, 180))
  expect_equal(at(60)$.status, c(2, 0, 1, 2, 1, 0, 0))

  # s = 90: id 10 has T - s = 180 exactly, so its event is kept.
  expect_equal(at(90)$.id, 6:10)
  expect_equal(at(90)$.time, c(20, 55, 100, 160, 180))
  expect_equal(at(90)$.status, c(1, 2, 1, 0, 1))

  expect_true(all(lm$.time > 0))
  expect_true(all(lm$.horizon == 180))
})

test_that("an event or censoring exactly at the landmark makes a subject ineligible", {
  lm <- lm_spec(landmarks = c(30, 85))
  expect_false(2 %in% lm$.id[lm$.landmark == 30]) # event at 30
  expect_false(5 %in% lm$.id[lm$.landmark == 85]) # censored at 85
})

test_that("horizon = NULL keeps the full follow-up and original status", {
  lm <- lm_spec(landmarks = 30, horizon = NULL)
  expect_equal(lm$.time, subjects$followup_days[3:10] - 30)
  expect_equal(lm$.status, subjects$event[3:10])
  expect_true(all(is.na(lm$.horizon)))
  expect_type(lm$.horizon, "double")
})

test_that("output has the documented structure and types", {
  lm <- lm_spec()
  expect_s3_class(lm, c("landmark_tbl", "tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_named(lm, c(".id", ".landmark", ".time", ".status", ".horizon", "age", "sex"))
  expect_type(lm$.status, "integer")
  expect_type(lm$.time, "double")
  expect_type(lm$.landmark, "double")
  expect_equal(unique(lm$.landmark), c(30, 60, 90))
})

test_that("ID type is preserved", {
  d <- subjects
  d$patient_id <- paste0("P", d$patient_id)
  expect_type(lm_spec(d)$.id, "character")

  d$patient_id <- factor(d$patient_id)
  expect_s3_class(lm_spec(d)$.id, "factor")
})

test_that("landmark order in the input does not matter", {
  expect_identical(lm_spec(landmarks = c(90, 30, 60)), lm_spec())
})

test_that("strings and all_of() select the same columns as bare names", {
  ref <- lm_spec()
  by_string <- landmark_data(
    subjects, "patient_id", "followup_days", "event",
    landmarks = c(30, 60, 90), horizon = 180
  )
  col <- "patient_id"
  by_all_of <- landmark_data(
    subjects, tidyselect::all_of(col), followup_days, event,
    landmarks = c(30, 60, 90), horizon = 180
  )
  expect_identical(by_string, ref)
  expect_identical(by_all_of, ref)
})

test_that("plain data.frame input works", {
  expect_equal(
    tibble::as_tibble(lm_spec(as.data.frame(subjects))),
    tibble::as_tibble(lm_spec())
  )
})

test_that("missing baseline predictors never drop subjects", {
  d <- subjects
  d$age[10] <- NA
  lm <- lm_spec(d)
  expect_equal(sum(lm$.id == 10), 3)
})

test_that("factor status: first level is censoring, later levels become 1, 2, ...", {
  d <- subjects
  d$event <- factor(
    c("censored", "relapse", "death")[subjects$event + 1],
    levels = c("censored", "relapse", "death")
  )
  lm_fac <- lm_spec(d)
  expect_equal(lm_fac$.status, lm_spec()$.status)
  expect_equal(
    landmark_config(lm_fac)$status_labels,
    c(`0` = "censored", `1` = "relapse", `2` = "death")
  )
})

test_that("a landmark with nobody at risk warns, adds no rows and is still counted", {
  expect_warning(lm <- lm_spec(landmarks = c(30, 500)), "No subjects")
  expect_false(500 %in% lm$.landmark)

  counts <- lm_counts(lm)
  expect_equal(counts$landmark, c(30, 500))
  expect_equal(counts$n_at_risk, c(8, 0))
})

test_that("counts match the risk sets", {
  counts <- lm_counts(lm_spec())
  expect_equal(counts$n_at_risk, c(8, 7, 5))
  expect_equal(counts$n_censored, c(3, 3, 1))
  expect_equal(counts$n_event_1, c(3, 2, 3))
  expect_equal(counts$n_event_2, c(2, 2, 1))
})

test_that("config records the mappings and sample sizes", {
  cfg <- landmark_config(lm_spec())
  expect_equal(cfg[c("id", "time", "status")], list(id = "patient_id", time = "followup_days", status = "event"))
  expect_equal(cfg$landmarks, c(30, 60, 90))
  expect_equal(cfg$horizon, 180)
  expect_equal(cfg$event_codes, 1:2)
  expect_equal(cfg$n_subjects, 10)
  expect_equal(unname(cfg$n_eligible), c(8, 7, 5))
})

test_that("print shows the summary and as_tibble() strips the class", {
  lm <- lm_spec()
  expect_output(print(lm), "<landmark_tbl> 3 landmark")
  expect_output(print(lm), "n_event_2")

  plain <- tibble::as_tibble(lm)
  expect_s3_class(plain, c("tbl_df", "tbl", "data.frame"), exact = TRUE)
  expect_null(attr(plain, "landmark_config"))
})

test_that("filtering keeps the class and the summary describes the remaining rows", {
  lm <- lm_spec()
  only_30 <- dplyr::filter(lm, .landmark == 30)

  expect_s3_class(only_30, "landmark_tbl")
  expect_equal(lm_counts(only_30)$landmark, 30)
  expect_equal(lm_counts(only_30)$n_at_risk, 8)

  header <- lm_header(only_30, landmark_config(only_30))
  expect_match(header[1], "1 landmark\\(s\\): 30$")
  expect_match(header[3], "8 at risk in current rows \\(10 in source data\\)")
})

test_that("base `[` subsetting gives the same summary as filter()", {
  lm <- lm_spec()
  expect_equal(
    lm_counts(lm[lm$.landmark == 30, ]),
    lm_counts(dplyr::filter(lm, .landmark == 30))
  )
})

test_that("landmarks that were empty at construction stay visible until a populated landmark is removed", {
  lm <- suppressWarnings(lm_spec(landmarks = c(30, 60, 500)))

  # Row filtering that keeps every populated landmark: 500 still shown as zeros
  older <- dplyr::filter(lm, age > 10)
  expect_equal(lm_counts(older)$landmark, c(30, 60, 500))
  expect_equal(lm_counts(older)$n_at_risk[3], 0)

  # Removing a populated landmark: only what remains is shown
  only_30 <- dplyr::filter(lm, .landmark == 30)
  expect_equal(lm_counts(only_30)$landmark, 30)
})

test_that("invalid `data` and column arguments error", {
  expect_error(landmark_data(list(a = 1), a, a, a, 30), "data frame")
  expect_error(lm_spec(landmarks = NULL), "non-empty numeric")
  expect_error(
    landmark_data(subjects, patient_id, followup_days, event),
    "must be supplied"
  )
  expect_error(
    landmark_data(subjects, nope, followup_days, event, 30),
    "exist"
  )
  expect_error(
    landmark_data(subjects, patient_id, followup_days, followup_days, 30),
    "different columns"
  )

  d <- subjects
  d$.time <- 1
  expect_error(lm_spec(d), "reserved")
})

test_that("id, time and status contents are validated", {
  d <- subjects; d$patient_id[2] <- 1L
  expect_error(lm_spec(d), "uniquely identify")

  d <- subjects; d$patient_id[2] <- NA
  expect_error(lm_spec(d), "missing values")

  d <- subjects; d$followup_days[1] <- NA
  expect_error(lm_spec(d), "missing values")

  d <- subjects; d$followup_days[1] <- -1
  expect_error(lm_spec(d), "non-negative")

  d <- subjects; d$followup_days[1] <- Inf
  expect_error(lm_spec(d), "non-negative")

  d <- subjects; d$followup_days <- as.character(d$followup_days)
  expect_error(lm_spec(d), "numeric column")

  d <- subjects; d$event[1] <- NA
  expect_error(lm_spec(d), "missing values")

  d <- subjects; d$event[1] <- 1.5
  expect_error(lm_spec(d), "positive integer")

  d <- subjects; d$event[1] <- -1
  expect_error(lm_spec(d), "positive integer")

  d <- subjects; d$event <- as.character(d$event)
  expect_error(lm_spec(d), "numeric or factor")

  d <- subjects; d$event <- d$event > 0
  expect_error(lm_spec(d), "numeric or factor")

  d <- subjects; d$event <- factor(rep("censored", 10))
  expect_error(lm_spec(d), "at least two levels")
})

test_that("landmarks are validated", {
  expect_error(lm_spec(landmarks = numeric(0)), "non-empty")
  expect_error(lm_spec(landmarks = "30"), "non-empty numeric")
  expect_error(lm_spec(landmarks = c(30, NA_real_)), "finite")
  expect_error(lm_spec(landmarks = c(30, Inf)), "finite")
  expect_error(lm_spec(landmarks = -1), "non-negative")
  expect_error(lm_spec(landmarks = c(30, 30)), "duplicate")
})

test_that("horizon is validated", {
  expect_error(lm_spec(horizon = c(1, NA)), "single positive")
  expect_error(lm_spec(horizon = NA_real_), "single positive")
  expect_error(lm_spec(horizon = 0), "single positive")
  expect_error(lm_spec(horizon = -5), "single positive")
  expect_error(lm_spec(horizon = Inf), "single positive")
  expect_error(lm_spec(horizon = "180"), "single positive")
})
