# Fixtures: lm_cr() (competing risks, landmarks 50/200/400, n = 400) and
# sim_lme() (mixed-model data, landmarks 100/200, n = 30) from the helper
# files; lm_spec() / `subjects` (10 subjects) from helper-fixtures.R.

fit_q  <- function(...) fit_landmark(..., progress = FALSE)
boot_q <- function(...) bootstrap_landmark(..., progress = FALSE)
cox_pipe <- \(d) fit_q(d, ~ age + sex)

# Run bootstrap_landmark() while recording every landmark_tbl the pipeline
# receives. The first element is the full data, the rest are the replicates.
boot_seen <- function(baseline, n_boot, formula = ~ age + sex) {
  seen <- list()
  pipe <- function(d) {
    seen[[length(seen) + 1L]] <<- d
    fit_q(d, formula)
  }
  out <- boot_q(baseline, pipe, n_boot = n_boot)
  list(out = out, seen = seen)
}

# Validation ---------------------------------------------------------------------

test_that("bootstrap_landmark() validates its arguments", {
  lm <- lm_cr()

  expect_error(boot_q(tibble::as_tibble(lm), cox_pipe), "landmark_tbl")
  expect_error(boot_q(lm, "not a function"), "pipeline")
  expect_error(boot_q(lm, cox_pipe, n_boot = 1), "n_boot")
  expect_error(boot_q(lm, cox_pipe, n_boot = 2.5), "n_boot")
  expect_error(boot_q(lm, cox_pipe, n_boot = NA), "n_boot")
  expect_error(boot_q(lm, cox_pipe, exponentiate = NA), "exponentiate")
  expect_error(boot_q(lm, cox_pipe, conf.level = 1.5), "conf.level")
  expect_error(boot_q(lm, cox_pipe, alternative = "foo"), "alternative")
  expect_error(bootstrap_landmark(lm, cox_pipe, progress = "yes"), "progress")
  expect_error(boot_q(lm, \(d) "nope"), "landmark_fit")
})

test_that("validation errors are reported for bootstrap_landmark()", {
  err <- expect_error(boot_q(lm_cr(), cox_pipe, exponentiate = "yes"))
  expect_equal(rlang::call_name(err$call), "bootstrap_landmark")
})

# Output structure -----------------------------------------------------------------

test_that("output has one row per full-data term and the documented columns", {
  lm <- lm_cr()
  set.seed(1)
  out <- boot_q(lm, cox_pipe, n_boot = 5)
  full <- tidy(cox_pipe(lm), conf.int = TRUE)

  expect_s3_class(out, "tbl_df")
  expect_named(out, c(
    ".landmark", ".model", ".cause", "origin", "term",
    "estimate__full", "conf.low__full", "conf.high__full", "p.value__full",
    "estimate__boot", "conf.low__boot", "conf.high__boot", "p.value__boot",
    "n_boot_ok",
    ".fit_status__full", ".message__full",
    ".estimates__boot", ".fit_status__boot", ".message__boot"
  ))
  expect_equal(out$.landmark, full$.landmark)
  expect_equal(out$term, full$term)

  # The __full columns are exactly tidy() on the full data
  expect_equal(out$estimate__full, full$estimate)
  expect_equal(out$conf.low__full, full$conf.low)
  expect_equal(out$conf.high__full, full$conf.high)
  expect_equal(out$p.value__full, full$p.value)
  expect_equal(out$.fit_status__full, full$.fit_status)

  expect_equal(out$n_boot_ok, rep(5L, nrow(out)))
  expect_true(all(lengths(out$.estimates__boot) == 5L))
  expect_true(all(lengths(out$.fit_status__boot) == 5L))
  expect_true(all(lengths(out$.message__boot) == 5L))
  expect_true(all(out$conf.low__boot <= out$estimate__boot))
  expect_true(all(out$estimate__boot <= out$conf.high__boot))
  expect_true(all(out$p.value__boot >= 0 & out$p.value__boot <= 1))
})

test_that("bootstrap summaries are the median, percentiles and p-value of the estimates", {
  set.seed(2)
  out <- boot_q(lm_cr(), cox_pipe, n_boot = 20, conf.level = 0.9)

  for (i in seq_len(nrow(out))) {
    est <- out$.estimates__boot[[i]]
    expect_equal(out$estimate__boot[i], stats::median(est))
    expect_equal(out$conf.low__boot[i], stats::quantile(est, 0.05, names = FALSE))
    expect_equal(out$conf.high__boot[i], stats::quantile(est, 0.95, names = FALSE))
    expect_equal(out$p.value__boot[i], min(1, 2 * min(mean(est <= 0), mean(est >= 0))))
  }
})

test_that("results are reproducible with set.seed()", {
  lm <- lm_cr()
  set.seed(3)
  a <- boot_q(lm, cox_pipe, n_boot = 4)
  set.seed(3)
  b <- boot_q(lm, cox_pipe, n_boot = 4)
  expect_identical(a, b)
})

test_that("exponentiate = TRUE transforms estimates but not p-values", {
  lm <- lm_cr()
  set.seed(4)
  log_scale <- boot_q(lm, cox_pipe, n_boot = 10)
  set.seed(4)
  hr_scale <- boot_q(lm, cox_pipe, n_boot = 10, exponentiate = TRUE)

  expect_equal(hr_scale$estimate__full, exp(log_scale$estimate__full))
  expect_equal(hr_scale$estimate__boot, exp(log_scale$estimate__boot))
  expect_equal(hr_scale$conf.low__boot, exp(log_scale$conf.low__boot))
  expect_equal(hr_scale$conf.high__boot, exp(log_scale$conf.high__boot))
  expect_equal(hr_scale$.estimates__boot, purrr::map(log_scale$.estimates__boot, exp))
  expect_equal(hr_scale$p.value__boot, log_scale$p.value__boot)
})

test_that("alternative picks the one-sided tail", {
  lm <- lm_cr()
  set.seed(5)
  two <- boot_q(lm, cox_pipe, n_boot = 10)
  set.seed(5)
  greater <- boot_q(lm, cox_pipe, n_boot = 10, alternative = "greater")
  set.seed(5)
  less <- boot_q(lm, cox_pipe, n_boot = 10, alternative = "less")

  expect_equal(greater$p.value__boot, purrr::map_dbl(greater$.estimates__boot, \(e) mean(e <= 0)))
  expect_equal(less$p.value__boot, purrr::map_dbl(less$.estimates__boot, \(e) mean(e >= 0)))
  expect_equal(two$p.value__boot, pmin(1, 2 * pmin(greater$p.value__boot, less$p.value__boot)))
})

# Resampling ---------------------------------------------------------------------

test_that("each replicate resamples rows within landmarks, with replacement", {
  lm <- lm_cr()
  set.seed(6)
  res <- boot_seen(lm, n_boot = 3)
  expect_length(res$seen, 4L)
  expect_identical(res$seen[[1]], lm)

  ref_key <- paste(lm$.id, lm$.landmark)
  for (d in res$seen[-1]) {
    expect_s3_class(d, "landmark_tbl")
    expect_identical(landmark_config(d), landmark_config(lm))
    # Same risk-set size at every landmark
    expect_equal(table(d$.landmark), table(lm$.landmark))
    # Every row is a real row of the baseline at the same landmark
    key <- paste(d$.id, d$.landmark)
    expect_true(all(key %in% ref_key))
    expect_equal(d$.time, lm$.time[match(key, ref_key)])
    expect_equal(d$.status, lm$.status[match(key, ref_key)])
    # Sampling is with replacement
    expect_gt(anyDuplicated(key), 0L)
  }
})

test_that("a filtered baseline is resampled as it now stands", {
  lm <- dplyr::filter(lm_cr(), .landmark == 200, age > 40)
  set.seed(7)
  res <- boot_seen(lm, n_boot = 2)

  expect_equal(unique(res$out$.landmark), 200)
  for (d in res$seen[-1]) {
    expect_equal(unique(d$.landmark), 200)
    expect_true(all(d$.id %in% lm$.id))
    expect_false(anyNA(d$.time))
    expect_false(anyNA(d$.status))
  }
})

test_that("a landmark with a single subject keeps that subject", {
  # Only subject 10 (follow-up 270) is still at risk at day 260. sample(10L, 1)
  # would draw from 1:10, so this catches the classic sample() pitfall.
  lm <- landmark_data(
    subjects, id = patient_id, time = followup_days, status = event,
    landmarks = c(30, 260)
  )
  expect_equal(lm$.id[lm$.landmark == 260], 10L)

  set.seed(8)
  res <- suppressWarnings(boot_seen(lm, n_boot = 5, formula = ~ age))
  for (d in res$seen[-1]) {
    expect_equal(d$.id[d$.landmark == 260], 10L)
  }

  # The model at 260 is unestimable (one subject), so the bootstrap has
  # nothing to summarise there.
  at_260 <- res$out[res$out$.landmark == 260, ]
  expect_equal(at_260$n_boot_ok, 0L)
  expect_true(is.na(at_260$estimate__boot))
})

# Failed and noisy replicates -----------------------------------------------------

test_that("an error in one replicate is recorded instead of stopping the bootstrap", {
  calls <- 0L
  pipe <- function(d) {
    calls <<- calls + 1L
    # Call 1 is the full-data fit; replicates 1, 3 and 5 (calls 2, 4, 6) fail
    if (calls %% 2L == 0L) stop("boom")
    cox_pipe(d)
  }
  set.seed(9)
  expect_no_error(out <- boot_q(lm_cr(), pipe, n_boot = 6))

  expect_equal(out$n_boot_ok, rep(3L, nrow(out)))
  st  <- out$.fit_status__boot[[1]]
  msg <- out$.message__boot[[1]]
  expect_length(st, 6L)
  expect_equal(st == "failed", c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
  expect_equal(unique(msg[st == "failed"]), "boom")
})

test_that("an error in the full-data fit is not caught", {
  expect_error(boot_q(lm_cr(), \(d) stop("full-data failure")), "full-data failure")
})

test_that("replicates are silent: only the full-data fit's conditions surface", {
  pipe <- function(d) {
    message("hello")
    warning("careful")
    cox_pipe(d)
  }
  set.seed(10)
  # Exactly one message and one warning (from the full-data fit); any from the
  # replicates would make this test fail with an unexpected condition.
  expect_message(
    expect_warning(out <- boot_q(lm_cr(), pipe, n_boot = 3), "careful"),
    "hello"
  )
  expect_equal(out$n_boot_ok, rep(3L, nrow(out)))
})

test_that("replicates run with the progress option turned off", {
  rlang::local_options(tidylandmark.progress = TRUE)
  seen <- logical(0)
  pipe <- function(d) {
    seen <<- c(seen, getOption("tidylandmark.progress"))
    cox_pipe(d)
  }
  set.seed(11)
  boot_q(lm_cr(), pipe, n_boot = 3)
  expect_equal(seen, c(TRUE, FALSE, FALSE, FALSE))
  # The option is restored afterwards
  expect_true(getOption("tidylandmark.progress"))
})

# bootstrap_p() --------------------------------------------------------------------

test_that("bootstrap_p() computes percentile p-values", {
  est <- c(-1, 1, 2, 3)
  expect_equal(bootstrap_p(est, 0, "greater"), 0.25)
  expect_equal(bootstrap_p(est, 0, "less"), 0.75)
  expect_equal(bootstrap_p(est, 0, "two.sided"), 0.5)
  expect_equal(bootstrap_p(c(1, 2, 3), 0, "two.sided"), 0)
  # Capped at 1
  expect_equal(bootstrap_p(c(-1, 1), 0, "two.sided"), 1)
  # NA estimates (failed replicates) are dropped
  expect_equal(bootstrap_p(c(NA, -1, 1, 2, 3), 0, "greater"), 0.25)
  expect_identical(bootstrap_p(c(NA_real_, NA_real_), 0, "two.sided"), NA_real_)
  expect_identical(bootstrap_p(numeric(0), 0, "two.sided"), NA_real_)
})

# Mixed models with resampled subjects ----------------------------------------------

lme_q <- function(lm, sim) {
  suppressWarnings(suppressMessages(
    add_landmark_covariates(lm, sim$visits, id = pid, time = day, lme = marker, progress = FALSE)
  ))
}
lme_model_at <- function(x, s) {
  fits <- landmark_lme_fits(x)
  fits$model[[which(fits$.landmark == s)]]
}

test_that("a subject drawn twice gets a separate random effect per copy", {
  sim <- sim_lme()
  s <- 100
  lm <- sim$lm
  dup_row <- which(lm$.landmark == s & lm$.id == 5)
  lm_dup <- lm[c(seq_len(nrow(lm)), dup_row), ]

  m_orig <- lme_model_at(lme_q(lm, sim), s)
  out_dup <- lme_q(lm_dup, sim)
  m_dup <- lme_model_at(out_dup, s)

  # One extra group, carrying a second copy of subject 5's history
  n5 <- sum(sim$visits$pid == 5 & sim$visits$day <= s & !is.na(sim$visits$marker))
  expect_equal(lme4::ngrps(m_dup)[[".id"]], lme4::ngrps(m_orig)[[".id"]] + 1)
  expect_equal(stats::nobs(m_dup), stats::nobs(m_orig) + n5)
  expect_true(all(c("5", "5#1") %in% rownames(lme4::ranef(m_dup)$.id)))

  # Identical copies get identical predictions, one per row
  pred5 <- out_dup$marker[out_dup$.landmark == s & out_dup$.id == 5]
  expect_length(pred5, 2L)
  expect_equal(pred5[1], pred5[2])
  expect_equal(sum(is.na(out_dup$marker)), sum(is.na(lme_q(lm, sim)$marker)))
})

test_that("without duplicates, random effects are still labelled by subject id", {
  sim <- sim_lme()
  m <- lme_model_at(lme_q(sim$lm, sim), 100)
  expect_false(any(grepl("#", rownames(lme4::ranef(m)$.id), fixed = TRUE)))
  expect_true(all(rownames(lme4::ranef(m)$.id) %in% as.character(sim$subjects$pid)))
})

test_that("variance components are reported for bootstrap-style duplicated data", {
  sim <- sim_lme()
  lm <- sim$lm
  lm_dup <- lm[c(seq_len(nrow(lm)), which(lm$.id == 5)), ]
  fits <- landmark_lme_fits(lme_q(lm_dup, sim))
  ok <- fits$status != "failed"
  expect_false(anyNA(fits$var_intercept[ok]))
})

test_that("bootstrap_landmark() runs an LME pipeline end to end", {
  sim <- sim_lme()
  pipe <- function(d) {
    d |>
      add_landmark_covariates(sim$visits, id = pid, time = day, lme = marker, progress = FALSE) |>
      fit_q(~ age + marker)
  }
  set.seed(12)
  out <- suppressWarnings(suppressMessages(boot_q(sim$lm, pipe, n_boot = 3)))

  expect_true(any(grepl("lme", out$origin, fixed = TRUE)))
  expect_true(all(out$n_boot_ok[!is.na(out$term)] > 0L))
  expect_false(anyNA(out$estimate__boot[out$term == "marker"]))
})

test_that("factor coefficients keep their full-data contrasts", {
  d <- sim_cr(100)
  d$grp <- factor(c("a", rep(c("b", "c"), length.out = 99)))
  x <- lm_cr(d, landmarks = 0, horizon = NULL)
  # Both draws omit the only reference subject. coxph still returns a finite
  # grpb coefficient under its singular-fit convention.
  for (model in c("cox", "finegray")) {
    set.seed(31)
    pipe <- function(z) fit_q(z, ~ grp + age, model = model)
    result <- boot_q(x, pipe, n_boot = 2)
    expect_equal(result$n_boot_ok, c(0L, 0L, 2L))
    expect_equal(result$.fit_status__boot[[1]], rep("unestimable", 2))
    expect_true(all(is.na(result$.estimates__boot[[1]])))
  }

  # Losing a non-reference level leaves b versus a estimable, including
  # when the pipeline drops unused levels and the factor coding changes.
  d$grp <- factor(c("c", rep(c("a", "b"), length.out = 99)))
  x <- lm_cr(d, landmarks = 0, horizon = NULL)
  set.seed(31)
  result <- boot_q(x, function(z) {
    z$grp <- droplevels(z$grp)
    fit_q(z, ~ grp + age)
  }, n_boot = 2)
  expect_equal(result$n_boot_ok, c(2L, 0L, 2L))
  set.seed(31)
  for (b in 1:2) {
    z <- x[sample.int(nrow(x), replace = TRUE), ]
    direct <- ref_cox(z, 0, "I(grp == 'b') + age")
    expect_equal(result$.estimates__boot[[1]][b], unname(stats::coef(direct)[1]))
  }
})

test_that("coefficient intervals require a common spline basis", {
  x <- lm_cr(sim_cr(100), landmarks = 0, horizon = NULL)
  set.seed(13)
  expect_error(boot_q(x, \(z) fit_q(z, ~ splines::ns(age, 3)), n_boot = 3),
               "fix spline knots")
  set.seed(13)
  fixed <- boot_q(x, \(z) fit_q(z, ~ splines::ns(age, knots = c(44, 57),
                                                Boundary.knots = c(30, 70))), n_boot = 3)
  expect_equal(fixed$n_boot_ok, rep(3L, 3))
})

test_that("precomputed LME covariates are checked only when used", {
  sim <- sim_lme()
  x <- lme_q(sim$lm, sim)
  expect_error(boot_q(x, \(z) fit_q(z, ~ marker), n_boot = 2),
               "precomputed LME covariate.*marker")
  expect_no_error(boot_q(x, \(z) fit_q(z, ~ age), n_boot = 2))
})

test_that("failed draws keep their positions at all full-data landmarks", {
  x <- suppressWarnings(lm_cr(landmarks = c(50, 2000)))
  calls <- 0L
  pipe <- function(z) {
    calls <<- calls + 1L
    if (calls %in% c(2L, 4L)) stop("failed draw")
    fit_q(z, ~ age)
  }
  out <- suppressWarnings(boot_q(x, pipe, n_boot = 3))
  expect_equal(out$n_boot_ok, c(1L, 0L))
  expect_equal(is.na(out$.estimates__boot[[1]]), c(TRUE, FALSE, TRUE))
  expect_equal(out$.estimates__boot[[2]], rep(NA_real_, 3))
  expect_equal(out$.fit_status__boot[[2]], c("failed", "unestimable", "failed"))
  expect_equal(out$.message__boot[[2]][c(1, 3)], rep("failed draw", 2))
  expect_equal(out$estimate__boot[1], out$.estimates__boot[[1]][2])
  expect_true(all(is.na(out$conf.low__boot)))
  expect_true(all(is.na(out$conf.high__boot)))
  expect_true(all(is.na(out$p.value__boot)))
})

test_that("LME adjustment variables can share the helper's preferred name", {
  sim <- sim_lme()
  x <- sim$lm[c(seq_len(nrow(sim$lm)), which(sim$lm$.id == 5)), ]
  x$.re_key <- x$age
  derive <- function(adjust) suppressWarnings(suppressMessages(
    add_landmark_covariates(x, sim$visits, pid, day, lme = marker,
                            lme_args = list(adjust = adjust), progress = FALSE)
  ))
  expect_equal(derive(~ .re_key)$marker, derive(~ age)$marker)
})
