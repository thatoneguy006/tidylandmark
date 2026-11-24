# `sim_cr()`, `lm_cr()`, `ref_cox()` and `ref_fg()` are in the helper file.
# Landmarks 50, 200, 400 with a 730-day horizon; status 0/1/2.

fit_q <- function(...) fit_landmark(..., progress = FALSE)
coefs_at <- function(fit, s) fit$coefficients[fit$coefficients$.landmark == s, ]

# Agreement with direct survival fits -------------------------------------------

test_that("Cox fits match coxph() on each landmark's risk set", {
  lm <- lm_cr()
  expect_no_warning(fit <- fit_q(lm, ~ age + sex))

  for (s in c(50, 200, 400)) {
    ref <- ref_cox(lm, s, "age + sex")
    m <- fit$models[[as.character(s)]]
    expect_equal(coef(m), coef(ref))
    expect_equal(m$var, ref$var)
    co <- coefs_at(fit, s)
    expect_equal(co$term, c("age", "sexM"))
    expect_equal(co$estimate, unname(coef(ref)))
    expect_equal(co$std.error, unname(sqrt(diag(ref$var))))
  }

  # Stored fits work with survival's diagnostics outside fit_landmark()
  expect_no_error(survival::cox.zph(fit$models[["50"]]))
  expect_no_error(stats::residuals(fit$models[["50"]], type = "schoenfeld"))
})

test_that("Fine-Gray fits match finegray() + weighted, clustered coxph()", {
  lm <- lm_cr()
  expect_no_warning(fit <- fit_q(lm, ~ age + sex, model = "finegray"))

  for (s in c(50, 200, 400)) {
    ref <- ref_fg(lm, s, "age + sex")
    m <- fit$models[[as.character(s)]]
    expect_equal(coef(m), coef(ref))
    expect_equal(m$var, ref$var)               # robust variance
    expect_equal(m$naive.var, ref$naive.var)
    expect_equal(coefs_at(fit, s)$std.error, unname(sqrt(diag(ref$var))))
  }
})

test_that("Fine-Gray estimates agree with cmprsk::crr()", {
  skip_if_not_installed("cmprsk")
  lm <- lm_cr()
  fit <- fit_q(lm, ~ age + sex, model = "finegray")

  d <- at_landmark(lm, 50)
  ref <- cmprsk::crr(
    d$.time, d$.status,
    cov1 = cbind(age = d$age, sexM = as.numeric(d$sex == "M")),
    failcode = 1, cencode = 0
  )
  co <- coefs_at(fit, 50)
  expect_equal(co$estimate, unname(ref$coef), tolerance = 1e-3)
  # The robust coxph SE and crr's SE differ slightly (crr includes a term for
  # estimating the censoring distribution), so only rough agreement.
  expect_equal(co$std.error, unname(sqrt(diag(ref$var))), tolerance = 0.05)
})

test_that("formula features work as in coxph()", {
  lm <- lm_cr()
  for (rhs in c("splines::ns(age, df = 2) + sex * grp", "age + strata(sex)")) {
    fit <- fit_q(lm, stats::as.formula(paste("~", rhs)))
    expect_equal(coef(fit$models[["200"]]), coef(ref_cox(lm, 200, rhs)))
  }

  # User-defined functions and variables resolve in the formula's environment
  sq <- function(x) x^2
  fit_sq <- fit_q(lm, ~ age + sq(age))
  expect_equal(
    coefs_at(fit_sq, 200)$estimate,
    unname(coef(ref_cox(lm, 200, "age + I(age^2)")))
  )
  k <- 3
  fit_k <- fit_q(lm, ~ splines::ns(age, df = k))
  expect_equal(nrow(coefs_at(fit_k, 200)), 3L)
})

# Models, counts and missing data --------------------------------------------------

test_that("the model label reflects the status coding", {
  expect_equal(fit_q(lm_cr(), ~ age)$model_info$.model, rep("cause-specific cox", 3))
  expect_equal(
    fit_q(lm_cr(), ~ age, model = "finegray")$model_info$.model,
    rep("finegray", 3)
  )

  d <- sim_cr()
  d$status[d$status == 2L] <- 0L
  expect_equal(fit_q(lm_cr(d), ~ age)$model_info$.model, rep("cox", 3))
})

test_that("counts describe each landmark's model sample", {
  d <- sim_cr()
  d$age[c(5, 10, 15, 20, 25)] <- NA
  lm <- lm_cr(d)
  fit <- fit_q(lm, ~ age + sex)
  info <- fit$model_info

  rows <- tibble::as_tibble(lm)
  expect_equal(info$n_eligible, as.vector(table(rows$.landmark)))
  expect_equal(info$n_eligible - info$n_model, as.vector(tapply(is.na(rows$age), rows$.landmark, sum)))

  cc <- rows[!is.na(rows$age), ]
  expect_equal(info$n_event, as.vector(tapply(cc$.status == 1L, cc$.landmark, sum)))
  expect_equal(info$n_competing, as.vector(tapply(cc$.status == 2L, cc$.landmark, sum)))

  # Complete-case fit equals coxph()'s default na.omit fit
  expect_equal(coef(fit$models[["50"]]), coef(ref_cox(lm, 50, "age + sex")))
})

test_that("a factor status accepts a level name as cause", {
  d <- sim_cr()
  d$status <- factor(d$status, levels = 0:2, labels = c("censored", "relapse", "death"))
  lm <- lm_cr(d)

  by_label <- fit_q(lm, ~ age, cause = "death")
  by_code  <- fit_q(lm, ~ age, cause = 2)
  expect_equal(by_label$coefficients, by_code$coefficients)
  expect_equal(by_label$model_info$.cause, rep(2L, 3))
  expect_output(print(by_label), "cause 2 \\(death\\)")

  expect_error(fit_q(lm, ~ age, cause = "censored"), "relapse")
})

test_that("fitting a filtered landmark_tbl fits only the remaining landmarks", {
  lm <- dplyr::filter(lm_cr(), .landmark == 200)
  fit <- fit_q(lm, ~ age)
  expect_equal(fit$model_info$.landmark, 200)
  expect_named(fit$models, "200")
})

# Unestimable and failed landmarks ---------------------------------------------

test_that("landmarks without target or competing events are handled", {
  d <- sim_cr()
  d$status[d$status == 2L & d$time > 400] <- 0L  # no competing events after day 400
  lm <- lm_cr(d, landmarks = c(50, 400))

  expect_warning(
    fit_fg <- fit_q(lm, ~ age, model = "finegray"),
    "1 of 2 landmark models did not fit cleanly"
  )
  info <- fit_fg$model_info
  expect_equal(info$.fit_status, c("ok", "warning"))
  expect_equal(info$n_competing[2], 0L)
  expect_match(info$.message[2], "no competing events")
  # With no competing events, Fine-Gray equals the Cox model
  expect_equal(
    coef(fit_fg$models[["400"]]),
    coef(fit_q(lm, ~ age)$models[["400"]]),
    tolerance = 1e-6
  )

  # No events of cause 2 after day 400: nothing to estimate
  expect_warning(fit2 <- fit_q(lm, ~ age, cause = 2), "did not fit cleanly")
  expect_equal(fit2$model_info$.fit_status, c("ok", "unestimable"))
  expect_match(fit2$model_info$.message[2], "no target events \\(cause 2\\)")
  expect_null(fit2$models[["400"]])
  expect_equal(nrow(coefs_at(fit2, 400)), 0L)
})

test_that("a covariate with one observed value makes a landmark unestimable", {
  d <- sim_cr()
  d$time[d$sex == "M"] <- pmin(d$time[d$sex == "M"], 100)  # no men at risk at 200
  lm <- lm_cr(d, landmarks = c(50, 200))

  expect_warning(fit <- fit_q(lm, ~ age + sex), "did not fit cleanly")
  expect_equal(fit$model_info$.fit_status, c("ok", "unestimable"))
  expect_match(fit$model_info$.message[2], "sex has a single observed value")

  # Inside strata() a single level is fine
  expect_no_warning(fit_s <- fit_q(lm, ~ age + strata(sex)))
  expect_equal(fit_s$model_info$.fit_status, c("ok", "ok"))
})

test_that("a landmark with nobody at risk is reported, not dropped", {
  lm <- suppressWarnings(lm_cr(landmarks = c(50, 1e6)))
  expect_warning(fit <- fit_q(lm, ~ age), "did not fit cleanly")
  expect_equal(fit$model_info$.landmark, c(50, 1e6))
  expect_equal(fit$model_info$n_eligible[2], 0L)
  expect_equal(fit$model_info$.fit_status[2], "unestimable")
  expect_match(fit$model_info$.message[2], "no subjects at risk")
})

test_that("an error at one landmark doesn't stop the others", {
  lm <- lm_cr()
  sizes <- as.vector(table(tibble::as_tibble(lm)$.landmark))
  threshold <- mean(range(sizes))
  boom <- function(x) {
    if (length(x) < threshold) stop("boom")
    x
  }

  expect_warning(fit <- fit_q(lm, ~ boom(age)), "did not fit cleanly")
  info <- fit$model_info
  expect_equal(info$.fit_status == "failed", sizes < threshold)
  expect_true(all(grepl("boom", info$.message[info$.fit_status == "failed"])))
  expect_true(all(info$.fit_status[sizes >= threshold] == "ok"))
})

test_that("problems with LME covariates are carried into the fit status", {
  sim <- sim_lme()
  lm <- suppressWarnings(add_landmark_covariates(
    sim$lm, sim$visits, pid, day,
    lme = marker, lme_args = list(min_subjects = 1000), progress = FALSE
  ))
  expect_warning(fit <- fit_q(lm, ~ age + marker), "did not fit cleanly")
  expect_equal(fit$model_info$.fit_status, c("unestimable", "unestimable"))
  expect_match(fit$model_info$.message, "no complete cases")
  expect_match(fit$model_info$.message, "lme covariate marker \\(failed\\)")
})

# Input validation -------------------------------------------------------------------

test_that("fit_landmark() validates its inputs", {
  lm <- lm_cr()
  expect_error(fit_q(tibble::as_tibble(lm), ~ age), "landmark_tbl")
  expect_error(fit_q(lm), "formula")
  expect_error(fit_q(lm, "~ age"), "one-sided formula")
  expect_error(fit_q(lm, Surv(.time, .status) ~ age), "one-sided")
  expect_error(fit_q(lm, ~ .), "explicitly")
  expect_error(fit_q(lm, ~ age + cluster(.id)), "cluster")
  expect_error(fit_q(lm, ~ age + tt(age)), "tt")
  expect_error(fit_q(lm, ~ age + .time), "reserved")
  expect_error(fit_q(lm, ~ age + not_a_column), "not_a_column")
  expect_error(fit_q(lm, ~ 1), "at least one column")
  expect_error(fit_q(lm, ~ age, model = "aft"), "model")
  expect_error(fit_q(lm, ~ age, cause = 3), "event codes")
  expect_error(fit_q(lm, ~ age, cause = c(1, 2)), "event codes")
  expect_error(fit_q(lm, ~ age, cause = "death"), "factor")
  expect_error(fit_landmark(lm, ~ age, progress = NA), "progress")

  # finegray() column names can't be used as covariates in Fine-Gray models
  d <- sim_cr()
  d$fgwt <- stats::rnorm(nrow(d))
  lm_fg <- lm_cr(d)
  expect_error(fit_q(lm_fg, ~ age + fgwt, model = "finegray"), "fgwt")
  expect_no_error(fit_q(lm_fg, ~ age + fgwt))
})

# tidy(), print() and summary() -------------------------------------------------------

test_that("tidy() matches coxph() inference and has a stable schema", {
  lm <- lm_cr()
  fit <- fit_q(lm, ~ age + sex)
  td <- tidy(fit, conf.int = TRUE)

  expect_named(td, c(
    ".landmark", "term", "estimate", "std.error", "statistic", "p.value",
    "conf.low", "conf.high", "origin", "n_eligible", "n_model", "n_event",
    "n_competing", ".model", ".cause", ".fit_status", ".message"
  ))
  expect_named(tidy(fit), setdiff(names(td), c("conf.low", "conf.high")))
  expect_equal(nrow(td), 6L)

  m  <- fit$models[["200"]]
  sm <- summary(m)$coefficients
  row <- td[td$.landmark == 200, ]
  expect_equal(row$statistic, unname(sm[, "z"]))
  expect_equal(row$p.value, unname(sm[, "Pr(>|z|)"]))
  expect_equal(row$conf.low, unname(stats::confint(m)[, 1]))
  expect_equal(row$conf.high, unname(stats::confint(m)[, 2]))

  td90 <- tidy(fit, conf.int = TRUE, conf.level = 0.90)
  expect_equal(td90$conf.low[td90$.landmark == 200], unname(stats::confint(m, level = 0.90)[, 1]))

  tde <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  expect_equal(tde$estimate, exp(td$estimate))
  expect_equal(tde$conf.low, exp(td$conf.low))
  expect_equal(tde$std.error, td$std.error)

  expect_equal(unique(td$origin), "baseline")
  expect_error(tidy(fit, conf.level = 95), "conf.level")
  expect_error(tidy(fit, exponentiate = "yes"), "exponentiate")
})

test_that("tidy() uses the robust standard error for Fine-Gray models", {
  fit <- fit_q(lm_cr(), ~ age + sex, model = "finegray")
  sm  <- summary(fit$models[["50"]])$coefficients
  td  <- tidy(fit)
  expect_equal(td$std.error[td$.landmark == 50], unname(sm[, "robust se"]))
  expect_equal(td$p.value[td$.landmark == 50], unname(sm[, "Pr(>|z|)"]))
})

test_that("tidy() keeps a row for landmarks without estimates", {
  d <- sim_cr()
  d$status[d$status == 2L & d$time > 400] <- 0L
  fit <- suppressWarnings(fit_q(lm_cr(d, landmarks = c(50, 400)), ~ age, cause = 2))
  td <- tidy(fit)
  at400 <- td[td$.landmark == 400, ]
  expect_equal(nrow(at400), 1L)
  expect_true(is.na(at400$term))
  expect_true(is.na(at400$estimate))
  expect_equal(at400$.fit_status, "unestimable")
})

test_that("tidy() records where each term came from", {
  sim <- sim_lme()
  lm <- suppressWarnings(
    sim$lm |>
      add_landmark_covariates(
        sim$visits, pid, day, locf = marker, names = "{.col}_locf", progress = FALSE
      ) |>
      add_landmark_covariates(sim$visits, pid, day, lme = marker, progress = FALSE)
  )
  fit <- suppressWarnings(fit_q(lm, ~ age * grp + marker_locf + marker + marker:age))

  rlang::local_options(rlib_message_verbosity = "verbose")
  expect_message(tidy(fit), "LME-derived covariates")
  td <- suppressMessages(tidy(fit))
  origin <- rlang::set_names(td$origin, td$term)

  expect_equal(unname(origin[td$.landmark == 100][c("age", "grpb", "age:grpb")]), rep("baseline", 3))
  expect_equal(unname(origin[td$.landmark == 100]["marker_locf"]), "locf")
  expect_equal(unname(origin[td$.landmark == 100]["marker"]), "lme")
  expect_true("baseline:lme" %in% td$origin)
  expect_equal(fit$config$lme_columns, "marker")

  # No message without LME covariates
  fit_locf <- suppressWarnings(fit_q(lm, ~ age + marker_locf))
  expect_no_message(tidy(fit_locf))
})

test_that("print() and summary() describe the fit", {
  fit <- fit_q(lm_cr(), ~ age + sex, model = "finegray")
  expect_output(print(fit), "<landmark_fit> Fine-Gray model for cause 1 at 3 landmark\\(s\\)")
  expect_output(print(fit), "Formula: ~age \\+ sex")
  expect_output(print(fit), "Horizon: 730")
  expect_output(print(fit), "n_competing")
  expect_output(expect_invisible(print(fit)))

  cox <- fit_q(lm_cr(), ~ age)
  expect_output(print(cox), "Cause-specific Cox")

  s <- summary(fit)
  expect_equal(s, tidy(fit, exponentiate = TRUE, conf.int = TRUE))
})
