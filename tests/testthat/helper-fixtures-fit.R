# Fixtures for fit_landmark() tests. testthat sources every helper-*.R file
# before running tests, so these are available alongside helper-fixtures.R.

# Competing-risks data with continuous times (no ties), two event types and
# baseline covariates with real effects on cause 1.
#   status: 0 = censored, 1 = target event, 2 = competing event
sim_cr <- function(n = 400, seed = 123) {
  set.seed(seed)
  age <- round(stats::runif(n, 30, 70))
  sex <- factor(sample(c("F", "M"), n, replace = TRUE))
  grp <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  t1   <- stats::rexp(n, 0.0015 * exp(0.03 * (age - 50) + 0.4 * (sex == "M")))
  t2   <- stats::rexp(n, 0.0008 * exp(-0.02 * (age - 50)))
  cens <- stats::runif(n, 300, 1500)
  time <- pmin(t1, t2, cens)
  tibble::tibble(
    pid    = seq_len(n),
    time   = time,
    status = ifelse(time == t1, 1L, ifelse(time == t2, 2L, 0L)),
    age    = age,
    sex    = sex,
    grp    = grp
  )
}

lm_cr <- function(d = sim_cr(), landmarks = c(50, 200, 400), horizon = 730) {
  landmark_data(d, pid, time, status, landmarks = landmarks, horizon = horizon)
}

# Rows of one landmark as a plain data frame.
at_landmark <- function(lm, s) {
  as.data.frame(tibble::as_tibble(lm)[lm$.landmark == s, , drop = FALSE])
}

# Reference fits written directly against survival, independent of the
# package code. `rhs` is the covariate part of the formula as text.
# survival isn't attached in tests, so Surv() and strata() are bound locally;
# as.formula() uses this function's environment, where they are found.
ref_cox <- function(lm, s, rhs, cause = 1) {
  Surv <- survival::Surv
  strata <- survival::strata
  d <- at_landmark(lm, s)
  f <- stats::as.formula(paste0("Surv(.time, .status == ", cause, ") ~ ", rhs))
  survival::coxph(f, data = d)
}

ref_fg <- function(lm, s, rhs, cause = 1) {
  Surv <- survival::Surv
  strata <- survival::strata
  d <- at_landmark(lm, s)
  d$ev  <- factor(d$.status, levels = 0:2)
  d$cid <- seq_len(nrow(d))
  fg <- survival::finegray(Surv(.time, ev) ~ ., data = d, etype = as.character(cause))
  f <- stats::as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", rhs))
  survival::coxph(f, data = fg, weights = fgwt, cluster = cid)
}
