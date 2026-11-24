#' Attach pre-landmark covariates to a landmark dataset
#'
#' @description
#' `add_landmark_covariates()` adds one column per selected longitudinal
#' variable to a `landmark_tbl`. For every subject-landmark row, the value is
#' derived only from that subject's measurements taken at or before the
#' landmark, so no information from after the landmark leaks into the
#' covariates. Each variable is derived in one of two ways:
#'
#' * **`locf`**: last observation carried forward, optionally limited to a
#'   `lookback` window.
#' * **`lme`**: the subject-specific value at the landmark predicted from a
#'   linear mixed-effects model fitted separately for every variable and
#'   landmark.
#'
#' Rows are never added or removed, and `.time`/`.status` are untouched.
#' Subjects with no usable measurement keep their row with a missing value.
#'
#' @param landmark_data A `landmark_tbl`, as returned by [landmark_data()] (or
#'   by a previous call to `add_landmark_covariates()`).
#' @param longitudinal_data A data frame with **one row per measurement
#'   occasion** (multiple rows per subject). Rows whose `id` does not appear
#'   in `landmark_data` are ignored, and their number is reported in a
#'   message.
#' @param id <[`tidy-select`][tidyselect::language]> The subject identifier
#'   column in `longitudinal_data`. It must be non-missing and have the same
#'   type as the `id` used in [landmark_data()]; types are never coerced
#'   implicitly.
#' @param time <[`tidy-select`][tidyselect::language]> Numeric, finite,
#'   non-missing measurement time. It must be on the **same time scale and
#'   origin** as the `time` given to [landmark_data()]; this cannot be checked
#'   automatically. Measurements before time zero (e.g. screening visits) are
#'   allowed.
#' @param locf,lme <[`tidy-select`][tidyselect::language]> Columns of
#'   `longitudinal_data` to derive by last observation carried forward and by
#'   linear mixed-effects model, e.g. `locf = c(treatment, stage)` and
#'   `lme = c(bili, albumin)`. Either may be `NULL` (the default), but not
#'   both. Rename inside the selection to name a derived column, e.g.
#'   `lme = c(bili_lme = bili)`. `lme` variables must be continuous numeric.
#' @param include_at_landmark If `TRUE` (default), a measurement taken exactly
#'   at the landmark (\eqn{t = s}) may be used. If `FALSE`, only measurements
#'   strictly before the landmark (\eqn{t < s}) are used. Applies to both
#'   methods.
#' @param lookback `locf` only. `NULL` (default) carries values forward
#'   without a limit. A single non-negative number only uses measurements
#'   within that window: \eqn{s - t \le} `lookback`. Older values give `NA`.
#' @param lme_args `lme` only. A named list of model options:
#'   * `random_slope`: `FALSE` (default) fits a random intercept; `TRUE`
#'     adds a random slope for time.
#'   * `adjust`: `NULL` (default) or a one-sided formula of fixed-effect
#'     adjustment terms, e.g. `~ age + sex`. Variables come from
#'     `landmark_data`, so they are baseline values (or covariates added by an
#'     earlier call). Use `.t` for the time since the landmark, e.g.
#'     `~ sex * .t` for sex-specific slopes.
#'   * `REML`: passed to [lme4::lmer()]. Default `TRUE`.
#'   * `control`: an [lme4::lmerControl()] object. Default
#'     `lme4::lmerControl()`.
#'   * `min_subjects`: the fewest subjects with measurements needed to fit a
#'     model (default `5`). With `random_slope = TRUE`, this many subjects
#'     must have at least two measurements.
#' @param names `NULL` (default) or a template for derived column names, in
#'   which `"{.col}"` is replaced by the source variable and `"{.method}"` by
#'   `"locf"` or `"lme"`, e.g. `"{.col}_{.method}"`. Columns renamed in the
#'   selection keep their new name. A derived name that already exists in
#'   `landmark_data` is an error; existing columns are never overwritten.
#' @param progress Show a progress bar while fitting mixed models? Defaults to
#'   the `tidylandmark.progress` option, or [interactive()] if it is unset.
#'
#' @details
#' ## LOCF
#' For landmark \eqn{s}, subject \eqn{i} and variable \eqn{v}, the value is
#' the last *non-missing* measurement of \eqn{v} at time \eqn{t \le s} (or
#' \eqn{t < s} when `include_at_landmark = FALSE`), provided it is within
#' `lookback`. Each variable is handled separately, so a visit where one
#' variable is missing still contributes the others. Selection depends only
#' on `id` and `time`, never on row order. Any column type works, and the type
#' (including factor levels) is preserved. The time since each
#' carried-forward measurement is stored in the metadata for auditing.
#'
#' If a subject has two *different* non-missing values of an `locf` variable
#' at the same time, the function errors and lists the affected subjects and
#' times rather than choosing one arbitrarily. Exact duplicates are collapsed.
#'
#' ## LME
#' For every `lme` variable \eqn{v} and landmark \eqn{s}, an independent
#' Gaussian model is fitted with [lme4::lmer()] to the non-missing
#' measurements, taken at or before \eqn{s}, of the subjects in the risk set at
#' \eqn{s}:
#'
#' \deqn{Y_{ij} = \beta_0 + \beta_1 (t_{ij} - s) + \gamma^\top x_i + b_{0i} [+ b_{1i}(t_{ij} - s)] + \epsilon_{ij}}
#'
#' where \eqn{x_i} are the optional `adjust` terms. Time is centred at the
#' landmark and rescaled internally for numerical stability, which does not
#' change the predictions. The derived value is the conditional prediction at
#' \eqn{t = s} including the subject's random effects
#' (`predict(fit, re.form = NULL)`), i.e. \eqn{\hat\beta_0 + \hat\gamma^\top
#' x_i + \hat b_{0i}}. A subject with a single measurement is kept; their
#' prediction is shrunk towards the population trajectory. A subject with no
#' measurement before the landmark, or with missing `adjust` values, gets
#' `NA` rather than a population-only prediction. Different values recorded
#' at the same time are treated as replicate measurements.
#'
#' Each fit gets a status: `"ok"`; `"warning"` when [lme4::lmer()] warned
#' (e.g. about convergence) or the fit is singular ([lme4::isSingular()]); or
#' `"failed"` when there were too few subjects or [lme4::lmer()] gave an
#' error. Failed fits give `NA` for that variable at that landmark without
#' affecting other fits, and a warning lists every fit that was not `"ok"`.
#' There is no automatic fallback from a random slope to a random intercept.
#'
#' Standard errors in a second-stage survival model treat LME-derived values
#' as known, so they understate uncertainty.
#'
#' @returns `landmark_data` with one new column per selected variable (`locf`
#'   columns first, then `lme`), appended after the existing columns. The
#'   class and configuration are kept, and the configuration records how each
#'   derived column was made: method, source variable, options, the LOCF audit
#'   trail, and for LME the fitted models and their diagnostics. The printed
#'   header lists the derived columns and their method.
#'
#' @seealso [landmark_data()] to build the risk sets first.
#' @export
#'
#' @examples
#' subjects <- data.frame(
#'   patient_id    = 1:10,
#'   followup_days = c(20, 30, 48, 65, 85, 110, 145, 190, 250, 270),
#'   event         = c(1, 2, 1, 2, 0, 1, 2, 1, 0, 1)
#' )
#' visits <- data.frame(
#'   patient_id = c(6L, 6L, 6L, 7L, 7L, 7L, 10L, 10L),
#'   visit_day  = c(0, 45, 90, 0, 50, 120, 60, 140),
#'   treatment  = c(0, 1, 1, 0, 0, 1, 1, 1),
#'   biomarker  = c(25, 20, 17, 23, NA, 18, 24, 20)
#' )
#'
#' lm_dat <- subjects |>
#'   landmark_data(
#'     patient_id, followup_days, event,
#'     landmarks = c(30, 60, 90), horizon = 180
#'   ) |>
#'   add_landmark_covariates(
#'     visits, id = patient_id, time = visit_day,
#'     locf = c(treatment, biomarker)
#'   )
#' lm_dat
#'
#' # Only carry values forward for up to 30 days, under new names
#' landmark_data(subjects, patient_id, followup_days, event, landmarks = 60) |>
#'   add_landmark_covariates(
#'     visits, patient_id, visit_day,
#'     locf = biomarker, lookback = 30, names = "{.col}_recent"
#'   )
#'
#' # Mixed models with the PBC data: log bilirubin by LME, adjusted
#' # for age, alongside LOCF for albumin
#' visits_pbc <- transform(pbcseq_visits, log_bili = log(bili))
#' pbc_lm <- pbcseq_baseline |>
#'   landmark_data(id, futime, status, landmarks = c(365, 730), horizon = 1825) |>
#'   add_landmark_covariates(
#'     visits_pbc, id, day,
#'     locf = albumin,
#'     lme = log_bili,
#'     lme_args = list(adjust = ~ age)
#'   )
#' pbc_lm
add_landmark_covariates <- function(landmark_data,
                                    longitudinal_data,
                                    id,
                                    time,
                                    locf = NULL,
                                    lme = NULL,
                                    include_at_landmark = TRUE,
                                    lookback = NULL,
                                    lme_args = list(),
                                    names = NULL,
                                    progress = getOption("tidylandmark.progress", interactive())) {
  call <- rlang::current_env()

  # ---- 1. Validate containers and scalar arguments -------------------------
  lm_check_intact(landmark_data, "landmark_data", call)
  if (!is.data.frame(longitudinal_data)) {
    cli::cli_abort(
      "{.arg longitudinal_data} must be a data frame, not {.cls {class(longitudinal_data)}}.",
      call = call
    )
  }
  rlang::check_required(id)
  rlang::check_required(time)
  lm_check_flag(include_at_landmark, "include_at_landmark", call)
  lm_check_flag(progress, "progress", call)
  lm_check_names_template(names, call = call)

  locf_quo <- rlang::enquo(locf)
  lme_quo  <- rlang::enquo(lme)
  if (rlang::quo_is_null(locf_quo) && rlang::quo_is_null(lme_quo)) {
    cli::cli_abort(
      c(
        "Select variables with {.arg locf}, {.arg lme}, or both.",
        "i" = "For example {.code locf = c(treatment, stage), lme = bili}."
      ),
      call = call
    )
  }

  long <- tibble::as_tibble(longitudinal_data)

  # ---- 2. Resolve columns and output names ---------------------------------
  id_col   <- lm_select_one(long, rlang::enquo(id), "id", call)
  time_col <- lm_select_one(long, rlang::enquo(time), "time", call)
  if (identical(id_col, time_col)) {
    cli::cli_abort(
      "{.arg id} and {.arg time} must refer to two different columns.",
      call = call
    )
  }
  role_cols <- c(id_col, time_col)
  locf_sel <- lm_select_vars(long, locf_quo, role_cols, "locf", call = call)
  lme_sel  <- lm_select_vars(long, lme_quo, role_cols, "lme", call = call)

  locf_names <- lm_output_names(locf_sel, "locf", names, call = call)
  lme_names  <- lm_output_names(lme_sel, "lme", names, call = call)
  out_names  <- c(locf_names, lme_names)
  lm_check_unique_names(out_names, call)
  lm_check_collisions(out_names, base::names(landmark_data), call)

  lookback <- lm_check_lookback(lookback, length(locf_sel) > 0L, call = call)
  lme_opts <- lm_check_lme_args(
    lme_args, length(lme_sel) > 0L, base::names(landmark_data), call = call
  )

  # ---- 3. Validate column contents -----------------------------------------
  lm_check_no_na(long[[id_col]], "id", call)
  lm_check_no_na(long[[time_col]], "time", call)
  lm_check_measure_time(long[[time_col]], "time", call)
  lm_check_id_type(landmark_data$.id, long[[id_col]], call)
  purrr::walk(unname(lme_sel), \(v) lm_check_lme_var(long[[v]], v, call))

  # ---- 4. Keep only measurements of subjects in the landmark data -----------
  in_lm <- long[[id_col]] %in% landmark_data$.id
  n_unmatched <- sum(!in_lm)
  if (n_unmatched > 0L) {
    cli::cli_inform(
      c("i" = "Ignored {n_unmatched} row{?s} of {.arg longitudinal_data} whose {.arg id} is not in {.arg landmark_data}.")
    )
  }
  long <- long[in_lm, , drop = FALSE]
  ids   <- long[[id_col]]
  times <- as.double(long[[time_col]])

  # ---- 5. Derive each variable independently --------------------------------
  keys <- tibble::tibble(.id = landmark_data$.id, .landmark = landmark_data$.landmark)

  locf_out <- purrr::map(unname(locf_sel), \(v) {
    lm_locf_one(
      keys, id = ids, time = times, value = long[[v]], var = v,
      include_at_landmark = include_at_landmark, lookback = lookback,
      call = call
    )
  })

  lme_out <- list()
  if (length(lme_sel) > 0L) {
    adjust_vars <- setdiff(all.vars(lme_opts$adjust %||% ~1), ".t")
    lm_rows <- tibble::as_tibble(landmark_data)[c(".id", ".landmark", adjust_vars)]
    lme_out <- lm_lme_all(
      lm_rows, ids = ids, times = times,
      values = purrr::map(unname(lme_sel), \(v) long[[v]]),
      include_at_landmark = include_at_landmark,
      opts = lme_opts,
      progress = progress
    )
    lm_warn_lme_issues(lme_out, unname(lme_sel))
  }

  derived <- c(locf_out, lme_out)

  # ---- 6. Assemble output and record source ---------------------------------
  values <- rlang::set_names(purrr::map(derived, "value"), out_names)
  # mutate() on a landmark_tbl keeps its class and attributes.
  out <- dplyr::mutate(landmark_data, !!!values)

  n_locf <- length(locf_sel)
  n_lme  <- length(lme_sel)
  new_info <- tibble::tibble(
    column              = out_names,
    source              = c(unname(locf_sel), unname(lme_sel)),
    method              = rep(c("locf", "lme"), c(n_locf, n_lme)),
    include_at_landmark = include_at_landmark,
    lookback            = rep(c(lookback %||% NA_real_, NA_real_), c(n_locf, n_lme)),
    random_slope        = rep(c(NA, lme_opts$random_slope), c(n_locf, n_lme)),
    adjust              = rep(
      c(NA_character_, lm_deparse_adjust(lme_opts$adjust)),
      c(n_locf, n_lme)
    ),
    n_missing           = purrr::map_int(unname(values), \(x) sum(is.na(x))),
    details             = purrr::map(derived, "details")
  )

  cfg <- landmark_config(landmark_data)
  # Drop records for derived columns the user has since removed with select(),
  # so re-adding a column of the same name doesn't leave a stale record.
  old_info <- cfg$covariates
  if (!is.null(old_info)) {
    old_info <- old_info[old_info$column %in% base::names(landmark_data), , drop = FALSE]
  }
  cfg$covariates <- dplyr::bind_rows(old_info, new_info)
  attr(out, "landmark_config") <- cfg

  out
}

# LOCF ---------------------------------------------------------------------------

#' LOCF for one variable
#'
#' `keys` holds the (.id, .landmark) rows of the landmark data, in order.
#' `id`, `time` and `value` are parallel vectors from the longitudinal data.
#' Returns a list with `value` (one entry per row of `keys`, same type as the
#' input) and `details` (keys plus the time of the measurement used and its
#' lag behind the landmark).
#' @noRd
lm_locf_one <- function(keys, id, time, value, var,
                        include_at_landmark = TRUE,
                        lookback = NULL,
                        call = rlang::caller_env()) {
  # Built from vectors so internal column names can't clash with user columns.
  meas <- tibble::tibble(.id = id, .obs_time = as.double(time), .value = value)
  meas <- meas[!is.na(meas$.value), , drop = FALSE]

  lm_check_locf_conflicts(meas, var, call)
  meas <- dplyr::distinct(meas)

  # closest(.landmark >= .obs_time) picks, for each landmark row, the latest
  # measurement at or before the landmark. After the conflict check each
  # (.id, .obs_time) is unique, so at most one row matches.
  by <- if (include_at_landmark) {
    dplyr::join_by(".id", closest(".landmark" >= ".obs_time"))
  } else {
    dplyr::join_by(".id", closest(".landmark" > ".obs_time"))
  }
  matched <- dplyr::left_join(keys, meas, by = by, relationship = "many-to-one")

  # The match is the *latest* qualifying measurement, so if it is outside the
  # window every other candidate is too: blank it out.
  if (!is.null(lookback)) {
    too_old <- which(matched$.landmark - matched$.obs_time > lookback)
    matched$.value[too_old]    <- NA
    matched$.obs_time[too_old] <- NA_real_
  }

  list(
    value = matched$.value,
    details = tibble::tibble(
      .id       = matched$.id,
      .landmark = matched$.landmark,
      .obs_time = matched$.obs_time,
      .lag      = matched$.landmark - matched$.obs_time
    )
  )
}

# LME ----------------------------------------------------------------------------

#' Fit every (variable, landmark) model
#'
#' `lm_rows` holds .id, .landmark and any adjustment columns, one row per
#' landmark-data row, in order. `values` is a list of the lme variables'
#' vectors. Returns one list(value, details) per variable, like lm_locf_one().
#' @noRd
lm_lme_all <- function(lm_rows, ids, times, values, include_at_landmark,
                       opts, progress = FALSE) {
  landmarks <- sort(unique(lm_rows$.landmark))
  # Every (variable, landmark) pair, variable varying slowest.
  grid_var <- rep(seq_along(values), each = length(landmarks))
  grid_s   <- rep(landmarks, times = length(values))

  fits <- purrr::map2(
    grid_var, grid_s,
    \(k, s) {
      lm_lme_fit_one(
        lm_rows, ids = ids, times = times, value = values[[k]], s = s,
        include_at_landmark = include_at_landmark, opts = opts
      )
    },
    .progress = lm_progress_arg(progress, "Fitting mixed models")
  )

  purrr::map(seq_along(values), \(k) {
    mine <- fits[grid_var == k]
    value <- rep(NA_real_, nrow(lm_rows))
    for (f in mine) {
      value[f$rows] <- f$pred
    }
    list(
      value = value,
      details = purrr::list_rbind(purrr::map(mine, "diag"))
    )
  })
}

#' One mixed model: variable `value`, landmark `s`
#'
#' Returns list(rows, pred, diag): the landmark-data row numbers that got a
#' prediction, the predictions, and a one-row diagnostics tibble.
#' @noRd
lm_lme_fit_one <- function(lm_rows, ids, times, value, s,
                           include_at_landmark, opts) {
  landmark_rows <- which(lm_rows$.landmark == s)
  at_s <- lm_rows[landmark_rows, , drop = FALSE]
  adjust_vars <- setdiff(base::names(lm_rows), c(".id", ".landmark"))

  # Random-effect key, one per landmark row. Normally a subject has one row per
  # landmark, so the key is just its id and the model is unchanged. In a
  # bootstrap sample (bootstrap_landmark()) a subject drawn k times has k rows;
  # make.unique() gives them "12", "12#1", ... so each copy gets its own random
  # effect instead of being pooled into one subject with duplicated data.
  re_key <- make.unique(as.character(at_s$.id), sep = "#")
  has_dups <- anyDuplicated(at_s$.id) > 0L

  # History of subjects at risk at s, up to the landmark.
  in_window <- if (include_at_landmark) times <= s else times < s
  keep <- in_window & !is.na(value) & ids %in% at_s$.id
  obs <- tibble::tibble(.id = ids[keep], .time = times[keep], .value = value[keep])
  obs <- dplyr::distinct(obs)

  # Copy each subject's history once per landmark row, attaching that row's
  # random-effect key and adjustment values. Without duplicated subjects this
  # is many-to-one, which keeps the original safety check. Subjects missing an
  # adjustment value can't be predicted, so they are left out of the fit too.
  key_col <- make.unique(c(names(at_s), ".re_key"))[ncol(at_s) + 1L]
  at_s[[key_col]] <- re_key
  obs <- dplyr::left_join(
    obs, at_s[c(".id", key_col, adjust_vars)], by = ".id",
    relationship = if (has_dups) "many-to-many" else "many-to-one"
  )
  if (length(adjust_vars) > 0L) {
    obs <- obs[stats::complete.cases(obs[adjust_vars]), , drop = FALSE]
  }

  # Centre at the landmark and rescale. Predictions at .t = 0 don't depend on
  # the scale, but lmer converges more reliably on unit-scale time.
  t_centred <- obs$.time - s
  time_scale <- if (length(t_centred) > 1L) stats::sd(t_centred) else NA_real_
  if (!is.finite(time_scale) || time_scale <= 0) time_scale <- 1
  obs$.t <- t_centred / time_scale
  obs$.id <- obs[[key_col]]

  per_subject <- table(obs$.id)
  n_subjects  <- length(per_subject)
  n_multi     <- sum(per_subject >= 2L)

  diag <- tibble::tibble(
    .landmark     = s,
    status        = "ok",
    message       = NA_character_,
    n_eligible    = nrow(at_s),
    n_subjects    = n_subjects,
    n_obs         = nrow(obs),
    median_obs    = if (n_subjects > 0L) stats::median(as.numeric(per_subject)) else NA_real_,
    time_min      = if (nrow(obs) > 0L) min(obs$.time) else NA_real_,
    time_max      = if (nrow(obs) > 0L) max(obs$.time) else NA_real_,
    time_scale    = time_scale,
    slope         = NA_real_,
    var_intercept = NA_real_,
    var_slope     = NA_real_,
    var_residual  = NA_real_,
    singular      = NA,
    model         = list(NULL)
  )
  failed <- function(msg) {
    diag$status  <- "failed"
    diag$message <- msg
    list(rows = integer(0), pred = numeric(0), diag = diag)
  }

  if (n_subjects < opts$min_subjects) {
    return(failed(sprintf(
      "%d subject(s) with measurements; at least %d needed",
      n_subjects, opts$min_subjects
    )))
  }
  if (opts$random_slope && n_multi < opts$min_subjects) {
    return(failed(sprintf(
      "%d subject(s) with 2+ measurements; at least %d needed for a random slope",
      n_multi, opts$min_subjects
    )))
  }

  formula <- lm_lme_formula(opts$adjust, opts$random_slope)
  res <- lm_quiet_lmer(formula, obs, REML = opts$REML, control = opts$control)
  if (is.null(res$fit)) {
    return(failed(res$error))
  }
  fit <- res$fit

  fitted_rows <- re_key %in% obs$.id
  newdata <- at_s[fitted_rows, , drop = FALSE]
  newdata$.id <- re_key[fitted_rows]
  rows <- landmark_rows[fitted_rows]
  newdata$.t <- 0
  pred <- unname(stats::predict(fit, newdata = newdata, re.form = NULL))

  vc <- as.data.frame(lme4::VarCorr(fit))
  pick <- function(var1) {
    v <- vc$vcov[vc$grp == ".id" & vc$var1 %in% var1 & is.na(vc$var2)]
    if (length(v) == 1L) v else NA_real_
  }
  singular <- lme4::isSingular(fit)
  issues <- unique(c(res$warnings, if (singular) "singular fit"))

  diag$status        <- if (length(issues) > 0L) "warning" else "ok"
  diag$message       <- if (length(issues) > 0L) paste(issues, collapse = "; ") else NA_character_
  # Back on the original time scale.
  diag$slope         <- unname(lme4::fixef(fit)[".t"]) / time_scale
  diag$var_intercept <- pick("(Intercept)")
  diag$var_slope     <- pick(".t") / time_scale^2
  diag$var_residual  <- vc$vcov[vc$grp == "Residual"]
  diag$singular      <- singular
  diag$model         <- list(fit)

  list(rows = rows, pred = pred, diag = diag)
}

#' Build `.value ~ .t + <adjust> + (1 | .id)` or `(1 + .t | .id)`
#' @noRd
lm_lme_formula <- function(adjust = NULL, random_slope = FALSE) {
  fixed <- quote(.t)
  if (!is.null(adjust)) {
    fixed <- call("+", fixed, rlang::f_rhs(adjust))
  }
  re <- if (random_slope) quote((1 + .t | .id)) else quote((1 | .id))
  # The adjust formula's environment lets terms like splines::ns(age, 3) or a
  # user-defined function resolve as they would in a normal lmer() call.
  env <- if (is.null(adjust)) baseenv() else rlang::f_env(adjust)
  rlang::new_formula(quote(.value), call("+", fixed, re), env = env)
}

#' Fit lmer, collecting warnings and messages instead of printing them
#'
#' Returns list(fit, warnings, error); `fit` is NULL on error.
#' @noRd
lm_quiet_lmer <- function(formula, data, REML, control) {
  caught <- character(0)
  fit <- tryCatch(
    withCallingHandlers(
      lme4::lmer(formula, data = data, REML = REML, control = control),
      warning = function(w) {
        caught <<- c(caught, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        # lme4 reports singular fits as a message; isSingular() records it.
        invokeRestart("muffleMessage")
      }
    ),
    error = function(e) {
      structure(conditionMessage(e), class = "lm_lmer_error")
    }
  )
  if (inherits(fit, "lm_lmer_error")) {
    return(list(fit = NULL, warnings = caught, error = unclass(fit)))
  }
  list(fit = fit, warnings = trimws(caught), error = NULL)
}

#' One warning listing every LME fit that was not "ok"
#' @noRd
lm_warn_lme_issues <- function(lme_out, vars) {
  diags <- purrr::map2(lme_out, vars, \(o, v) dplyr::mutate(o$details, .var = v)) |>
    purrr::list_rbind()
  bad <- diags[diags$status != "ok", , drop = FALSE]
  if (nrow(bad) == 0L) {
    return(invisible())
  }
  lines <- paste0(
    bad$.var, " at landmark ", bad$.landmark, " (", bad$status, "): ", bad$message
  )
  # Messages from lme4 are inserted verbatim, so escape cli's glue braces.
  lines <- gsub("}", "}}", gsub("{", "{{", lines, fixed = TRUE), fixed = TRUE)
  n_failed <- sum(bad$status == "failed")
  cli::cli_warn(
    c(
      "{nrow(bad)} mixed model{?s} did not fit cleanly.",
      rlang::set_names(lines, rep("*", length(lines))),
      "i" = if (n_failed > 0L) "Failed fits give {.code NA} for that variable at that landmark, so complete-case models will use fewer subjects.",
      "i" = "Inspect the fits with {.fn landmark_lme_fits}."
    )
  )
}

#' Text form of the adjust formula for source ("~age + sex")
#' @noRd
lm_deparse_adjust <- function(adjust) {
  if (is.null(adjust)) NA_character_ else paste(deparse(adjust), collapse = " ")
}

#' Source of derived covariates
#'
#' Returns the tibble of derived-covariate records (one row per derived
#' column), or NULL if none have been added.
#' @noRd
lm_covariate_info <- function(x) {
  landmark_config(x)$covariates
}
