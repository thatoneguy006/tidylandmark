#' Fit a survival model at every landmark
#'
#' @description
#' `fit_landmark()` fits one independent regression model per landmark of a
#' `landmark_tbl`: a Cox proportional hazards model ([survival::coxph()]) or a
#' Fine-Gray subdistribution hazards model ([survival::finegray()] followed by
#' a weighted [survival::coxph()]). The outcome is built internally from
#' `.time` and `.status`, so `formula` only lists the covariates.
#'
#' @param data A `landmark_tbl`, as returned by [landmark_data()] or
#'   [add_landmark_covariates()].
#' @param formula A one-sided formula of covariates, e.g. `~ age + sex + bili`.
#'   Any right-hand side that [survival::coxph()] accepts can be used,
#'   including factors, interactions, transformations (`log(bili)`), splines
#'   (`splines::ns(age, 3)`) and `strata()`. Leave out `Surv()`. The formula
#'   can't use `.`, reserved columns (`.time`, `.status`, ...), `tt()` or
#'   `cluster()`.
#' @param model `"cox"` (default) or `"finegray"`.
#' @param cause The event of interest: one of the event codes of the `status`
#'   given to [landmark_data()] (e.g. `1`), or, if `status` was a factor, one
#'   of its event levels (e.g. `"relapse"`). Default `1`.
#' @param progress Show a progress bar while fitting? Defaults to the
#'   `tidylandmark.progress` option, or [interactive()] if it is unset.
#'
#' @details
#' ## Models
#' With `model = "cox"`, each landmark gets
#' `coxph(Surv(.time, .status == cause) ~ <formula>)` with Efron's method for
#' ties. If the status coding has a single event type, this is an ordinary
#' Cox model. If it has other event types, they are treated as censored at
#' their event time, so the model estimates **cause-specific** hazard ratios;
#' the `.model` column then reads `"cause-specific cox"`.
#'
#' With `model = "finegray"`, [survival::finegray()] expands each landmark's
#' data with time-dependent censoring weights for `cause`, and a weighted
#' [survival::coxph()] is fitted to the expanded data with a robust variance
#' clustered on subject, as the [survival::finegray()] documentation
#' recommends. Coefficients are log subdistribution hazard ratios. The
#' censoring distribution is estimated on each landmark's complete-case model
#' sample and is not stratified; `strata()` terms in `formula` stratify the
#' weighted Cox model only.
#'
#' ## Missing covariates and failed fits
#' Each landmark's model uses the rows with no missing values in the variables
#' of `formula`; `n_eligible` and `n_model` in the results show how many rows
#' that leaves. A landmark is not fitted (status `"unestimable"`) when it has
#' no subjects at risk, no complete cases, no target events in the model
#' sample, or a covariate with a single observed value. A fit that warns (for
#' example about convergence, or a coefficient that may be infinite) gets
#' status `"warning"`, and a fit that errors gets `"failed"`. Neither affects
#' the other landmarks, and one warning at the end lists every landmark that
#' did not fit cleanly. A Fine-Gray model at a landmark with no competing
#' events is still fitted (it then equals the Cox model) with status
#' `"warning"`. If a covariate made by an `lme` model in
#' [add_landmark_covariates()] did not fit cleanly at a landmark, that
#' landmark's status is at least `"warning"` and its message names the
#' covariate.
#'
#' ## Interpretation
#' Each model is conditional on being event-free and under follow-up at its
#' landmark, so different landmarks describe different populations. Fits at
#' different landmarks share subjects and are not independent, so differences
#' between landmarks are not tested. Standard errors treat covariates derived
#' by [add_landmark_covariates()] as known values; for `lme` covariates they
#' are therefore too small, because the mixed models were estimated from the
#' same subjects.
#'
#' @returns A `landmark_fit` object: a list with
#'   * `data`: the `landmark_tbl` that was fitted.
#'   * `models`: one fitted [survival::coxph()] object per landmark, named by
#'     landmark (`NULL` where no model was fitted). For Fine-Gray models this
#'     is the weighted fit to the expanded data.
#'   * `model_info`: a tibble with one row per landmark: `.landmark`,
#'     `.model`, `.cause`, `n_eligible` (rows at risk), `n_model` (complete
#'     cases used), `n_event` and `n_competing` (target and competing events
#'     in the model sample), `.fit_status` (`"ok"`, `"warning"`,
#'     `"unestimable"` or `"failed"`) and `.message` (`NA` when `"ok"`).
#'   * `coefficients`: a tibble with one row per landmark and model term:
#'     `.landmark`, `term`, `estimate` and `std.error` (log scale) and
#'     `origin` (`"baseline"`, `"locf"`, `"lme"`, or e.g. `"baseline:lme"`
#'     for an interaction).
#'   * `warnings`: a tibble of the warnings captured from each fit.
#'   * `config`: the formula, model, cause and related settings.
#'
#'   Use [tidy()][tidy.landmark_fit] for a results table and `print()` or
#'   `summary()` for an overview.
#'
#' @seealso [tidy.landmark_fit()], [landmark_data()],
#'   [add_landmark_covariates()]
#' @export
#'
#' @examples
#' # PBC: status 0 = censored, 1 = transplant, 2 = death
#' pbc_lm <- pbcseq_baseline |>
#'   landmark_data(id, futime, status, landmarks = c(365, 730), horizon = 1825) |>
#'   add_landmark_covariates(pbcseq_visits, id, day, locf = c(bili, albumin))
#'
#' # Cause-specific Cox model for death (transplant is a competing event)
#' fit <- fit_landmark(pbc_lm, ~ age + log(bili) + albumin, cause = 2)
#' fit
#' tidy(fit, exponentiate = TRUE, conf.int = TRUE)
#'
#' # Fine-Gray model for death
#' fit_fg <- fit_landmark(
#'   pbc_lm, ~ age + log(bili) + albumin, model = "finegray", cause = 2
#' )
#' summary(fit_fg)
#'
#' # The native coxph() fit for one landmark
#' fit$models[["365"]]
fit_landmark <- function(data,
                         formula,
                         model = c("cox", "finegray"),
                         cause = 1,
                         progress = getOption("tidylandmark.progress", interactive())) {
  call <- rlang::current_env()

  # ---- 1. Validate ----------------------------------------------------------
  # Re-check the reserved columns rather than trusting the stored config:
  # dplyr verbs between pipeline steps keep the attribute but can break the
  # data.
  lm_check_intact(data, "data", call)
  rlang::check_required(formula)
  model <- rlang::arg_match(model, error_call = call)
  lm_check_flag(progress, "progress", call)

  cfg   <- landmark_config(data)
  vars  <- lm_check_fit_formula(formula, data, model, call)
  cause <- lm_check_cause(cause, cfg, call = call)

  model_label <- if (model == "finegray") {
    "finegray"
  } else if (any(cfg$event_codes != cause)) {
    "cause-specific cox"
  } else {
    "cox"
  }

  # ---- 2. Covariate provenance ----------------------------------------------
  cov_info   <- lm_present_covariates(data)
  origins    <- lm_term_origins(formula, vars$data_vars, cov_info)
  lme_issues <- lm_lme_issues(cov_info, vars$data_vars)

  # ---- 3. One model per landmark --------------------------------------------
  # Same landmark list as print(): after filter(.landmark == 30) only 30 is
  # fitted; landmarks built with no subjects at risk are kept so they are
  # reported as unestimable.
  landmarks <- lm_landmarks_shown(data$.landmark, cfg)
  rows <- tibble::as_tibble(data)

  fits <- purrr::map(
    landmarks,
    \(s) {
      lm_fit_one(
        rows[rows$.landmark == s, , drop = FALSE],
        s           = s,
        formula     = formula,
        vars        = vars,
        model       = model,
        cause       = cause,
        event_codes = cfg$event_codes,
        lme_issues  = lme_issues[lme_issues$.landmark == s, , drop = FALSE],
        keep_model  = TRUE
      )
    },
    .progress = lm_progress_arg(progress, "Fitting survival models")
  )

  # ---- 4. Assemble ----------------------------------------------------------
  model_info <- purrr::list_rbind(purrr::map(fits, "info")) |>
    tibble::add_column(.model = model_label, .cause = cause, .after = ".landmark")

  coefficients <- purrr::list_rbind(purrr::map(fits, "coefs"))
  coefficients$origin <- unname(origins[coefficients$.term_label])
  coefficients$.term_label <- NULL

  warn_list <- purrr::map(fits, "warnings")
  warnings <- tibble::tibble(
    .landmark = rep(landmarks, lengths(warn_list)),
    message   = as.character(unlist(warn_list))
  )

  lm_warn_fit_issues(model_info)

  cause_label <- if (is.null(cfg$status_labels)) NULL else cfg$status_labels[[as.character(cause)]]
  new_landmark_fit(
    data         = data,
    models       = rlang::set_names(purrr::map(fits, "model"), as.character(landmarks)),
    model_info   = model_info,
    coefficients = coefficients,
    warnings     = warnings,
    config       = list(
      formula     = formula,
      model       = model,
      model_label = model_label,
      cause       = cause,
      cause_label = cause_label,
      data_vars   = vars$data_vars,
      lme_columns = intersect(vars$data_vars, cov_info$column[cov_info$method == "lme"]),
      origins     = origins
    )
  )
}

#' Low-level constructor (no validation)
#' @noRd
new_landmark_fit <- function(data, models, model_info, coefficients, warnings, config) {
  structure(
    list(
      data         = data,
      models       = models,
      model_info   = model_info,
      coefficients = coefficients,
      warnings     = warnings,
      config       = config
    ),
    class = "landmark_fit"
  )
}

# One landmark -------------------------------------------------------------------

#' Fit the model for one landmark
#'
#' `d` holds that landmark's rows (a plain tibble). Never errors or prints:
#' problems become a status and message, so one landmark can't stop the
#' others, and the future bootstrap can call this once per replicate with
#' `keep_model = FALSE`.
#'
#' Returns list(info, coefs, model, warnings): a one-row info tibble, a
#' coefficient tibble (.landmark, term, .term_label, estimate, std.error), the
#' coxph fit (or NULL), and the captured warning messages.
#' @noRd
lm_fit_one <- function(d, s, formula, vars, model, cause, event_codes,
                       lme_issues = NULL, keep_model = TRUE) {
  data_vars <- vars$data_vars
  cc   <- if (nrow(d) > 0L) stats::complete.cases(d[data_vars]) else logical(0)
  d_cc <- d[cc, , drop = FALSE]

  info <- tibble::tibble(
    .landmark   = s,
    n_eligible  = nrow(d),
    n_model     = nrow(d_cc),
    n_event     = sum(d_cc$.status == cause),
    n_competing = sum(d_cc$.status > 0L & d_cc$.status != cause),
    .fit_status = "ok",
    .message    = NA_character_
  )

  # Problems with lme covariates at this landmark are reported whatever
  # happens to the survival model.
  upstream <- character(0)
  if (!is.null(lme_issues) && nrow(lme_issues) > 0L) {
    upstream <- paste0(
      "lme covariate ", lme_issues$column, " (", lme_issues$status, "): ",
      lme_issues$message
    )
  }

  finish <- function(status, messages = character(0), fit = NULL,
                     warnings = character(0)) {
    messages <- unique(c(messages, upstream))
    if (status == "ok" && length(messages) > 0L) {
      status <- "warning"
    }
    info$.fit_status <- status
    info$.message    <- if (length(messages) > 0L) paste(messages, collapse = "; ") else NA_character_
    coefs <- if (is.null(fit)) lm_empty_coefs() else lm_coef_table(fit)
    coefs <- tibble::add_column(coefs, .landmark = rep(s, nrow(coefs)), .before = 1)
    list(
      info     = info,
      coefs    = coefs,
      model    = if (keep_model) fit else NULL,
      warnings = warnings
    )
  }

  # ---- Checks that make a model meaningless --------------------------------
  if (info$n_eligible == 0L) {
    return(finish("unestimable", "no subjects at risk"))
  }
  if (info$n_model == 0L) {
    return(finish("unestimable", "no complete cases for the formula variables"))
  }
  if (info$n_event == 0L) {
    return(finish(
      "unestimable",
      sprintf("no target events (cause %d) in the model sample", cause)
    ))
  }
  constant <- purrr::keep(
    vars$check_vars,
    \(v) dplyr::n_distinct(d_cc[[v]], na.rm = TRUE) < 2L
  )
  if (length(constant) > 0L) {
    return(finish(
      "unestimable",
      paste0(constant, " has a single observed value", collapse = "; ")
    ))
  }

  # ---- Fit -------------------------------------------------------------------
  res <- lm_capture(
    lm_fit_backend(d_cc, formula, data_vars, model, cause, event_codes,
                   keep_frame = keep_model)
  )
  if (!is.null(res$error)) {
    return(finish("failed", res$error, warnings = res$warnings))
  }

  notes <- res$warnings
  if (model == "finegray" && info$n_competing == 0L) {
    notes <- c(notes, "no competing events: estimates equal the Cox model")
  }
  finish("ok", notes, fit = res$value, warnings = res$warnings)
}

#' Call the survival backend
#'
#' Calls are built with call2() and evaluated here, so the stored coxph call
#' shows the actual formula (not `formula = f`) and no bare column names
#' appear in package code.
#'
#' `keep_frame = TRUE` stores the model frame in the fit (`model = TRUE`).
#' The data only exist inside this function, so without it functions that
#' rebuild the model frame from the call, such as survival::cox.zph() and
#' some residuals(), would fail on the returned fit.
#' @noRd
lm_fit_backend <- function(d, formula, data_vars, model, cause, event_codes,
                           keep_frame = TRUE) {
  rhs <- rlang::f_rhs(formula)
  # A child of the user's formula environment that also holds Surv() and
  # strata(). Terms like log(bili), splines::ns(age, df) or the user's own
  # functions then resolve as they would in a direct coxph() call, and the
  # outcome can be written as plain `Surv(...)` (which survival's special-term
  # handling expects) without survival being attached.
  surv_fns <- list(Surv = survival::Surv, strata = survival::strata)
  env <- rlang::new_environment(surv_fns, parent = rlang::f_env(formula))

  if (model == "cox") {
    # Inline the cause value, e.g. Surv(.time, .status == 1L) ~ age.
    lhs <- rlang::call2("Surv", quote(.time), rlang::call2("==", quote(.status), cause))
    f  <- rlang::new_formula(lhs, rhs, env = env)
    cl <- rlang::call2(quote(survival::coxph), f, data = quote(d), model = keep_frame)
    return(rlang::eval_bare(cl, rlang::current_env()))
  }

  # Fine-Gray. finegray() needs a factor status whose first level is
  # censoring; `.fg_id` identifies subjects in the expanded data for the
  # robust variance. Only the outcome and formula columns are passed, so
  # `~ .` carries exactly the covariates.
  fg_in <- as.data.frame(d[c(".time", data_vars)])
  fg_in$.fg_status <- factor(d$.status, levels = sort(unique(c(0L, event_codes, d$.status))))
  fg_in$.fg_id     <- seq_len(nrow(fg_in))
  fg_formula <- rlang::new_formula(
    quote(Surv(.time, .fg_status)), quote(.),
    env = rlang::new_environment(surv_fns, parent = rlang::current_env())
  )
  fg_data <- survival::finegray(fg_formula, data = fg_in, etype = as.character(cause))

  f  <- rlang::new_formula(quote(Surv(fgstart, fgstop, fgstatus)), rhs, env = env)
  cl <- rlang::call2(
    quote(survival::coxph), f,
    data = quote(fg_data), weights = quote(fgwt), cluster = quote(.fg_id),
    model = keep_frame
  )
  rlang::eval_bare(cl, rlang::current_env())
}

#' Evaluate `expr`, collecting warnings instead of printing them
#'
#' Returns list(value, warnings, error); `value` is NULL and `error` the
#' message when `expr` errors. Same pattern as lm_quiet_lmer().
#' @noRd
lm_capture <- function(expr) {
  caught <- character(0)
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        caught <<- c(caught, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) structure(conditionMessage(e), class = "lm_capture_error")
  )
  if (inherits(value, "lm_capture_error")) {
    return(list(value = NULL, warnings = trimws(caught), error = unclass(value)))
  }
  list(value = value, warnings = trimws(caught), error = NULL)
}

# Coefficients -------------------------------------------------------------------

#' Log-scale estimates and standard errors from a coxph fit
#'
#' `fit$var` is the robust variance when the model was clustered (Fine-Gray)
#' and the model-based variance otherwise. Aliased coefficients (NA) get an NA
#' standard error. `.term_label` is the formula term each design column comes
#' from, used to look up its origin.
#' @noRd
lm_coef_table <- function(fit) {
  beta <- stats::coef(fit)
  p <- length(beta)
  if (p == 0L) {
    return(lm_empty_coefs())
  }
  ok <- !is.na(beta)
  se <- rep(NA_real_, p)
  v  <- fit$var
  if (is.matrix(v)) {
    # Depending on the survival version, the variance covers either all
    # columns or only the non-aliased ones.
    if (nrow(v) == p) {
      se <- sqrt(diag(v))
    } else if (nrow(v) == sum(ok)) {
      se[ok] <- sqrt(diag(v))
    }
  }
  se[!ok] <- NA_real_

  tibble::tibble(
    term        = names(beta),
    .term_label = lm_column_terms(fit, p),
    estimate    = unname(beta),
    std.error   = unname(se)
  )
}

#' Zero-row coefficient table
#' @noRd
lm_empty_coefs <- function() {
  tibble::tibble(
    term = character(), .term_label = character(),
    estimate = double(), std.error = double()
  )
}

#' Formula term label for each design column of a coxph fit
#'
#' coxph stores `assign` as a list named by term label, giving the design
#' columns of each term (e.g. list(age = 1, sex = 2, `ns(age, 2)` = 3:4)).
#' @noRd
lm_column_terms <- function(fit, p) {
  out <- rep(NA_character_, p)
  a <- fit$assign
  if (is.list(a) && length(a) > 0L) {
    idx  <- unlist(a, use.names = FALSE)
    labs <- rep(names(a), lengths(a))
    keep <- idx >= 1L & idx <= p
    out[idx[keep]] <- labs[keep]
  } else if (is.numeric(a) && length(a) == p) {
    # Fallback for model.matrix-style integer assign vectors.
    term_labels <- attr(stats::terms(fit), "term.labels")
    out <- term_labels[a]
  }
  out
}

#' Origin of each formula term: "baseline", "locf", "lme" or a combination
#'
#' Returns a character vector named by term label. A term's origin comes from
#' the data columns it uses: derived covariates have the method recorded by
#' add_landmark_covariates(), everything else is "baseline". Interactions of
#' different origins join them, e.g. "baseline:lme".
#' @noRd
lm_term_origins <- function(formula, data_vars, cov_info) {
  fac <- attr(stats::terms(formula), "factors")
  if (length(fac) == 0L) {
    return(character(0))
  }
  # Rows of `fac` are variable expressions ("ns(bili, 3)"), columns are terms.
  row_vars <- lapply(rownames(fac), \(e) intersect(all.vars(str2lang(e)), data_vars))
  method <- rlang::set_names(cov_info$method, cov_info$column)

  labs <- colnames(fac)
  out <- vapply(labs, \(lab) {
    used <- unique(unlist(row_vars[fac[, lab] > 0L]))
    if (length(used) == 0L) {
      return("baseline")
    }
    o <- ifelse(used %in% names(method), method[used], "baseline")
    paste(sort(unique(o)), collapse = ":")
  }, character(1))
  rlang::set_names(unname(out), labs)
}

#' LME fits that were not "ok", for lme covariates used in the formula
#'
#' One row per (column, landmark) problem: column, .landmark, status, message.
#' @noRd
lm_lme_issues <- function(cov_info, data_vars) {
  empty <- tibble::tibble(
    column = character(), .landmark = double(),
    status = character(), message = character()
  )
  lme <- cov_info[cov_info$method == "lme" & cov_info$column %in% data_vars, , drop = FALSE]
  if (nrow(lme) == 0L) {
    return(empty)
  }
  out <- purrr::map2(lme$column, lme$details, \(col, det) {
    tibble::tibble(
      column    = rep(col, nrow(det)),
      .landmark = det$.landmark,
      status    = det$status,
      message   = det$message
    )
  }) |>
    purrr::list_rbind()
  out[out$status != "ok", , drop = FALSE]
}

#' One warning listing every landmark that did not fit cleanly
#' @noRd
lm_warn_fit_issues <- function(model_info) {
  bad <- model_info[model_info$.fit_status != "ok", , drop = FALSE]
  if (nrow(bad) == 0L) {
    return(invisible())
  }
  lines <- paste0(
    "Landmark ", bad$.landmark, " (", bad$.fit_status, "): ", bad$.message
  )
  # Backend messages are inserted verbatim, so escape cli's glue braces.
  lines <- gsub("}", "}}", gsub("{", "{{", lines, fixed = TRUE), fixed = TRUE)
  n_bad <- nrow(bad)
  n_all <- nrow(model_info)
  cli::cli_warn(
    c(
      "{n_bad} of {n_all} landmark model{?s} did not fit cleanly.",
      rlang::set_names(lines, rep("*", length(lines))),
      "i" = "Unestimable and failed landmarks have no estimates. See {.code $model_info} for details."
    )
  )
}
