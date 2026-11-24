# Checks used by fit_landmark() ----------------------------------------------
#
# Same conventions as utils-checks.R: internal, `arg` names the user-facing
# argument and `call` is the user-facing function's environment.

# Column names that survival::finegray() adds, plus the two helper columns
# fit_landmark() adds before calling it. A formula variable with one of these
# names would be overwritten in the expanded data.
lm_finegray_cols <- c("fgstart", "fgstop", "fgstatus", "fgwt", ".fg_status", ".fg_id")

#' Validate the covariate formula of fit_landmark()
#'
#' Returns a list with
#'   * `data_vars`: formula variables that are columns of `data` (these define
#'     the complete cases), and
#'   * `check_vars`: the subset used outside `strata()` terms, which must vary
#'     within each landmark for the model to be estimable.
#' Variables that are not columns of `data` must exist in the formula's
#' environment (e.g. a spline `df` set by the user); anything else is an error
#' here rather than one failed fit per landmark.
#' @noRd
lm_check_fit_formula <- function(formula, data, model, call = rlang::caller_env()) {
  if (!rlang::is_formula(formula)) {
    cli::cli_abort(
      "{.arg formula} must be a one-sided formula, such as {.code ~ age + sex}, not {.obj_type_friendly {formula}}.",
      call = call
    )
  }
  if (rlang::is_formula(formula, lhs = TRUE)) {
    cli::cli_abort(
      c(
        "{.arg formula} must be one-sided, such as {.code ~ age + sex}.",
        "i" = "Leave out {.code Surv()}: the outcome is built from {.field .time} and {.field .status}."
      ),
      call = call
    )
  }
  rhs  <- rlang::f_rhs(formula)
  vars <- all.vars(rhs)

  if ("." %in% vars) {
    cli::cli_abort(
      "{.arg formula} can't use {.code .}; list the covariates explicitly.",
      call = call
    )
  }

  # Function names are found in the deparsed text: `tt(` or `cluster(`,
  # optionally namespaced (`survival::cluster(`), but not e.g. `my_tt(`.
  rhs_text <- paste(deparse(rhs), collapse = " ")
  banned <- c("tt", "cluster")
  bad_calls <- banned[vapply(
    banned,
    \(fn) grepl(paste0("(^|[^.[:alnum:]_])", fn, "\\s*\\("), rhs_text),
    logical(1)
  )]
  if (length(bad_calls) > 0L) {
    cli::cli_abort(
      c(
        "{.arg formula} can't contain {.fn {bad_calls}} terms.",
        "i" = "Each landmark model has one row per subject, so there is nothing to cluster or vary over time."
      ),
      call = call
    )
  }

  reserved <- intersect(vars, lm_reserved_cols)
  if (length(reserved) > 0L) {
    cli::cli_abort(
      c(
        "{.arg formula} can't use reserved column{?s} {.field {reserved}}.",
        "i" = "The outcome is built from {.field .time} and {.field .status} internally."
      ),
      call = call
    )
  }

  env <- rlang::f_env(formula)
  not_cols <- setdiff(vars, names(data))
  unknown  <- not_cols[!vapply(not_cols, exists, logical(1), envir = env)]
  if (length(unknown) > 0L) {
    cli::cli_abort(
      c(
        "{.arg formula} uses variable{?s} not found in {.arg data}: {.field {unknown}}.",
        "i" = "Covariates from longitudinal data must be added first with {.fn add_landmark_covariates}."
      ),
      call = call
    )
  }

  data_vars <- intersect(vars, names(data))
  if (length(data_vars) == 0L) {
    cli::cli_abort(
      "{.arg formula} must use at least one column of {.arg data}.",
      call = call
    )
  }

  if (model == "finegray") {
    clash <- intersect(data_vars, lm_finegray_cols)
    if (length(clash) > 0L) {
      cli::cli_abort(
        c(
          "With {.code model = \"finegray\"}, {.arg formula} can't use column{?s} named {.field {clash}}.",
          "i" = "{.fn survival::finegray} creates columns with these names; rename yours first."
        ),
        call = call
      )
    }
  }

  # Variables that appear only inside strata() don't need to vary: a single
  # stratum is a valid model. attr(, "variables") is list(age, strata(sex), ...)
  # and specials$strata indexes into it (there is no response to skip).
  tt <- stats::terms(formula, specials = "strata")
  var_exprs  <- as.list(attr(tt, "variables"))[-1]
  strata_idx <- attr(tt, "specials")$strata
  non_strata <- if (is.null(strata_idx)) var_exprs else var_exprs[-strata_idx]
  check_vars <- intersect(unique(unlist(lapply(non_strata, all.vars))), data_vars)

  list(data_vars = data_vars, check_vars = check_vars)
}

#' Validate `cause` and return it as an integer event code
#'
#' `cause` may be an event code, or a level name when `status` was a factor.
#' @noRd
lm_check_cause <- function(cause, cfg, arg = "cause", call = rlang::caller_env()) {
  codes  <- cfg$event_codes
  labels <- cfg$status_labels
  # status_labels is named by code ("0", "1", ...); keep the event levels only.
  event_labels <- if (is.null(labels)) NULL else unname(labels[as.character(codes)])

  if (is.character(cause) && length(cause) == 1L && !is.na(cause)) {
    if (is.null(labels)) {
      cli::cli_abort(
        c(
          "{.arg {arg}} can only be a level name when {.arg status} was a factor.",
          "i" = "Use one of the event codes: {.val {codes}}."
        ),
        call = call
      )
    }
    if (!cause %in% event_labels) {
      cli::cli_abort(
        c(
          "{.arg {arg}} must be one of the event levels {.val {event_labels}}, not {.val {cause}}.",
          "i" = "The first level ({.val {labels[[1]]}}) is censoring."
        ),
        call = call
      )
    }
    return(codes[event_labels == cause])
  }

  if (!is.numeric(cause) || length(cause) != 1L || is.na(cause) || !cause %in% codes) {
    hint <- if (!is.null(labels)) c("i" = "Or one of the event levels: {.val {event_labels}}.")
    cli::cli_abort(
      c("{.arg {arg}} must be one of the event codes in the data: {.val {codes}}.", hint),
      call = call
    )
  }
  as.integer(cause)
}
