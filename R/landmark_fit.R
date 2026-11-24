# Methods for landmark_fit objects ----------------------------------------------

#' Tidy a landmark model fit
#'
#' Returns one row per landmark and model term, combining the coefficients of
#' every landmark model with that landmark's counts and fit status. Landmarks
#' that were not fitted (`"unestimable"` or `"failed"`) keep one row, with
#' `term` and the estimates `NA`, so they don't silently disappear from a
#' results table.
#'
#' Inference is Wald-based, as in `summary()` of a [survival::coxph()] fit:
#' `statistic` is `estimate / std.error`, `p.value` is two-sided and
#' confidence limits are `estimate` \eqn{\pm z_{1-\alpha/2}} `std.error` on the
#' log scale. Fine-Gray models use the robust standard error. Standard errors
#' treat covariates derived by [add_landmark_covariates()] as known values, so
#' they are too small for `lme` covariates.
#'
#' @param x A `landmark_fit`, as returned by [fit_landmark()].
#' @param exponentiate If `TRUE`, report hazard ratios (subdistribution hazard
#'   ratios for Fine-Gray models) instead of log hazard ratios. `std.error`
#'   stays on the log scale. Default `FALSE`.
#' @param conf.int If `TRUE`, add `conf.low` and `conf.high`. Default `FALSE`.
#' @param conf.level Confidence level for the interval. Default `0.95`.
#' @param ... Must be empty. Present for compatibility with the generic.
#'
#' @returns A tibble with columns
#'   * `.landmark`, `term`
#'   * `estimate`, `std.error`, `statistic`, `p.value`, and `conf.low`,
#'     `conf.high` when `conf.int = TRUE`
#'   * `origin`: how the term's variables were obtained: `"baseline"`,
#'     `"locf"`, `"lme"`, or a combination such as `"baseline:lme"` for an
#'     interaction
#'   * `n_eligible`, `n_model`, `n_event`, `n_competing`
#'   * `.model`, `.cause`, `.fit_status`, `.message`
#'
#' @seealso [fit_landmark()]
#' @examples
#' pbc_lm <- pbcseq_baseline |>
#'   landmark_data(id, futime, status, landmarks = c(365, 730), horizon = 1825) |>
#'   add_landmark_covariates(pbcseq_visits, id, day, locf = bili)
#' fit <- fit_landmark(pbc_lm, ~ age + log(bili), cause = 2)
#'
#' tidy(fit)
#' tidy(fit, exponentiate = TRUE, conf.int = TRUE)
#' @exportS3Method generics::tidy
tidy.landmark_fit <- function(x,
                              exponentiate = FALSE,
                              conf.int = FALSE,
                              conf.level = 0.95,
                              ...) {
  call <- rlang::current_env()
  rlang::check_dots_empty(call = call)
  lm_check_flag(exponentiate, "exponentiate", call)
  lm_check_flag(conf.int, "conf.int", call)
  lm_check_conf_level(conf.level, call = call)

  co <- x$coefficients
  co$statistic <- co$estimate / co$std.error
  co$p.value   <- 2 * stats::pnorm(-abs(co$statistic))
  if (conf.int) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    co$conf.low  <- co$estimate - z * co$std.error
    co$conf.high <- co$estimate + z * co$std.error
  }
  if (exponentiate) {
    co$estimate <- exp(co$estimate)
    if (conf.int) {
      co$conf.low  <- exp(co$conf.low)
      co$conf.high <- exp(co$conf.high)
    }
  }

  # A placeholder row for each landmark without coefficients; bind_rows()
  # fills the other columns with NA.
  info <- x$model_info
  placeholder <- tibble::tibble(.landmark = setdiff(info$.landmark, co$.landmark))
  out <- dplyr::bind_rows(co, placeholder) |>
    dplyr::left_join(info, by = ".landmark", relationship = "many-to-one") |>
    dplyr::arrange(.data$.landmark)

  cols <- c(
    ".landmark", "term", "estimate", "std.error", "statistic", "p.value",
    if (conf.int) c("conf.low", "conf.high"),
    "origin", "n_eligible", "n_model", "n_event", "n_competing",
    ".model", ".cause", ".fit_status", ".message"
  )
  out <- out[cols]

  if (any(grepl("lme", out$origin, fixed = TRUE))) {
    cli::cli_inform(
      c(
        "i" = "Standard errors treat LME-derived covariates as known, so they understate uncertainty.",
        " " = "Use {.fn bootstrap_landmark} for intervals that account for this."
      ),
      .frequency = "once",
      .frequency_id = "tidylandmark_lme_se"
    )
  }
  out
}

#' Print a landmark model fit
#'
#' Shows the model, formula and horizon, then one line per landmark with its
#' counts and fit status.
#'
#' @param x A `landmark_fit`, as returned by [fit_landmark()].
#' @param ... Ignored.
#' @returns `x`, invisibly.
#' @export
print.landmark_fit <- function(x, ...) {
  cfg <- x$config
  model_txt <- switch(
    cfg$model_label,
    "cox"                = "Cox",
    "cause-specific cox" = "Cause-specific Cox",
    "finegray"           = "Fine-Gray"
  )
  cause_txt <- if (is.null(cfg$cause_label)) {
    cfg$cause
  } else {
    paste0(cfg$cause, " (", cfg$cause_label, ")")
  }
  horizon <- landmark_config(x$data)$horizon
  hz_txt <- if (is.null(horizon)) {
    "none (follow-up to observed time)"
  } else {
    paste0(horizon, " (relative to each landmark)")
  }

  info <- x$model_info
  writeLines(c(
    paste0(
      "<landmark_fit> ", model_txt, " model for cause ", cause_txt,
      " at ", nrow(info), " landmark(s)"
    ),
    paste0("Formula: ", paste(deparse(cfg$formula), collapse = " ")),
    paste0("Horizon: ", hz_txt)
  ))

  shown <- c(".landmark", "n_eligible", "n_model", "n_event",
             if (cfg$model_label != "cox") "n_competing", ".fit_status")
  table <- as.data.frame(info[shown])
  names(table) <- sub("^\\.", "", names(table))
  print(table, row.names = FALSE)

  notes <- character(0)
  if (any(info$.fit_status != "ok")) {
    notes <- c(notes, "Some landmarks did not fit cleanly; see `$model_info$.message`.")
  }
  if (length(cfg$lme_columns) > 0L) {
    notes <- c(notes, paste0(
      "Standard errors treat LME-derived covariates (",
      paste(cfg$lme_columns, collapse = ", "), ") as known."
    ))
  }
  notes <- c(notes, "Use tidy() for coefficients and summary() for hazard ratios.")
  cat("\n")
  writeLines(notes)
  invisible(x)
}

#' Summarise a landmark model fit
#'
#' A shortcut for `tidy(object, exponentiate = TRUE, conf.int = TRUE)`: hazard
#' ratios (or subdistribution hazard ratios) with confidence intervals for
#' every landmark and term.
#'
#' @param object A `landmark_fit`, as returned by [fit_landmark()].
#' @param conf.level Confidence level. Default `0.95`.
#' @param ... Must be empty.
#' @returns A tibble; see [tidy.landmark_fit()].
#' @export
summary.landmark_fit <- function(object, conf.level = 0.95, ...) {
  rlang::check_dots_empty(call = rlang::current_env())
  tidy(object, exponentiate = TRUE, conf.int = TRUE, conf.level = conf.level)
}

#' Validate a confidence level: one number strictly between 0 and 1
#' @noRd
lm_check_conf_level <- function(x, arg = "conf.level", call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    cli::cli_abort(
      "{.arg {arg}} must be a single number between 0 and 1, such as {.code 0.95}.",
      call = call
    )
  }
  invisible(x)
}
