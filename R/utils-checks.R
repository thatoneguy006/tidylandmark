# Internal input-validation helpers ------------------------------------------
#
# None of these are exported. Each one takes:
#   * `arg`  - the user-facing argument name, used in error messages, and
#   * `call` - the environment of the user-facing function (e.g.
#              `landmark_data()`), so errors read "Error in `landmark_data()`"
#              instead of naming the internal helper.
# They either return a cleaned value or throw a classed rlang error via cli.

# Column names that tidylandmark creates and needs to function.
lm_reserved_cols <- c(".id", ".landmark", ".time", ".status", ".horizon")

#' Resolve a tidy-select argument to exactly one column name
#'
#' Uses tidyselect so `id`, `time` and `status` accept the same syntax as
#' `vars` in `add_landmark_covariates()`: bare names, strings, or
#' `all_of(x)`. A non-existent column triggers tidyselect's own error.
#' @noRd
lm_select_one <- function(data, quo, arg, call = rlang::caller_env()) {
  pos <- tidyselect::eval_select(
    quo,
    data = data,
    allow_rename = FALSE,
    error_call = call
  )
  if (length(pos) != 1L) {
    cli::cli_abort(
      "{.arg {arg}} must select exactly one column, not {length(pos)}.",
      call = call
    )
  }
  names(pos)
}

#' Refuse user input if columns `.id`, `.time`, ..., already exist.
#' Necessary so we don't overwrite user data. Will also prevent
#' running `landmark_data()` on an already created `landmark_tbl`.
#' @noRd
lm_check_reserved <- function(data, call = rlang::caller_env()) {
  clash <- intersect(names(data), lm_reserved_cols)
  if (length(clash) > 0L) {
    cli::cli_abort(
      c(
        "{.arg data} already contains reserved column name{?s} {.field {clash}}.",
        "i" = "Rename the conflicting column(s) before calling {.fn landmark_data}."
      ),
      call = call
    )
  }
  invisible(data)
}

#' Check for missing data in user-supplied data.
#' @noRd
lm_check_no_na <- function(x, arg, call = rlang::caller_env()) {
  n_na <- sum(is.na(x))
  if (n_na > 0L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must not contain missing values.",
        "x" = "Found {n_na} missing value{?s}."
      ),
      call = call
    )
  }
  invisible(x)
}

#' Check that only one row per subject exists. Longitudinal data will be added
#' via `add_landmark_covariates()`.
#' @noRd
lm_check_unique <- function(x, arg, call = rlang::caller_env()) {
  dup <- unique(x[duplicated(x)])
  if (length(dup) > 0L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must uniquely identify the rows of {.arg data}.",
        "x" = "Duplicated value{?s}: {.val {dup}}."
      ),
      call = call
    )
  }
  invisible(x)
}

#' Ensure time variable is valid.
#' @noRd
lm_check_time <- function(x, arg, call = rlang::caller_env()) {
  if (!is.numeric(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be a numeric column, not {.cls {class(x)}}.",
      call = call
    )
  }
  n_bad <- sum(!is.finite(x) | x < 0)
  if (n_bad > 0L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must contain non-negative, finite values.",
        "x" = "Found {n_bad} negative or infinite value{?s}."
      ),
      call = call
    )
  }
  invisible(x)
}

#' Validate and recode the status column to integer codes
#'
#' Returns a list with `codes` (integer: 0 = censored, 1, 2, ... = event
#' types), `labels` (NULL for numeric input; a named character vector mapping
#' code -> factor level for factor input) and `event_codes` (the event codes
#' that exist in the source coding).
#' @noRd
lm_code_status <- function(x, arg, call = rlang::caller_env()) {
  if (is.factor(x)) {
    lev <- levels(x)
    if (length(lev) < 2L) {
      cli::cli_abort(
        c(
          "A factor {.arg {arg}} must have at least two levels.",
          "i" = "The first level is censoring; later levels are event types."
        ),
        call = call
      )
    }
    labels <- lev
    names(labels) <- seq_along(lev) - 1L
    return(list(
      codes = as.integer(x) - 1L,
      labels = labels,
      event_codes = seq_len(length(lev) - 1L)
    ))
  }

  if (!is.numeric(x)) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be a numeric or factor column, not {.cls {class(x)}}.",
        "i" = "Use 0 for censoring and 1, 2, ... for event types."
      ),
      call = call
    )
  }

  is_bad <- !is.finite(x) | x < 0 | x != trunc(x)
  if (any(is_bad)) {
    bad <- unique(x[is_bad])
    cli::cli_abort(
      c(
        "{.arg {arg}} must be coded as 0 (censored) or a positive integer (event type).",
        "x" = "Invalid value{?s}: {.val {bad}}."
      ),
      call = call
    )
  }

  codes <- as.integer(x)
  list(
    codes = codes,
    labels = NULL,
    event_codes = sort(unique(codes[codes > 0L]))
  )
}

#' Validate and normalize the landmark times.
#' @noRd
lm_check_landmarks <- function(x, arg = "landmarks", call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) == 0L) {
    cli::cli_abort(
      "{.arg {arg}} must be a non-empty numeric vector.",
      call = call
    )
  }
  if (!all(is.finite(x))) {
    cli::cli_abort("{.arg {arg}} must be finite and non-missing.", call = call)
  }
  if (any(x < 0)) {
    cli::cli_abort("{.arg {arg}} must be non-negative.", call = call)
  }
  dup <- unique(x[duplicated(x)])
  if (length(dup) > 0L) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must not contain duplicate values.",
        "x" = "Duplicated: {.val {dup}}."
      ),
      call = call
    )
  }
  sort(as.double(x))
}

#' Validate user-supplied horizon.
#'
#' Will change to allow multiple horizons in a future update.
#' @noRd
lm_check_horizon <- function(x, arg = "horizon", call = rlang::caller_env()) {
  if (is.null(x)) {
    return(NULL)
  }
  # `||` short-circuits left to right, so `is.finite(x)` and `x <= 0` are only
  # evaluated once we know `x` is a single number; an NA therefore fails the
  # `is.finite()` test instead of reaching `if (NA)`.
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    cli::cli_abort(
      c(
        "{.arg {arg}} must be {.code NULL} or a single positive, finite number.",
        "i" = "Use {.code horizon = NULL} for unbounded follow-up."
      ),
      call = call
    )
  }
  as.double(x)
}

#' Validate a single TRUE/FALSE argument (e.g. `progress`).
#'
#' `isTRUE()`/`isFALSE()` accept only a length-1, non-missing logical, so
#' `NA`, `c(TRUE, FALSE)`, `1` and `"yes"` all fail here instead of
#' reaching an `if ()` later.
#' @noRd
lm_check_flag <- function(x, arg, call = rlang::caller_env()) {
  if (!isTRUE(x) && !isFALSE(x)) {
    cli::cli_abort(
      "{.arg {arg}} must be {.code TRUE} or {.code FALSE}.",
      call = call
    )
  }
  invisible(x)
}

#' Turn a `progress` flag into purrr's `.progress` argument.
#' @noRd
lm_progress_arg <- function(progress, label) {
  if (progress) label else FALSE
}
