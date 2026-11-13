#' Construct landmark risk sets
#'
#' @description
#' `landmark_data()` builds the outcome side of a landmark analysis. For each
#' landmark time \eqn{s} it keeps every subject who is still under observation
#' and event-free at \eqn{s}. It moves the time origin to \eqn{s} and,
#' optionally, administratively censors follow-up at a fixed post-landmark
#' `horizon`. The result stacks all landmarks into one long table with one row
#' per eligible subject-landmark pair.
#'
#' Longitudinal covariates are not handled here. Attach them afterwards with
#' `add_landmark_covariates()`.
#'
#' @param data A data frame with **one row per subject**, holding follow-up
#'   time, event status and any baseline predictors.
#' @param id <[`tidy-select`][tidyselect::language]> The subject identifier
#'   column. It must be unique and non-missing. Its type is preserved in `.id`.
#' @param time <[`tidy-select`][tidyselect::language]> Numeric, non-negative,
#'   finite time from the common origin (time zero) to the first event or
#'   censoring.
#' @param status <[`tidy-select`][tidyselect::language]> Event indicator. It is
#'   either:
#'   * **numeric**: `0` for censoring and positive integers `1, 2, ...` for
#'     event types. Codes other than the target cause are treated as competing
#'     events by `fit_landmark()`.
#'   * **factor**: the first level is censoring and later levels are event
#'     types. They are recoded to `0, 1, 2, ...` in level order. The mapping is
#'     stored in the object's configuration.
#' @param landmarks A non-empty numeric vector of unique, non-negative landmark
#'   times on the original time scale. Landmarks are processed in ascending
#'   order regardless of input order.
#' @param horizon `NULL` (default) for follow-up to each subject's observed
#'   time, or a single positive number giving the follow-up window *relative to
#'   each landmark*. For example, `horizon = 180` at landmark `90` ends follow-up
#'   at original time `270`.
#'
#' @details
#' **Eligibility.** A subject is in the risk set at landmark \eqn{s} exactly when
#' \eqn{T_i > s}. A subject whose event or censoring is recorded *at* \eqn{s} is
#' not eligible.
#'
#' **Time and status.** For eligible subjects, `.time` is \eqn{T_i - s} when
#' `horizon` is `NULL`, and \eqn{\min(T_i - s, w)} otherwise. `.status` keeps
#' the original code when \eqn{T_i - s \le w}, so an event exactly at the
#' horizon is retained. Otherwise `.status` becomes `0`.
#'
#' **Columns.** The original `id`, `time` and `status` columns are replaced by
#' the reserved columns below. All other columns are carried over unchanged,
#' and subjects are never dropped for missing baseline predictors. If `data`
#' already contains a reserved column name, the function throws an error.
#'
#' A landmark with no eligible subjects triggers a warning. It contributes zero
#' rows but still appears in the configuration and in the printed counts.
#'
#' @returns A `landmark_tbl`, which is a tibble subclass. Its leading columns are:
#'   * `.id`: subject identifier, with the original type.
#'   * `.landmark`: landmark time on the original scale.
#'   * `.time`: follow-up since the landmark. It is always strictly positive.
#'   * `.status`: integer, `0` for censored or `1, 2, ...` for the event type.
#'   * `.horizon`: the configured horizon, or `NA` when follow-up is unbounded.
#'
#'   The remaining baseline columns follow. Rows are ordered by landmark, then
#'   by the original row order of `data`. Use [tibble::as_tibble()] to drop the
#'   class and metadata.
#'
#' @export
#'
#' @examples
#' subjects <- data.frame(
#'   patient_id    = 1:10,
#'   followup_days = c(20, 30, 48, 65, 85, 110, 145, 190, 250, 270),
#'   event         = c(1, 2, 1, 2, 0, 1, 2, 1, 0, 1),
#'   age           = c(12, 15, 11, 17, 9, 14, 13, 16, 10, 15)
#' )
#'
#' lm_dat <- landmark_data(
#'   subjects,
#'   id = patient_id,
#'   time = followup_days,
#'   status = event,
#'   landmarks = c(30, 60, 90),
#'   horizon = 180
#' )
#' lm_dat
#'
#' # Without a horizon, follow-up runs to each subject's observed time
#' landmark_data(subjects, patient_id, followup_days, event, landmarks = 60)
landmark_data <- function(data,
                          id,
                          time,
                          status,
                          landmarks,
                          horizon = NULL) {
  call <- rlang::current_env()

  # ---- 1. Validate the container and required arguments ---------------------
  if (!is.data.frame(data)) {
    cli::cli_abort(
      "{.arg data} must be a data frame, not {.cls {class(data)}}.",
      call = call
    )
  }
  rlang::check_required(id)
  rlang::check_required(time)
  rlang::check_required(status)
  rlang::check_required(landmarks)

  # as_tibble() drops dplyr grouping and row names but keeps column types.
  data <- tibble::as_tibble(data)
  lm_check_reserved(data, call = call)

  # ---- 2. Resolve tidy-select arguments to column names ---------------------
  # imap_chr() passes each quosure and its name ("id", "time", "status"), so
  # error messages name the right argument.
  role_cols <- purrr::imap_chr(
    list(
      id     = rlang::enquo(id),
      time   = rlang::enquo(time),
      status = rlang::enquo(status)
    ),
    \(quo, arg) lm_select_one(data, quo, arg, call)
  )
  if (anyDuplicated(role_cols) > 0L) {
    cli::cli_abort(
      "{.arg id}, {.arg time}, and {.arg status} must refer to three different columns.",
      call = call
    )
  }
  id_col     <- role_cols[["id"]]
  time_col   <- role_cols[["time"]]
  status_col <- role_cols[["status"]]

  # ---- 3. Validate column contents and scalar arguments ---------------------
  purrr::iwalk(role_cols, \(col, arg) lm_check_no_na(dplyr::pull(data, !!col), arg, call))
  lm_check_unique(dplyr::pull(data, !!id_col), "id", call)
  lm_check_time(dplyr::pull(data, !!time_col), "time", call)
  status_info <- lm_code_status(dplyr::pull(data, !!status_col), "status", call)

  landmarks <- lm_check_landmarks(landmarks, "landmarks", call)
  horizon   <- lm_check_horizon(horizon, "horizon", call)

  # ---- 4. Build one risk set per landmark -----------------------------------
  # Standardize the subject-level columns once: `.id` replaces the id column
  # and `.status` holds the validated integer codes.
  subject_data <- data |>
    dplyr::rename(.id = dplyr::all_of(id_col)) |>
    dplyr::mutate(.status = status_info$codes) |>
    dplyr::select(-dplyr::all_of(status_col))

  # `.env$s` forces the landmark value from the function environment. A bare
  # `s` would silently use a column called `s` if `data` happened to have one.
  risk_set_list <- purrr::map(landmarks, \(s) {
    subject_data |>
      dplyr::filter(.data[[time_col]] > .env$s) |>          # T > s: T == s is ineligible
      dplyr::mutate(
        .landmark = .env$s,
        .time     = as.double(.data[[time_col]] - .env$s)   # shift origin to the landmark
      )
  })

  n_eligible <- purrr::map_int(risk_set_list, nrow) |> rlang::set_names(landmarks)

  empty <- landmarks[n_eligible == 0L]
  if (length(empty) > 0L) {
    cli::cli_warn(
      c(
        "No subjects are at risk at landmark{?s} {.val {empty}}.",
        "i" = "Check that your supplied landmarks are correct."
      )
    )
  }

  risk_sets <- purrr::list_rbind(risk_set_list)

  # Administrative censoring. Both expressions read the *unclipped* `.time`
  # because mutate() evaluates `.status` before `.time` is overwritten.
  # Strict `>` keeps an event that falls exactly on the horizon.
  if (!is.null(horizon)) {
    risk_sets <- risk_sets |>
      dplyr::mutate(
        .status = dplyr::if_else(.data$.time > .env$horizon, 0L, .data$.status),
        .time   = pmin(.data$.time, .env$horizon)
      )
  }

  # ---- 5. Assemble output ---------------------------------------------------
  out <- risk_sets |>
    dplyr::mutate(.horizon = .env$horizon %||% NA_real_) |>
    dplyr::select(-dplyr::all_of(time_col)) |>
    dplyr::relocate(dplyr::all_of(lm_reserved_cols))

  config <- list(
    id            = id_col,
    time          = time_col,
    status        = status_col,
    landmarks     = landmarks,
    horizon       = horizon,
    status_labels = status_info$labels,
    event_codes   = status_info$event_codes,
    n_subjects    = nrow(data),
    n_eligible    = n_eligible
  )

  new_landmark_tbl(out, config = config)
}
