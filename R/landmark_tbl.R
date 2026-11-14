# The landmark_tbl class ------------------------------------------------------
#
# A landmark_tbl is a tibble with:
#   * class c("landmark_tbl", "tbl_df", "tbl", "data.frame")
#   * a "landmark_config" attribute (list) recording how it was built.
#
# The tibble `[` method and dplyr's default dplyr_reconstruct() both copy
# attributes from the input, so the class and config survive filter(),
# mutate(), etc. Those verbs can still invalidate the data, so downstream
# functions (fit_landmark()) must re-validate the reserved columns instead of
# trusting the config blindly.

#' Low-level constructor (no validation)
#' @noRd
new_landmark_tbl <- function(x = tibble::tibble(), config = list()) {
  out <- tibble::new_tibble(as.list(x), nrow = nrow(x), class = "landmark_tbl")
  attr(out, "landmark_config") <- config
  out
}

#' Accessor
#'
#' Returns the config list.
#' @noRd
landmark_config <- function(x) {
  attr(x, "landmark_config", exact = TRUE)
}

#' Which landmarks the summary should describe
#'
#' Once the user has removed a populated landmark (e.g. with
#' `filter(.landmark == 30)`), the summary shows only what remains.
#' @noRd
lm_landmarks_shown <- function(landmark_col, cfg){
  present <- sort(unique(landmark_col))

  built_populated <- cfg$landmarks[cfg$n_eligible > 0L]
  built_empty <- cfg$landmarks[cfg$n_eligible == 0L]

  if (all(built_populated %in% present)) {
    sort(c(present, built_empty))
  } else {
    present
  }
}

#' Per-landmark outcome counts, computed from the current rows
#'
#' Landmarks from the stored configuration are always included, so a landmark
#' with no eligible subjects shows up with zeros.
#' @noRd
lm_counts <- function(x) {
  cfg <- landmark_config(x)
  x   <- tibble::as_tibble(x)

  all_landmarks <- lm_landmarks_shown(x$.landmark, cfg)
  codes <- sort(unique(c(cfg$event_codes, x$.status[which(x$.status > 0L)])))

  # One expression per event code, e.g. n_event_2 = sum(.data$.status == 2L).
  # `!!k` inlines the code's value into the expression; `!!!` splices the
  # named list of expressions into summarise() as if typed out by hand.
  event_exprs <- purrr::map(
    rlang::set_names(codes, paste0("n_event_", codes)),
    \(k) rlang::expr(sum(.data$.status == !!k))
  )

  per_landmark <- x |>
    dplyr::summarise(
      n_at_risk  = dplyr::n(),
      n_censored = sum(.data$.status == 0L),
      !!!event_exprs,
      .by = ".landmark"
    )

  # Join onto the full landmark list so empty landmarks appear with zeros.
  tibble::tibble(.landmark = all_landmarks) |>
    dplyr::left_join(per_landmark, by = ".landmark") |>
    dplyr::mutate(dplyr::across(-".landmark", \(n) dplyr::coalesce(n, 0L))) |>
    dplyr::rename(landmark = ".landmark")
}

#' Headers for print functionality
#' @noRd
lm_header <- function(x, cfg) {
  shown <- lm_landmarks_shown(x$.landmark, cfg)
  hz <- if (is.null(cfg$horizon)) {
    "none (follow-up to observed time)"
  } else {
    paste0(cfg$horizon, " (relative to each landmark)")
  }
  lines <- c(
    paste0(
      "<landmark_tbl> ", length(shown), " landmark(s): ",
      paste(shown, collapse = ", ")
    ),
    paste0("Horizon: ", hz),
    paste0(
      "Subjects: ", dplyr::n_distinct(x$.id), " at risk in current rows (",
      cfg$n_subjects, " in source data)"
    ),
    paste0(
      "Source columns: id = ", cfg$id, ", time = ", cfg$time,
      ", status = ", cfg$status
    )
  )
  if (!is.null(cfg$status_labels)) {
    lines <- c(lines, paste0(
      "Status codes: ",
      paste0(names(cfg$status_labels), " = ", cfg$status_labels, collapse = ", ")
    ))
  }
  lines
}

#' Print a landmark dataset
#'
#' Prints a summary of the landmark configuration and per-landmark outcome
#' counts (at risk, censored, and events by cause), followed by the usual tibble
#' printout. Automatically updates after `dplyr::filter()` or `[` subsetting.
#'
#' @param x A `landmark_tbl`, as returned by [landmark_data()].
#' @param ... Passed on to the tibble print method
#'
#' @returns `x`, invisibly
#' @export
print.landmark_tbl <- function(x, ...) {
  cfg <- landmark_config(x)
  # Only show the summary when the metadata and reserved columns are intact
  # (e.g. not after select()-ing away `.status`).
  if (lm_is_intact(x)) {
    writeLines(lm_header(x, cfg))
    print(as.data.frame(lm_counts(x)), row.names = FALSE)
    cat("\n")
  }
  NextMethod()
  invisible(x)
}

#' Convert a landmark dataset to a plain tibble
#'
#' Drops the `landmark_tbl` class and its stored configuration, returning an
#' ordinary tibble with the same rows and columns.
#'
#' @param x A `landmark_tbl`, as returned by [landmark_data()].
#' @param ... Ignored.
#'
#' @returns A tibble
#' @exportS3Method tibble::as_tibble
as_tibble.landmark_tbl <- function(x, ...) {
  attr(x, "landmark_config") <- NULL
  class(x) <- c("tbl_df", "tbl", "data.frame")
  x
}

#' Tidy a landmark dataset
#'
#' Returns the per-landmark outcome counts that [print.landmark_tbl()] shows,
#' as a tibble with one row per landmark: the number at risk, the number
#' censored within the horizon, and the number of events of each type. Counts
#' are recomputed from the current rows, so they follow `dplyr::filter()` and
#' `[` subsetting.
#'
#' @param x A `landmark_tbl`, as returned by [landmark_data()].
#' @param status_labels If `TRUE` and `status` was supplied to
#'   [landmark_data()] as a factor, name the event columns after the factor
#'   levels (for example `n_relapse`, `n_death`) instead of the integer codes
#'   (`n_event_1`, `n_event_2`). Levels are made syntactic with [make.names()].
#'   When `status` was numeric there are no labels to use, so the columns keep
#'   their code names and a warning is given. Default `FALSE`.
#' @param show_proportion If `TRUE`, add `prop_at_risk`: the number at risk
#'   divided by the number of subjects in the source data. Default `FALSE`.
#' @param ... Must be empty. Present for compatibility with the generic.
#'
#' @returns A tibble with columns `landmark`, `prop_at_risk` (only when
#'   `show_proportion = TRUE`), `n_at_risk`, `n_censored`, and one event-count
#'   column per event type.
#'
#' @examples
#' subjects <- data.frame(
#'   patient_id    = 1:10,
#'   followup_days = c(20, 30, 48, 65, 85, 110, 145, 190, 250, 270),
#'   event         = c(1, 2, 1, 2, 0, 1, 2, 1, 0, 1)
#' )
#' lm_dat <- landmark_data(
#'   subjects, patient_id, followup_days, event,
#'   landmarks = c(30, 60, 90), horizon = 180
#' )
#' tidy(lm_dat)
#' tidy(lm_dat, show_proportion = TRUE)
#'
#' # With a factor status, event columns can be named after its levels
#' subjects$event <- factor(
#'   subjects$event, levels = 0:2, labels = c("censored", "relapse", "death")
#' )
#' lm_fac <- landmark_data(
#'   subjects, patient_id, followup_days, event,
#'   landmarks = c(30, 60, 90), horizon = 180
#' )
#' tidy(lm_fac, status_labels = TRUE)
#' @exportS3Method generics::tidy
tidy.landmark_tbl <- function(x,
                              status_labels = FALSE,
                              show_proportion = FALSE,
                              ...) {
  call <- rlang::current_env()
  rlang::check_dots_empty(call = call)
  lm_check_flag(status_labels, "status_labels", call)
  lm_check_flag(show_proportion, "show_proportion", call)
  # Unlike print(), tidy() errors when the object is broken: returning a
  # half-built table would be worse than failing loudly.
  lm_check_intact(x, call = call)
  cfg <- landmark_config(x)

  out <- lm_counts(x)

  if (status_labels) {
    if (is.null(cfg$status_labels)) {
      cli::cli_warn(
        c(
          "{.code status_labels = TRUE} has no effect because {.arg status} was numeric.",
          "i" = "Supply {.arg status} to {.fn landmark_data} as a factor to label the event columns."
        )
      )
    } else {
      # cfg$status_labels is named by code ("0", "1", ...), so index by code
      # rather than by position.
      old_names <- paste0("n_event_", cfg$event_codes)
      new_names <- paste0(
        "n_", make.names(cfg$status_labels[as.character(cfg$event_codes)])
      )
      lm_check_label_names(new_names, old_names, names(out), call)
      # rename() takes new_name = old_name pairs; set_names() builds them.
      out <- dplyr::rename(out, !!!rlang::set_names(old_names, new_names))
    }
  }

  if (show_proportion) {
    out <- out |>
      dplyr::mutate(prop_at_risk = .data$n_at_risk / cfg$n_subjects) |>
      dplyr::relocate("prop_at_risk", .after = "landmark")
  }

  out
}
