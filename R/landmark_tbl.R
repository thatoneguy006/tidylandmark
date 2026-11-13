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
  if (!is.null(cfg) && all(lm_reserved_cols %in% names(x))) {
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

#' Tidy a landmark dataset for additional summarization
#'
#' Returns a tibble with the summary of the landmark configuration and
#' per-landmark outcome counts (at risk, censored, and events by cause).
#' Automatically updates after `dplyr::filter()` or `[` subsetting.

#'
#' @param x A `landmark_tbl`, as returned by [landmark_data()].
#' @param status_labels Should status labels be shown in tibble if `status`
#'  was provided as factor? Default FALSE.
#' @param show_proportion Indicate proportion of full cohort at risk? Default
#'  FALSE.
#' @param ... Ignored.
#'
#' @returns A tibble
#' @exportS3Method generics::tidy
tidy.landmark_tbl <- function(
    x,
    status_labels = FALSE,
    show_proportion = FALSE,
    ...
  ){
  cfg <- landmark_config(x)

  # Only tidy when the metadata and reserved columns are intact
  # (e.g. not after select()-ing away `.status`).
  if (!is.null(cfg) && all(lm_reserved_cols %in% names(x))) {
    df <- lm_counts(x) |> tibble::as_tibble()
  }

  # If status_labels should be set as column labels, update
  if (status_labels && !is.null(cfg$status_labels)) {
    rename_map <- rlang::set_names(
      paste0('n_event_', cfg$event_codes),
      paste0('n_', make.names(cfg$status_labels[-1])))
    df <- df |> dplyr::rename(!!!rename_map)
  }

  # If proportion should be shown, include
  if (show_proportion) {
    df <- df |>
      dplyr::mutate(prop_at_risk = .data$n_at_risk / cfg$n_subjects) |>
      dplyr::relocate(.data$prop_at_risk, .after = .data$landmark)
  }

  df
}

