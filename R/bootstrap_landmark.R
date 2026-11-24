#' Two-stage bootstrap procedure for propagating uncertainty from LME predictions
#'
#' @description
#' `bootstrap_landmark()` performs a subject-level bootstrap of the complete
#' landmark modeling workflow. For each bootstrap sample, landmark risk sets are
#' re-sampled, longitudinal mixed-effects models are refit, predictions at each
#' landmark time are regenerated, and landmark survival models are refit using
#' the updated predictions.
#'
#' The resulting bootstrap distribution is used to quantify uncertainty arising
#' from estimation of the longitudinal trajectory and its propagation into the
#' landmark survival model. Bootstrap parameter estimates are presented as the
#' median of the bootstrap sample distribution and percentile confidence
#' intervals and bootstrap p-values are additionally reported.
#'
#' This approach is particularly useful when landmark covariates are derived
#' from predicted values obtained from mixed-effects models, where standard
#' errors from the survival model alone do not account for uncertainty in the
#' first-stage prediction process.
#'
#' @param baseline A `landmark_tbl`, as returned by [landmark_data()].
#' @param pipeline A function accepting exactly one argument, the `baseline`
#'  parameter. This function should define the rest of the landmark workflow,
#'  including calls to [add_landmark_covariates()] and [fit_landmark()] with
#'  the resulting output being a `landmark_fit` object.
#' @param n_boot Number of bootstrap samples to take.
#' @param exponentiate f `TRUE`, report hazard ratios (subdistribution hazard
#'   ratios for Fine-Gray models) instead of log hazard ratios.
#' @param conf.level Confidence level for the interval. Default `0.95`.
#' @param alternative Alternative hypothesis for bootstrap p-values. One of:
#'  "two.sided" (default), "greater", or "less".
#'
#' @returns A tibble structured similarly to the output of
#'  [tidy.landmark_fit()]. Estimates from the complete data are stored in
#'  columns with "__full" appended while bootstrap estimates have suffix
#'  "__boot". Estimates and diagnostics at the bootstrap-sample level are stored
#'  as lists.
#'
#' @seealso [landmark_data()], [add_landmark_covariates()], [fit_landmark()],
#'  [tidy.landmark_fit()]
#' @export
#'
#' @examples
#' # PBC: status 0 = censored, 1 = transplant, 2 = death
#' pbcseq_baseline |>
#'
#'   # Set up baseline landmark dataset with landmark_data()
#'   landmark_data(id, futime, status, landmarks = c(365, 730), horizon = 1825) |>
#'
#'   # Pipe to bootstrap function
#'   bootstrap_landmark(
#'
#'     # Continue landmark procedure pipeline
#'     pipeline = \(data){
#'       data |>
#'         add_landmark_covariates(pbcseq_visits, id, day, lme = c(bili, albumin)) |>
#'         fit_landmark(~ age + log(bili) + albumin, cause = 2)
#'     },
#'
#'     n_boot = 200,
#'     exponentiate = FALSE,
#'     conf.level = 0.95,
#'     alternative = 'two.sided'
#'   )

bootstrap_landmark <- function(
  baseline,
  pipeline,
  n_boot = 200,
  exponentiate = FALSE,
  conf.level = 0.95,
  alternative = 'two.sided'
  ){

  # Validate input
  lm_check_flag(exponentiate, "exponentiate", call)
  lm_check_conf_level(conf.level, call = call)
  rlang::check_string(alternative)
  alternative <- rlang::arg_match(
    alternative, c("two.sided", "greater", "less")
  )

  # Determine null based on exponentiate argument
  null <- ifelse(exponentiate, 1, 0)

  # Get fit in complete data
  fit_full <- baseline |>

    # Call user-defined pipeline on data
    pipeline() |>

    # Use tidy.landmark_fit to get model fits
    tidy(
      exponentiate = exponentiate,
      conf.int = TRUE,
      conf.level = conf.level
    ) |>

    # Select relevant variables
    dplyr::select(

      # Fit-invariant measures
      .data$.landmark,
      .data$.model,
      .data$.cause,
      .data$origin,
      .data$term,

      # Full fit measures
      .data$estimate,
      .data$conf.low,
      .data$conf.high,
      .data$p.value,
      .data$.fit_status,
      .data$.message
    ) |>

    dplyr::rename_with(
      .cols = c(
        .data$estimate,
        .data$conf.low,
        .data$conf.high,
        .data$p.value,
        .data$.fit_status,
        .data$.message
      ),
      .f = \(col) paste0(col, "__full")
    )

  # Iterate and get bootstrap samples and fits
  fit_boot <- purrr::map(

    # Iterate over bootstrap indices
    1:n_boot,

    function(b){

      # From baseline dataset, pull IDs at each landmark and update baseline dataset
      id_sample <-

        # Sample available IDs at each landmark
        purrr::map(
          attributes(baseline)$landmark_config$id_list,
          \(ids){
            sample(
              x = ids,
              size = length(ids),
              replace = TRUE
            )
          }
        ) |>

        # Create tibble of sampled IDs at landmark
        tibble::tibble(
          .id = _,
          .landmark = attributes(baseline)$landmark_config$landmarks,
        ) |>
        tidyr::unnest(cols = .data$.id) |>

        # Create new bootstrap ID (so  latter LME does not infer incorrect grouping)
        dplyr::mutate(.boot_id = 1:dplyr::n()) |>

        # Join to formatted landmark_tbl
        dplyr::left_join(
          baseline,
          by = c(".id", ".landmark"),
          relationship = "many-to-one"
        ) |>

        # Reset landmark_tbl class
        new_landmark_tbl(
          config = attributes(baseline)$landmark_config
        )

      # Set boot attribute to detect in next function
      attributes(id_sample)$boot <- TRUE

      id_sample |>

        # Call user-defined pipeline on data
        pipeline() |>

        # Use tidy.landmark_fit to get model fits
        tidy(exponentiate = exponentiate) |>

        # Select relevant variables
        dplyr::select(
          .data$.landmark,
          .data$term,
          .data$estimate,
          .data$.fit_status,
          .data$.message
        )

    } |>

      # Do not print messages to console for each bootstrap fit
      suppressMessages()

  ) |>

    # Join bootstrap fits and aggregate
    dplyr::bind_rows() |>

    # Summarize
    dplyr::summarize(
      estimate__boot = stats::median(.data$estimate, na.rm = TRUE),
      conf.low__boot = stats::quantile(.data$estimate, probs = (1 - conf.level) / 2, na.rm = TRUE),
      conf.high__boot = stats::quantile(.data$estimate, probs = 1 - (1 - conf.level) / 2, na.rm = TRUE),
      p.value__boot = bootstrap_p(.data$estimate, null = null, alternative = alternative),
      .estimates__boot = list(.data$estimate),
      .fit_status__boot = list(.data$.fit_status),
      .message__boot = list(.data$.message),
      .by = c(.data$.landmark, .data$term)
    )

  # Join full-data and bootstrap fits and return
  fit_full |>
    dplyr::left_join(
      fit_boot,
      by = c('.landmark', 'term'),
      relationship = 'one-to-one'
    ) |>

    # Intuitive ordering of columns
    dplyr::select(
      .data$.landmark,
      .data$.model,
      .data$.cause,
      .data$origin,
      .data$term,
      .data$estimate__full,
      .data$conf.low__full,
      .data$conf.high__full,
      .data$p.value__full,
      .data$estimate__boot,
      .data$conf.low__boot,
      .data$conf.high__boot,
      .data$p.value__boot,
      .data$.fit_status__full,
      .data$.message__full,
      .data$.estimates__boot,
      .data$.fit_status__boot,
      .data$.message__boot
    )

}

#' Function to calculate percentile-based bootstrap p-value
#'
#' `estimate` is the vector of parameter estimates from the bootstrap samples.
#' `null` defines the location of the null hypothesis and is automatically
#' determined based on if the estimates are presented on the log scale.
#' `alternative` defines the alternative hypothesis for the test. One of:
#'  "two.sided" (default), "greater", or "less".
#'
#' Returns the numeric p-value.
#' @noRd
bootstrap_p <- function(
    estimate,
    null,
    alternative
) {

  p_left <- mean(estimate <= null, na.rm = TRUE)
  p_right <- mean(estimate >= null, na.rm = TRUE)

  dplyr::case_when(
    alternative == "greater" ~ p_left,
    alternative == "less" ~ p_right,
    alternative == "two.sided" ~ 2 * min(p_left, p_right)
  ) |>
    min(1)
}
