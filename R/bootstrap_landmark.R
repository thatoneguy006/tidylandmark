#' Two-stage bootstrap for propagating uncertainty from LME predictions
#'
#' @description
#' `bootstrap_landmark()` bootstraps the complete landmark workflow. In each
#' replicate, the risk set at every landmark is resampled with replacement,
#' the user's `pipeline` is re-run on the resampled data (refitting any
#' longitudinal mixed-effects models, regenerating their predictions at each
#' landmark and refitting the landmark survival models), and the coefficients
#' are collected.
#'
#' The bootstrap distribution quantifies uncertainty from estimating the
#' longitudinal trajectories and its propagation into the landmark survival
#' models, which the standard errors from [fit_landmark()] ignore. Results are
#' reported as percentile confidence intervals, percentile bootstrap p-values
#' and the median of the bootstrap distribution, next to the full-data fit.
#'
#' @details
#' ## Resampling
#' Rows are resampled independently within each landmark, which matches
#' [fit_landmark()] fitting one independent model per landmark and
#' [add_landmark_covariates()] fitting one mixed model per landmark from that
#' landmark's risk set. The bootstrap therefore supports per-landmark
#' inference, but not comparisons between landmarks. A subject drawn more than
#' once is treated as distinct subjects, each with its own random effects in
#' the mixed models.
#'
#' Risk sets are read from the rows of `baseline` as they are now, so a
#' `landmark_tbl` filtered after [landmark_data()] is resampled as filtered.
#' Create any LME predictors used by the survival model inside `pipeline` so
#' their mixed models are refitted on each draw. For spline coefficients, set
#' the interior and boundary knots before resampling and pass them to the
#' spline function. The bootstrap checks the fitted basis at each landmark.
#'
#' ## Failed replicates
#' Each replicate runs quietly: progress bars (when `pipeline` uses the default
#' `progress`), messages and warnings are suppressed, and a replicate whose
#' `pipeline` errors is recorded as `"failed"` at every landmark instead of
#' stopping the bootstrap. `n_boot_ok` gives the number of replicates with an
#' estimate for each term; summaries use only those. Check
#' `.fit_status__boot` when `n_boot_ok` is well below `n_boot`. An unavailable
#' coefficient (for example, a comparison to an absent factor level) has an
#' `NA` estimate in that draw. Intervals and p-values require at least two
#' usable estimates.
#'
#' ## Reproducibility and p-values
#' Results depend on the random number generator; call [set.seed()] first for
#' reproducible results. A zero p-value means that no usable draw fell in the
#' relevant tail. Its Monte Carlo precision depends on `n_boot_ok`.
#'
#' @param baseline A `landmark_tbl`, as returned by [landmark_data()].
#' @param pipeline A function of one argument, a `landmark_tbl`, that runs the
#'   rest of the workflow (e.g. [add_landmark_covariates()] then
#'   [fit_landmark()]) and returns a `landmark_fit`.
#' @param n_boot Number of bootstrap replicates. Default `200`.
#' @param exponentiate If `TRUE`, report hazard ratios (subdistribution hazard
#'   ratios for Fine-Gray models) instead of log hazard ratios.
#' @param conf.level Confidence level for the intervals. Default `0.95`.
#' @param alternative Alternative hypothesis for the bootstrap p-values: one of
#'   `"two.sided"` (default), `"greater"` or `"less"`.
#' @param progress Show a progress bar over replicates? Defaults to the
#'   `tidylandmark.progress` option, or [interactive()] if it is unset.
#'
#' @returns A tibble with one row per landmark and term of the full-data fit:
#'   * `.landmark`, `.model`, `.cause`, `origin`, `term`: as in
#'     [tidy.landmark_fit()].
#'   * `estimate__full`, `conf.low__full`, `conf.high__full`,
#'     `p.value__full`: the full-data fit with model-based (Wald) inference.
#'   * `estimate__boot`, `conf.low__boot`, `conf.high__boot`,
#'     `p.value__boot`: median, percentile interval and percentile p-value of
#'     the bootstrap distribution.
#'   * `n_boot_ok`: replicates with an estimate for this term.
#'   * `.fit_status__full`, `.message__full`: fit status of the full-data fit.
#'   * `.estimates__boot`: one estimate per replicate, in draw order.
#'   * `.fit_status__boot`, `.message__boot`: list-columns with the fit status
#'     and message for the same replicate, landmark and term.
#'
#' @seealso [landmark_data()], [add_landmark_covariates()], [fit_landmark()],
#'   [tidy.landmark_fit()]
#' @export
#'
#' @examples
#' \donttest{
#' # PBC: status 0 = censored, 1 = transplant, 2 = death
#' set.seed(2026)
#' pbcseq_baseline |>
#'   landmark_data(id, futime, status, landmarks = c(365, 730), horizon = 1825) |>
#'   bootstrap_landmark(
#'     pipeline = \(data) {
#'       data |>
#'         add_landmark_covariates(pbcseq_visits, id, day, lme = c(bili, albumin)) |>
#'         fit_landmark(~ age + log(bili) + albumin, cause = 2)
#'     },
#'     n_boot = 50
#'   )
#' }
bootstrap_landmark <- function(baseline,
                               pipeline,
                               n_boot = 200,
                               exponentiate = FALSE,
                               conf.level = 0.95,
                               alternative = c("two.sided", "greater", "less"),
                               progress = getOption("tidylandmark.progress", interactive())) {
  call <- rlang::current_env()

  # ---- 1. Validate ----------------------------------------------------------
  lm_check_intact(baseline, "baseline", call)
  if (!is.function(pipeline)) {
    cli::cli_abort(
      "{.arg pipeline} must be a function taking one argument, not {.cls {class(pipeline)}}.",
      call = call
    )
  }
  if (!rlang::is_scalar_integerish(n_boot, finite = TRUE) || n_boot < 2) {
    cli::cli_abort("{.arg n_boot} must be a single whole number of at least 2.", call = call)
  }
  lm_check_flag(exponentiate, "exponentiate", call)
  lm_check_conf_level(conf.level, call = call)
  alternative <- rlang::arg_match(alternative, error_call = call)
  lm_check_flag(progress, "progress", call)

  # ---- 2. Full-data fit -----------------------------------------------------
  fit_full <- pipeline(baseline)
  if (!inherits(fit_full, "landmark_fit")) {
    cli::cli_abort(
      "{.arg pipeline} must return a {.cls landmark_fit}, not {.cls {class(fit_full)}}.",
      call = call
    )
  }
  original <- lm_present_covariates(baseline)
  fitted <- lm_present_covariates(fit_full$data)
  used <- intersect(original$column[original$method == "lme"], fit_full$config$lme_columns)
  inherited <- used[vapply(used, \(col) {
    identical(
      original$details[[match(col, original$column)]],
      fitted$details[[match(col, fitted$column)]]
    )
  }, logical(1))]
  if (length(inherited)) {
    cli::cli_abort(c(
      "The survival model uses precomputed LME covariate{?s} {.field {inherited}}.",
      "i" = "Move {.fn add_landmark_covariates} inside {.arg pipeline} to refit each mixed model."
    ), call = call)
  }
  full_cols <- c("estimate", "conf.low", "conf.high", "p.value", ".fit_status", ".message")
  full <- tidy(fit_full, exponentiate = exponentiate, conf.int = TRUE, conf.level = conf.level)
  full <- full[c(".landmark", ".model", ".cause", "origin", "term", full_cols)]
  full <- dplyr::rename_with(full, \(x) paste0(x, "__full"), .cols = dplyr::all_of(full_cols))

  # ---- 3. Bootstrap replicates ----------------------------------------------
  # Row numbers of each landmark's risk set, from the current rows.
  rows_by_landmark <- split(seq_len(nrow(baseline)), baseline$.landmark)
  template <- full[c(".landmark", "term")]

  reps <- purrr::map(
    seq_len(n_boot),
    \(b) lm_boot_replicate(baseline, pipeline, rows_by_landmark, template, fit_full),
    .progress = lm_progress_arg(progress, "Bootstrapping")
  ) |>
    purrr::list_rbind()

  # ---- 4. Summarise on the log scale ----------------------------------------
  # Transform after interpolation so reporting on either scale is consistent.
  probs <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
  boot <- reps |>
    dplyr::summarise(
      n_boot_ok        = sum(!is.na(.data$estimate)),
      estimate__boot   = stats::median(.data$estimate, na.rm = TRUE),
      conf.low__boot   = if (.data$n_boot_ok >= 2L) {
        stats::quantile(.data$estimate, probs[1], na.rm = TRUE, names = FALSE)
      } else NA_real_,
      conf.high__boot  = if (.data$n_boot_ok >= 2L) {
        stats::quantile(.data$estimate, probs[2], na.rm = TRUE, names = FALSE)
      } else NA_real_,
      p.value__boot    = bootstrap_p(.data$estimate, null = 0, alternative = alternative),
      .estimates__boot = list(.data$estimate),
      .fit_status__boot = list(.data$.fit_status),
      .message__boot    = list(.data$.message),
      .by = c(".landmark", "term")
    )
  if (exponentiate) {
    boot <- dplyr::mutate(
      boot,
      dplyr::across(c("estimate__boot", "conf.low__boot", "conf.high__boot"), exp),
      .estimates__boot = purrr::map(.data$.estimates__boot, exp)
    )
  }

  # ---- 5. Assemble ----------------------------------------------------------
  out <- full |>
    dplyr::left_join(boot, by = c(".landmark", "term"), relationship = "one-to-one")

  out[c(
    ".landmark", ".model", ".cause", "origin", "term",
    "estimate__full", "conf.low__full", "conf.high__full", "p.value__full",
    "estimate__boot", "conf.low__boot", "conf.high__boot", "p.value__boot",
    "n_boot_ok",
    ".fit_status__full", ".message__full",
    ".estimates__boot", ".fit_status__boot", ".message__boot"
  )]
}

#' Run one bootstrap replicate
#'
#' Resamples rows within each landmark, runs `pipeline` quietly and returns the
#' coefficients in the full fit's row order, including failed and absent terms.
#' @noRd
lm_boot_replicate <- function(baseline, pipeline, rows_by_landmark, template, reference) {
  # r[sample.int()] rather than sample(r): with one row, sample(57, 1) would
  # draw from 1:57.
  idx <- unlist(
    lapply(rows_by_landmark, \(r) r[sample.int(length(r), replace = TRUE)]),
    use.names = FALSE
  )
  # tibble's `[` keeps the landmark_tbl class and config.
  boot_data <- baseline[idx, , drop = FALSE]

  # Fit problems are already recorded in .fit_status and .message, so
  # per-replicate progress bars, messages and warnings would only be noise.
  rlang::local_options(tidylandmark.progress = FALSE)
  res <- withCallingHandlers(
    tryCatch(
      {
        fit <- pipeline(boot_data)
        if (!inherits(fit, "landmark_fit")) {
          stop("`pipeline` did not return a landmark_fit.", call. = FALSE)
        }
        fit
      },
      error = \(e) e
    ),
    warning = \(w) invokeRestart("muffleWarning"),
    message = \(m) invokeRestart("muffleMessage")
  )

  if (inherits(res, "error")) {
    return(dplyr::mutate(template, estimate = NA_real_, .fit_status = "failed",
                         .message = conditionMessage(res)))
  }
  info <- res$model_info[c(".landmark", ".fit_status", ".message")]
  out <- dplyr::left_join(template, info, by = ".landmark", relationship = "many-to-one")
  missing <- is.na(out$.fit_status)
  out$.fit_status[missing] <- "unestimable"
  out$.message[missing] <- "Landmark absent from this replicate."
  out$estimate <- rep(NA_real_, nrow(out))
  for (s in unique(out$.landmark)) {
    i <- which(out$.landmark == s & !is.na(out$term))
    model <- res$models[[as.character(s)]]
    if (!length(i) || is.null(model)) next
    beta <- lm_boot_coefficients(model, reference$models[[as.character(s)]], s)
    out$estimate[i] <- unname(beta[out$term[i]])
    bad <- i[!is.finite(out$estimate[i])]
    if (length(bad)) {
      out$estimate[bad] <- NA_real_
      out$.fit_status[bad] <- "unestimable"
      reason <- "Full-data coefficient is not estimable in this draw."
      out$.message[bad] <- ifelse(is.na(out$.message[bad]), reason,
                                  paste(out$.message[bad], reason))
    }
  }
  out
}

# Recover coefficients in the full-data parameterization. With a missing
# reference level, coxph can return a finite coefficient for a different
# comparison, even when the coefficient name and factor levels are unchanged.
lm_boot_coefficients <- function(model, reference, landmark) {
  if (!identical(attr(model$terms, "predvars"), attr(reference$terms, "predvars"))) {
    cli::cli_abort(c(
      "The fitted formula or predictor basis changed at landmark {landmark}.",
      "i" = "Use the same formula and fix spline knots and boundary knots outside the pipeline."
    ))
  }
  beta <- stats::coef(model)
  same_coding <- identical(model$xlevels, reference$xlevels) &&
    identical(model$contrasts, reference$contrasts)
  if (same_coding && !anyNA(beta)) return(beta)

  # Express the fitted linear predictor in the reference design. Rebuilding
  # the design is only needed if the replicate used different factor coding.
  x <- stats::model.matrix(model)
  active <- !is.na(beta)
  eta <- drop(x[, active, drop = FALSE] %*% beta[active])
  if (!same_coding) {
    frame <- stats::model.frame(model)
    for (v in names(reference$xlevels)) {
      frame[[v]] <- factor(frame[[v]], levels = reference$xlevels[[v]],
                           ordered = is.ordered(reference$model[[v]]))
    }
    x <- stats::model.matrix(reference, data = frame)
  }

  # A stratum-specific additive constant is absorbed by the baseline hazard.
  # Center before checking which coefficients the observed design identifies.
  strata <- attr(x, "strata")
  rows <- seq_len(nrow(x))
  groups <- if (is.null(strata)) list(rows) else split(rows, strata, drop = TRUE)
  for (rows in groups) {
    block <- x[rows, , drop = FALSE]
    x[rows, ] <- sweep(block, 2L, colMeans(block))
    eta[rows] <- eta[rows] - mean(eta[rows])
  }
  design <- qr(x)
  result <- qr.coef(design, eta)
  rank <- design$rank
  if (rank > 0L && rank < ncol(x)) {
    # X_dependent = X_independent %*% aliases. Any nonzero row of
    # aliases makes that independent column's coefficient non-unique too.
    independent <- seq_len(rank)
    dependent <- seq.int(rank + 1L, ncol(x))
    r <- qr.R(design)
    aliases <- backsolve(r[independent, independent, drop = FALSE],
                         r[independent, dependent, drop = FALSE])
    affected <- rowSums(abs(aliases) > 1e-7) > 0L
    result[design$pivot[independent[affected]]] <- NA_real_
  }
  if (same_coding) {
    beta[is.na(result)] <- NA_real_
    return(beta)
  }
  result
}

#' Percentile bootstrap p-value
#'
#' `estimate` holds the bootstrap estimates on the log scale; `null` is the
#' value under the null hypothesis on that scale. Returns `NA` when there are
#' fewer than two estimates.
#' @noRd
bootstrap_p <- function(estimate, null = 0, alternative = "two.sided") {
  estimate <- estimate[!is.na(estimate)]
  if (length(estimate) < 2L) {
    return(NA_real_)
  }
  p_left  <- mean(estimate <= null)
  p_right <- mean(estimate >= null)
  p <- switch(
    alternative,
    greater   = p_left,
    less      = p_right,
    two.sided = 2 * min(p_left, p_right)
  )
  min(p, 1)
}
