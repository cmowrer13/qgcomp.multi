# ------------------------------------------------------------------------------
# Diagnostic helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
build_qgcompmulti_bootstrap <- function(coef_draws,
                                        B_requested,
                                        failure_log = NULL) {
  if (is.null(coef_draws)) {
    B_success <- 0L
  } else {
    B_success <- nrow(as.matrix(coef_draws))
  }
  B_failed <- max(0L, as.integer(B_requested) - B_success)
  list(
    coef_draws = coef_draws,
    B_requested = as.integer(B_requested),
    B_success = B_success,
    B_failed = B_failed,
    failure_log = failure_log
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_validate_bootstrap <- function(bootstrap, labels) {
  coef_names <- labels$coefficient_names
  if (!is.null(bootstrap$coef_draws)) {
    if (!is.matrix(bootstrap$coef_draws) && !is.data.frame(bootstrap$coef_draws)) {
      stop("`bootstrap$coef_draws` must be `NULL`, a matrix, or a data frame.", call. = FALSE)
    }
    coef_draws <- as.matrix(bootstrap$coef_draws)
    if (is.null(colnames(coef_draws))) {
      stop("`bootstrap$coef_draws` must have column names.", call. = FALSE)
    }
    if (!is.null(coef_names) && !identical(colnames(coef_draws), coef_names)) {
      stop(
        "Column names of `bootstrap$coef_draws` must match `labels$coefficient_names`.",
        call. = FALSE
      )
    }
  }
  scalar_fields <- c("B_requested", "B_success", "B_failed")
  for (field in scalar_fields) {
    value <- bootstrap[[field]]
    if (!is.null(value) && length(value) != 1L) {
      stop(sprintf("`bootstrap$%s` must be length 1.", field), call. = FALSE)
    }
  }
  if (!is.null(bootstrap$coef_draws) &&
      !is.null(bootstrap$B_success) &&
      nrow(as.matrix(bootstrap$coef_draws)) != bootstrap$B_success) {
    stop("`bootstrap$B_success` must match the number of retained bootstrap draws.", call. = FALSE)
  }
  if (!is.null(bootstrap$B_requested) &&
      !is.null(bootstrap$B_success) &&
      !is.null(bootstrap$B_failed) &&
      bootstrap$B_failed != (bootstrap$B_requested - bootstrap$B_success)) {
    stop("`bootstrap$B_failed` must equal `B_requested - B_success`.", call. = FALSE)
  }
  if (!is.null(bootstrap$failure_log) && !is.data.frame(bootstrap$failure_log)) {
    stop("`bootstrap$failure_log` must be `NULL` or a data frame.", call. = FALSE)
  }
  invisible(bootstrap)
}

#' @keywords internal
#' @noRd
qgcompmulti_support_summary <- function(object) {
  validate_qgcompmulti(object)
  intervention_grid <- object$prediction$intervention_grid
  msm_grid <- object$prediction$msm_grid
  list(
    psi1_values = sort(unique(intervention_grid$psi1)),
    psi2_values = sort(unique(intervention_grid$psi2)),
    intervention_psi1_range = range(intervention_grid$psi1, na.rm = TRUE),
    intervention_psi2_range = range(intervention_grid$psi2, na.rm = TRUE),
    msm_psi1_range = range(msm_grid$psi1, na.rm = TRUE),
    msm_psi2_range = range(msm_grid$psi2, na.rm = TRUE),
    n_grid_points = nrow(intervention_grid)
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_support_flags <- function(object, support_summary) {
  is_quantized <- !is.null(object$mixtures$q)
  is_original_scale <- is.null(object$mixtures$q)
  pooled_percentile_grid <- is_original_scale
  list(
    is_quantized = is_quantized,
    is_original_scale = is_original_scale,
    pooled_percentile_grid = pooled_percentile_grid,
    centered_msm_grid = is_original_scale && identical(object$mixtures$centering, "median")
  )
}
#' @keywords internal
#' @noRd
build_qgcompmulti_support_diagnostic <- function(object) {
  validate_qgcompmulti(object)
  support_summary <- qgcompmulti_support_summary(object)
  flags <- qgcompmulti_support_flags(object, support_summary)
  notes <- if (isTRUE(flags$pooled_percentile_grid)) {
    c(
      "Intervention support is defined by pooled marginal percentile values for each mixture.",
      "For original-scale fits, these grid values may be less comparable when the two mixtures live on very different scales."
    )
  } else {
    "Intervention support is the quantized intervention grid used in fitting."
  }
  structure(
    list(
      diagnostic_type = "support",
      mode = if (flags$is_quantized) "quantized" else "original_scale",
      centering = object$mixtures$centering,
      intervention_grid = object$prediction$intervention_grid,
      msm_grid = object$prediction$msm_grid,
      support_summary = support_summary,
      flags = flags,
      notes = notes
    ),
    class = "qgcompmulti_support_diagnostic"
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_failure_summary <- function(failure_log) {
  if (is.null(failure_log) || nrow(failure_log) == 0L) {
    return(data.frame(message = character(), n = integer(), row.names = NULL))
  }
  if (!"message" %in% names(failure_log)) {
    return(data.frame(message = character(), n = integer(), row.names = NULL))
  }
  counts <- sort(table(failure_log$message), decreasing = TRUE)
  data.frame(
    message = names(counts),
    n = as.integer(counts),
    row.names = NULL,
    check.names = FALSE
  )
}
#' @keywords internal
#' @noRd
build_qgcompmulti_bootstrap_diagnostic <- function(object) {
  validate_qgcompmulti(object)
  B_requested <- object$bootstrap$B_requested
  B_success <- object$bootstrap$B_success
  B_failed <- object$bootstrap$B_failed
  success_rate <- if (is.null(B_requested) || B_requested == 0L) {
    NA_real_
  } else {
    B_success / B_requested
  }
  failure_summary <- qgcompmulti_failure_summary(object$bootstrap$failure_log)
  flags <- list(
    any_failures = !is.null(B_failed) && B_failed > 0L,
    high_failure_rate = is.finite(success_rate) && success_rate < 0.9
  )
  structure(
    list(
      diagnostic_type = "bootstrap",
      B_requested = B_requested,
      B_success = B_success,
      B_failed = B_failed,
      success_rate = success_rate,
      failure_log = object$bootstrap$failure_log,
      failure_summary = failure_summary,
      flags = flags
    ),
    class = "qgcompmulti_bootstrap_diagnostic"
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_adequacy_metrics <- function(surface_comparison) {
  residual <- surface_comparison$residual
  exact_mean <- surface_comparison$exact_mean
  msm_mean <- surface_comparison$msm_mean
  list(
    n_grid_points = nrow(surface_comparison),
    mae = mean(abs(residual)),
    rmse = sqrt(mean(residual^2)),
    max_abs_error = max(abs(residual)),
    mean_signed_error = mean(residual),
    correlation = stats::cor(exact_mean, msm_mean)
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_adequacy_flags <- function(metrics) {
  list(
    large_max_error = is.finite(metrics$max_abs_error) && metrics$max_abs_error > 1,
    nontrivial_average_misfit = is.finite(metrics$mae) && metrics$mae > 0.25
  )
}
#' @keywords internal
#' @noRd
build_qgcompmulti_adequacy_diagnostic <- function(object) {
  validate_qgcompmulti(object)
  comparison <- object$prediction$surface_comparison
  metrics <- qgcompmulti_adequacy_metrics(comparison)
  flags <- qgcompmulti_adequacy_flags(metrics)
  notes <- "Adequacy compares the exact fit-time counterfactual surface to the fitted MSM surface on the response scale."
  structure(
    list(
      diagnostic_type = "adequacy",
      comparison = comparison,
      summary_metrics = metrics,
      flags = flags,
      notes = notes
    ),
    class = "qgcompmulti_adequacy_diagnostic"
  )
}
