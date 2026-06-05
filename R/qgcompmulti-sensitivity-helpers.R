# ------------------------------------------------------------------------------
# Sensitivity helpers
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
qgcompmulti_sensitivity_fit_args <- function(fit) {
  validate_qgcompmulti(fit)
  list(
    f = fit$formula,
    data = fit$fits$outcome_fit$data,
    mix1 = fit$mixtures$mix1,
    mix2 = fit$mixtures$mix2,
    interaction = fit$analysis$interaction,
    family = fit$analysis$family,
    q = fit$mixtures$q,
    centering = fit$mixtures$centering,
    B = fit$analysis$B,
    id = fit$analysis$id,
    MCsize = fit$analysis$MCsize,
    seed = fit$analysis$seed
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_extract_fit_snapshot <- function(fit) {
  validate_qgcompmulti(fit)
  coef_vec <- coef(fit)
  adequacy_diag <- build_qgcompmulti_adequacy_diagnostic(fit)
  bootstrap_diag <- build_qgcompmulti_bootstrap_diagnostic(fit)
  out <- list(
    coefficients = coef_vec,
    std_error = fit$results$std_error,
    adequacy = adequacy_diag$summary_metrics,
    bootstrap = list(
      B_requested = bootstrap_diag$B_requested,
      B_success = bootstrap_diag$B_success,
      B_failed = bootstrap_diag$B_failed,
      success_rate = bootstrap_diag$success_rate
    )
  )
  if (!is.null(fit$prediction$msm_surface)) {
    out$surface_range <- diff(range(fit$prediction$msm_surface$msm_mean))
  } else {
    out$surface_range <- NA_real_
  }
  out
}
#' @keywords internal
#' @noRd
qgcompmulti_mcsize_values <- function(MCsize_values, n) {
  if (!is.numeric(MCsize_values) || length(MCsize_values) == 0L || anyNA(MCsize_values)) {
    stop("`MCsize_values` must be a non-empty numeric vector.", call. = FALSE)
  }
  MCsize_values <- as.integer(MCsize_values)
  if (any(MCsize_values <= 0L)) {
    stop("`MCsize_values` must contain positive integers.", call. = FALSE)
  }
  unique(pmin(MCsize_values, as.integer(n)))
}
#' @keywords internal
#' @noRd
qgcompmulti_q_values <- function(q_values) {
  if (!is.numeric(q_values) || length(q_values) == 0L || anyNA(q_values)) {
    stop("`q_values` must be a non-empty numeric vector.", call. = FALSE)
  }
  q_values <- as.integer(q_values)
  if (any(q_values < 2L)) {
    stop("All `q_values` must be integers greater than or equal to 2.", call. = FALSE)
  }
  unique(q_values)
}
#' @keywords internal
#' @noRd
qgcompmulti_q_comparability_note <- function() {
  paste(
    "Raw MSM coefficients are not directly comparable across different choices of q.",
    "A larger q implies a smaller one-quantile intervention step, so smaller coefficient magnitudes may be expected mechanically rather than indicating a weaker overall mixture effect."
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_build_sensitivity_table <- function(values,
                                                parameter_name,
                                                fits) {
  rows <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    snap <- qgcompmulti_extract_fit_snapshot(fit)
    coef_vec <- snap$coefficients
    se_vec <- snap$std_error
    row <- data.frame(
      value = values[[i]],
      psi1 = unname(coef_vec[["psi1"]]),
      psi2 = unname(coef_vec[["psi2"]]),
      psi1_se = unname(se_vec[["psi1"]]),
      psi2_se = unname(se_vec[["psi2"]]),
      adequacy_mae = snap$adequacy$mae,
      adequacy_rmse = snap$adequacy$rmse,
      adequacy_max_abs_error = snap$adequacy$max_abs_error,
      bootstrap_success = snap$bootstrap$B_success,
      bootstrap_failed = snap$bootstrap$B_failed,
      bootstrap_success_rate = snap$bootstrap$success_rate,
      surface_range = snap$surface_range,
      row.names = NULL
    )
    if ("psi1:psi2" %in% names(coef_vec)) {
      row$psi12 <- unname(coef_vec[["psi1:psi2"]])
      row$psi12_se <- unname(se_vec[["psi1:psi2"]])
    }
    row
  })
  out <- do.call(rbind, rows)
  names(out)[1] <- parameter_name
  row.names(out) <- NULL
  out
}
#' @keywords internal
#' @noRd
qgcompmulti_print_sensitivity_table <- function(results_table) {
  print(results_table, row.names = FALSE)
  invisible(results_table)
}
