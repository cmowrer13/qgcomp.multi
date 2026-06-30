#' Summarize a qgcompmulti fit
#'
#' Builds a structured summary object for a fitted `qgcompmulti` model. The
#' summary foregrounds the MSM results because those are the main reported
#' estimands, but those coefficients should still be interpreted in light of the
#' intervention coding and the adequacy of the MSM as a summary of the fitted
#' outcome-model surface.
#'
#' @param object A fitted `qgcompmulti` object.
#' @param ... Not used.
#'
#' @return An object of class `"summary.qgcompmulti"` with components for the
#' matched call, original formula, model overview, mixture definitions,
#' labeled MSM coefficient table, outcome-model context, and stored labels.
#'
#' @seealso [print.qgcompmulti()], [coef.qgcompmulti()],
#'   [confint.qgcompmulti()], [qgcomp.glm.multi()]
#' @export

summary.qgcompmulti <- function(object, ...) {
  validate_qgcompmulti(object)

  structure(
    build_qgcompmulti_summary(object),
    class = "summary.qgcompmulti"
  )
}

#' Print a summary.qgcompmulti object
#'
#' Prints a structured summary of a fitted `qgcompmulti` model, emphasizing the
#' marginal structural model results while reporting the outcome model as
#' supporting context.
#'
#' @param x A `"summary.qgcompmulti"` object produced by
#' [summary.qgcompmulti()].
#' @param ... Not used.
#'
#' @return The input object, invisibly.
#'
#' @details
#' The printed summary highlights the marginal structural model (MSM)
#' coefficient results while reporting the outcome-model family, formula,
#' exposure-mode setting, sample-size metadata, and compact outcome-model fit
#' statistics as supporting context.
#' @rdname summary.qgcompmulti
#' @export

print.summary.qgcompmulti <- function(x, ...) {
  validate_summary_qgcompmulti(x)

  fit_overview <- x$fit_overview
  outcome_info <- x$outcome_model_info

  cat("Summary of qgcompmulti fit\n\n")
  qgcompmulti_print_call(x$call)

  cat("\nModel overview:\n")
  cat("  Formula: ", qgcompmulti_deparse_one_line(x$formula), "\n", sep = "")
  cat("  Outcome: ", fit_overview$outcome, "\n", sep = "")
  cat("  Family: ", fit_overview$family, "\n", sep = "")
  cat(
    "  Estimand: ",
    fit_overview$estimand,
    if (isTRUE(fit_overview$estimand_scale_defaulted)) " (default)" else "",
    "
",
sep = ""
  )
  cat("  MSM fitting scale: ", fit_overview$msm_fitting_scale, "
", sep = "")
  cat("  Default interval method: ", fit_overview$default_interval_method, "
", sep = "")
  cat("  Observations used: ", fit_overview$n_used, "\n", sep = "")
  if (!is.null(fit_overview$n_input) && fit_overview$n_input != fit_overview$n_used) {
    cat("  Input observations: ", fit_overview$n_input, "\n", sep = "")
  }
  cat("  Exposure mode: ", fit_overview$exposure_mode, "\n", sep = "")
  cat("  ", fit_overview$interaction, "\n", sep = "")
  if (!is.null(fit_overview$bootstrap_requested)) {
    if (!is.null(fit_overview$bootstrap_success) &&
        fit_overview$bootstrap_success != fit_overview$bootstrap_requested) {
      cat(
        "  Bootstrap replications: ",
        fit_overview$bootstrap_requested,
        " requested, ",
        fit_overview$bootstrap_success,
        " retained\n",
        sep = ""
      )
    } else {
      cat(
        "  Bootstrap replications: ",
        fit_overview$bootstrap_requested,
        "\n",
        sep = ""
      )
    }
  }
  cat("  Monte Carlo size: ", fit_overview$MCsize, "\n", sep = "")
  if (!is.null(fit_overview$seed)) {
    cat("  Random seed: ", fit_overview$seed, "\n", sep = "")
  }
  if (isTRUE(fit_overview$has_clusters)) {
    cat(
      "  Clusters: ",
      fit_overview$n_clusters,
      " (id = ",
      fit_overview$cluster_var,
      ")\n",
      sep = ""
    )
  }

  cat("\n")
  qgcompmulti_print_mixtures(x$mixtures)

  cat("\n")
  qgcompmulti_print_msm_table(x$msm_table)

  cat("\nOutcome model context:\n")
  cat("  Model class: ", outcome_info$model_class, "\n", sep = "")
  cat("  Estimated parameters: ", outcome_info$n_parameters, "\n", sep = "")
  cat("  AIC: ", qgcompmulti_format_metric(outcome_info$aic), "\n", sep = "")
  cat(
    "  Null deviance: ",
    qgcompmulti_format_metric(outcome_info$null_deviance),
    "\n",
    sep = ""
  )
  cat(
    "  Residual deviance: ",
    qgcompmulti_format_metric(outcome_info$residual_deviance),
    "\n",
    sep = ""
  )

  invisible(x)
}
