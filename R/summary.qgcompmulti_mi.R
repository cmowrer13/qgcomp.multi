#' Summarize a pooled qgcompmulti multiple-imputation fit
#'
#' Builds a structured summary object for a fitted `qgcompmulti_mi` model. The
#' summary follows the same general presentation as `summary.qgcompmulti()`
#' while adding the pooled multiple-imputation metadata and Rubin diagnostics
#' that matter for downstream reporting.
#'
#' @param object A fitted `qgcompmulti_mi` object.
#' @param ... Not used.
#'
#' @return An object of class `"summary.qgcompmulti_mi"` with components for
#' the matched call, original formula, model overview, MI overview, mixture
#' definitions, labeled pooled MSM coefficient table, pooled Rubin diagnostics,
#' and stored labels.
#'
#' @seealso [print.qgcompmulti_mi()], [coef.qgcompmulti_mi()],
#'   [confint.qgcompmulti_mi()], [qgcomp.glm.multi.mi()]
#' @export
summary.qgcompmulti_mi <- function(object, ...) {
  validate_qgcompmulti_mi(object)

  structure(
    build_qgcompmulti_mi_summary(object),
    class = "summary.qgcompmulti_mi"
  )
}

#' Print a summary.qgcompmulti_mi object
#'
#' Prints a structured summary of a pooled `qgcompmulti_mi` fit, emphasizing
#' the pooled MSM results while reporting the shared analysis settings and
#' compact Rubin diagnostics as supporting context.
#'
#' @param x A `"summary.qgcompmulti_mi"` object produced by
#'   [summary.qgcompmulti_mi()].
#' @param ... Not used.
#'
#' @return The input object, invisibly.
#'
#' @rdname summary.qgcompmulti_mi
#' @export
print.summary.qgcompmulti_mi <- function(x, ...) {
  validate_summary_qgcompmulti_mi(x)

  fit_overview <- x$fit_overview
  mi_overview <- x$mi_overview

  cat("Summary of qgcompmulti multiple-imputation fit\n\n")
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
  cat("  Interval method: wald (pooled multiple imputation)
", sep = "")
  cat("  Observations used: ", fit_overview$n_used, "\n", sep = "")
  if (!is.null(fit_overview$n_input) && fit_overview$n_input != fit_overview$n_used) {
    cat("  Input observations: ", fit_overview$n_input, "\n", sep = "")
  }
  cat("  Exposure mode: ", fit_overview$exposure_mode, "\n", sep = "")
  cat("  ", fit_overview$interaction, "\n", sep = "")
  cat("  Bootstrap replications per imputation: ", fit_overview$bootstrap_requested, "\n", sep = "")
  cat("  Monte Carlo size: ", fit_overview$MCsize, "\n", sep = "")
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

  cat("\nMultiple imputation overview:\n")
  cat("  Imputed datasets: ", mi_overview$m, "\n", sep = "")
  cat("  Input type: ", mi_overview$input_type, "\n", sep = "")
  cat("  Retained per-imputation fits: ", if (isTRUE(mi_overview$keep_fits)) "yes" else "no", "\n", sep = "")
  if (!is.null(mi_overview$seed)) {
    cat("  Master seed: ", mi_overview$seed, "\n", sep = "")
  }
  cat("  Stored fit-specific seeds: ", mi_overview$fit_seed_count, "\n", sep = "")

  cat("\n")
  qgcompmulti_print_mixtures(x$mixtures)

  cat("\n")
  qgcompmulti_print_msm_table(x$msm_table)

  cat("\n")
  qgcompmulti_print_pooling_table(x$pooling_table)

  invisible(x)
}
