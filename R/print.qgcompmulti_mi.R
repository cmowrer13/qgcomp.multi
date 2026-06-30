#' Print a pooled qgcompmulti multiple-imputation fit
#'
#' Prints a moderately detailed overview of a pooled `qgcompmulti_mi` object
#' for ordinary interactive use. The display mirrors the single-fit
#' `qgcompmulti` print method while adding the key multiple-imputation metadata
#' needed to interpret the pooled inference.
#'
#' @param x A fitted `qgcompmulti_mi` object.
#' @param ... Not used.
#'
#' @return The input object, invisibly.
#'
#' @seealso [summary.qgcompmulti_mi()], [coef.qgcompmulti_mi()],
#'   [confint.qgcompmulti_mi()], [qgcomp.glm.multi.mi()]
#' @export
print.qgcompmulti_mi <- function(x, ...) {
  validate_qgcompmulti_mi(x)

  msm_table <- qgcompmulti_labeled_coef_table(x$results, x$labels, x$analysis)

  cat("qgcompmulti multiple-imputation fit\n\n")
  qgcompmulti_print_call(x$call)

  cat("\nModel:\n")
  cat("  Outcome: ", x$data_info$outcome, "\n", sep = "")
  cat("  Family: ", qgcompmulti_family_label(x$analysis), "\n", sep = "")
  qgcompmulti_print_scale_info(x$analysis, mi = TRUE)
  cat("  Observations used: ", x$data_info$n_used, "\n", sep = "")
  cat(
    "  Exposure mode: ",
    qgcompmulti_exposure_mode_label(x$mixtures),
    "\n",
    sep = ""
  )
  cat("  ", qgcompmulti_interaction_label(x$analysis), "\n", sep = "")
  cat("  Bootstrap replications per imputation: ", x$analysis$B, "\n", sep = "")
  cat("  Monte Carlo size: ", x$analysis$MCsize, "\n", sep = "")
  if (isTRUE(x$data_info$has_clusters)) {
    cat(
      "  Clusters: ",
      x$data_info$n_clusters,
      " (id = ",
      x$data_info$cluster_var,
      ")\n",
      sep = ""
    )
  }

  cat("\nMultiple imputation:\n")
  cat("  Imputed datasets: ", x$mi_info$m, "\n", sep = "")
  cat("  Input type: ", qgcompmulti_mi_input_label(x$mi_info$input_type), "\n", sep = "")
  cat("  Retained per-imputation fits: ", if (isTRUE(x$mi_info$keep_fits)) "yes" else "no", "\n", sep = "")
  if (!is.null(x$mi_info$seed)) {
    cat("  Master seed: ", x$mi_info$seed, "\n", sep = "")
  }

  cat("\n")
  qgcompmulti_print_mixtures(
    list(
      labels = x$labels$mixture_labels,
      mix1 = x$mixtures$mix1,
      mix2 = x$mixtures$mix2
    )
  )

  cat("\n")
  qgcompmulti_print_msm_table(msm_table)

  invisible(x)
}
