#' Print a qgcompmulti fit
#'
#' Prints a moderately detailed overview of a fitted `qgcompmulti` object for
#' ordinary interactive use. The display highlights the marginal structural
#' model (MSM) results while keeping the full outcome-model context for
#' [summary.qgcompmulti()].
#'
#' @param x A fitted `qgcompmulti` object.
#' @param ... Not used.
#'
#' @return The input object, invisibly.
#'
#' @details
#' `print.qgcompmulti()` is designed for ordinary interactive use. It gives a
#' moderately detailed overview of the fitted model, including the call,
#' analysis mode, mixture definitions, and labeled MSM coefficient table.
#' For a fuller report, use [summary.qgcompmulti()].
#'
#' @seealso [summary.qgcompmulti()], [coef.qgcompmulti()],
#'   [confint.qgcompmulti()], [qgcomp.glm.multi()]
#' @export

print.qgcompmulti <- function(x, ...) {
  validate_qgcompmulti(x)

  msm_table <- qgcompmulti_labeled_coef_table(x$results, x$labels, x$analysis)

  cat("qgcompmulti fit\n\n")
  qgcompmulti_print_call(x$call)

  cat("\nModel:\n")
  cat("  Outcome: ", x$data_info$outcome, "\n", sep = "")
  cat("  Family: ", qgcompmulti_family_label(x$analysis), "\n", sep = "")
  qgcompmulti_print_scale_info(x$analysis)
  cat("  Observations used: ", x$data_info$n_used, "\n", sep = "")
  cat(
    "  Exposure mode: ",
    qgcompmulti_exposure_mode_label(x$mixtures),
    "\n",
    sep = ""
  )
  cat("  ", qgcompmulti_interaction_label(x$analysis), "\n", sep = "")
  if (!is.null(x$analysis$seed)) {
    cat("  Random seed: ", x$analysis$seed, "\n", sep = "")
  }
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
