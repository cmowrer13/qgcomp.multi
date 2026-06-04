#' Diagnostic summaries for a qgcompmulti fit
#'
#' Returns structured diagnostics for intervention support, bootstrap behavior,
#' and MSM adequacy from a fitted [qgcomp.glm.multi()] object.
#'
#' @param object A fitted `"qgcompmulti"` object.
#' @param type Character string indicating which diagnostic to return.
#' Supported values are `"all"`, `"support"`, `"bootstrap"`, and
#' `"adequacy"`.
#' @param ... Unused.
#'
#' @return A structured diagnostic object, or a named list of all diagnostics
#' when `type = "all"`.
#'
#' @export
diagnostics <- function(object, type = c("all", "support", "bootstrap", "adequacy"), ...) {
  UseMethod("diagnostics")
}
#' @export
diagnostics.qgcompmulti <- function(object,
                                    type = c("all", "support", "bootstrap", "adequacy"),
                                    ...) {
  validate_qgcompmulti(object)
  type <- match.arg(type)
  support_diag <- build_qgcompmulti_support_diagnostic(object)
  bootstrap_diag <- build_qgcompmulti_bootstrap_diagnostic(object)
  adequacy_diag <- build_qgcompmulti_adequacy_diagnostic(object)
  switch(
    type,
    all = structure(
      list(
        support = support_diag,
        bootstrap = bootstrap_diag,
        adequacy = adequacy_diag
      ),
      class = "qgcompmulti_diagnostics"
    ),
    support = support_diag,
    bootstrap = bootstrap_diag,
    adequacy = adequacy_diag
  )
}
#' Extract intervention support diagnostics from a qgcompmulti fit
#'
#' @param object A fitted `"qgcompmulti"` object.
#' @param ... Unused.
#'
#' @return An object of class `"qgcompmulti_support_diagnostic"`.
#' @export
support <- function(object, ...) {
  UseMethod("support")
}
#' @export
support.qgcompmulti <- function(object, ...) {
  validate_qgcompmulti(object)
  build_qgcompmulti_support_diagnostic(object)
}
#' Extract MSM adequacy diagnostics from a qgcompmulti fit
#'
#' @param object A fitted `"qgcompmulti"` object.
#' @param ... Unused.
#'
#' @return An object of class `"qgcompmulti_adequacy_diagnostic"`.
#' @export
adequacy <- function(object, ...) {
  UseMethod("adequacy")
}
#' @export
adequacy.qgcompmulti <- function(object, ...) {
  validate_qgcompmulti(object)
  build_qgcompmulti_adequacy_diagnostic(object)
}
#' @export
print.qgcompmulti_support_diagnostic <- function(x, ...) {
  summary <- x$support_summary
  flags <- x$flags
  cat("qgcompmulti intervention support diagnostic\n\n")
  cat("Mode: ", x$mode, "\n", sep = "")
  if (!is.null(x$centering)) {
    cat("Centering: ", x$centering, "\n", sep = "")
  }
  cat("Grid points: ", summary$n_grid_points, "\n", sep = "")
  cat(
    "Intervention psi1 range: [",
    qgcompmulti_format_axis_values(summary$intervention_psi1_range[1]),
    ", ",
    qgcompmulti_format_axis_values(summary$intervention_psi1_range[2]),
    "]\n",
    sep = ""
  )
  cat(
    "Intervention psi2 range: [",
    qgcompmulti_format_axis_values(summary$intervention_psi2_range[1]),
    ", ",
    qgcompmulti_format_axis_values(summary$intervention_psi2_range[2]),
    "]\n",
    sep = ""
  )
  if (isTRUE(flags$pooled_percentile_grid)) {
    cat("Support note:\n")
    for (note in x$notes) {
      cat("  - ", note, "\n", sep = "")
    }
  }
  invisible(x)
}
#' @export
print.qgcompmulti_bootstrap_diagnostic <- function(x, ...) {
  cat("qgcompmulti bootstrap diagnostic\n\n")
  cat("Requested replications: ", x$B_requested, "\n", sep = "")
  cat("Successful replications: ", x$B_success, "\n", sep = "")
  cat("Failed replications: ", x$B_failed, "\n", sep = "")
  cat("Success rate: ", qgcompmulti_format_metric(100 * x$success_rate), "%\n", sep = "")
  if (nrow(x$failure_summary) > 0L) {
    cat("\nFailure summary:\n")
    print(x$failure_summary, row.names = FALSE)
  }
  invisible(x)
}
#' @export
print.qgcompmulti_adequacy_diagnostic <- function(x, ...) {
  metrics <- x$summary_metrics
  cat("qgcompmulti MSM adequacy diagnostic\n\n")
  cat("Grid points: ", metrics$n_grid_points, "\n", sep = "")
  cat("Mean absolute error: ", qgcompmulti_format_metric(metrics$mae), "\n", sep = "")
  cat("RMSE: ", qgcompmulti_format_metric(metrics$rmse), "\n", sep = "")
  cat("Maximum absolute error: ", qgcompmulti_format_metric(metrics$max_abs_error), "\n", sep = "")
  cat("Mean signed error: ", qgcompmulti_format_metric(metrics$mean_signed_error), "\n", sep = "")
  cat("Correlation: ", qgcompmulti_format_metric(metrics$correlation), "\n", sep = "")
  cat("\n", x$notes, "\n", sep = "")
  invisible(x)
}
#' @export
print.qgcompmulti_diagnostics <- function(x, ...) {
  cat("qgcompmulti diagnostics\n\n")
  print(x$support)
  cat("\n")
  print(x$bootstrap)
  cat("\n")
  print(x$adequacy)
  invisible(x)
}
