#' Diagnostic summaries for a qgcompmulti fit
#'
#' Returns structured diagnostics for intervention support, bootstrap behavior,
#' and MSM adequacy from a fitted [qgcomp.glm.multi()] object.
#'
#' These diagnostics are deliberately kept outside the main fitted-model
#' summary so that users can inspect model fit, support, and computational
#' behavior explicitly rather than having a very dense `summary()` method.
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
#' @details
#' The three main diagnostic families are:
#'
#' \itemize{
#'   \item \strong{Support diagnostics}, which summarize the fit-time
#'   intervention grid. For `q = NULL`, these diagnostics are especially
#'   important because the intervention grid is built from pooled percentile
#'   values within each mixture, and under each intervention every component in
#'   a mixture is set to the same pooled mixture-specific value.
#'   \item \strong{Bootstrap diagnostics}, which summarize how many bootstrap
#'   replications were requested, retained, or failed, together with any
#'   lightweight failure metadata stored in the fitted object.
#'   \item \strong{MSM adequacy diagnostics}, which compare the exact fit-time
#'   counterfactual surface to the fitted MSM surface. For transformed-scale
#'   fits, the primary adequacy comparison is made on the MSM fitting scale.
#' }
#'
#' Support diagnostics should not be read as a full positivity proof. They are
#' designed to help users see which intervention values define the stored
#' surface and to highlight when original-scale pooled interventions may deserve
#' extra scrutiny.
#'
#' MSM adequacy is one of the most method-specific diagnostics in the package.
#' It asks whether the fitted MSM is a reasonable low-dimensional summary of the
#' exact fit-time surface implied by the fitted outcome model. It is not a test
#' of whether the outcome model is true, and it is not a general check of causal
#' identification.
#'
#' Adequacy is evaluated on the stored fit-time grid. A good adequacy result
#' therefore means that the MSM tracks the exact fit-time surface well on that
#' grid. It does not imply that the surface is globally linear between or beyond
#' those intervention points.
#'
#' @examples
#' \dontrun{
#' dat <- sim_mixture_data(
#'   n = 400,
#'   pA = 3,
#'   pB = 3,
#'   rho_within_A = 0.3,
#'   rho_within_B = 0.3,
#'   rho_between = 0.2,
#'   psi1 = 0.5,
#'   psi2 = 0.3,
#'   psi12 = 0.2,
#'   seed = 123
#' )
#'
#' fit <- qgcomp.glm.multi(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   q = 4,
#'   B = 100,
#'   seed = 13
#' )
#'
#' # Full diagnostic bundle
#' diagnostics(fit)
#'
#' # Focused diagnostics
#' support(fit)
#' diagnostics(fit, type = "bootstrap")
#' adequacy(fit)
#' }
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
#'
#' @details
#' The support diagnostic is particularly informative for `q = NULL`, where the
#' stored intervention grid is defined by pooled percentile values rather than
#' quantile indices. In that setting, every component in a mixture is set to the
#' same pooled mixture-specific intervention value, so users should inspect the
#' resulting grid for scientific plausibility as well as numerical range.
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
#'
#' @details
#' The adequacy diagnostic evaluates how closely the fitted MSM reproduces the
#' exact fit-time counterfactual surface stored in the fitted object. This is a
#' diagnostic of MSM approximation on the stored fit-time grid, not a test of
#' whether the outcome model is correctly specified relative to the true
#' data-generating process.
#'
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
  cat("Comparison scale: ", x$comparison_scale, "\n", sep = "")
  cat("MSM fitting scale: ", x$msm_fitting_scale, "\n", sep = "")
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
