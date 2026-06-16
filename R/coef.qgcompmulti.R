#' Core extractor methods for qgcompmulti objects
#'
#' Extract core fitted-model components from a `qgcompmulti` object using
#' standard R generic functions. These methods give the primary marginal
#' structural model (MSM) results and the metadata needed for model inspection.
#'
#' @param object A fitted `qgcompmulti` object for `coef()`, `vcov()`,
#'   `nobs()`, `df.residual()`, and `residuals()`.
#' @param x A fitted `qgcompmulti` object for `formula()`.
#' @param ... Not used.
#'
#' @return
#' `coef()` returns a named numeric vector of MSM coefficients.
#'
#' `vcov()` returns the full covariance matrix aligned with `coef(object)`.
#'
#' `formula()` returns the original formula supplied to
#' [qgcomp.glm.multi()].
#'
#' `nobs()` returns the number of observations actually used in the
#' model-fitting workflow.
#'
#' `df.residual()` returns a coefficient-level residual degrees-of-freedom
#' approximation based on the stored MSM parameter count and the number of
#' observations used. This is primarily intended to support `mice::pool()`
#' complete-data degrees-of-freedom handling for single-fit
#' `qgcompmulti` analyses.
#'
#' `residuals()` returns the residual vector from the stored outcome-model fit.
#' This preserves access to ordinary GLM residual behavior for downstream tools
#' that expect a conventional fitted-model residual interface.
#'
#' @seealso [summary.qgcompmulti()], [print.qgcompmulti()],
#'   [confint.qgcompmulti()], [qgcomp.glm.multi()]
#' @name qgcompmulti-extractors
#'
#' @export

coef.qgcompmulti <- function(object, ...) {
  validate_qgcompmulti(object)
  object$results$coefficients
}
