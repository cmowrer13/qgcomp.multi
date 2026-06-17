#' Core extractor methods for pooled qgcompmulti multiple-imputation objects
#'
#' Extract the pooled marginal structural model (MSM) coefficient results from a
#' fitted `qgcompmulti_mi` object using standard R generics.
#'
#' @param object A fitted `qgcompmulti_mi` object.
#' @param parm Optional specification of which coefficients to include. May be
#'   `NULL` for all coefficients, an integer vector of coefficient positions, or
#'   a character vector of coefficient names.
#' @param level Confidence level for the returned intervals. Must be strictly
#'   between 0 and 1.
#' @param ... Not used.
#'
#' @return
#' `coef()` returns a named numeric vector of pooled MSM coefficients.
#'
#' `vcov()` returns the pooled covariance matrix aligned with `coef(object)`.
#'
#' `confint()` returns pooled Wald-style confidence intervals that use
#' term-specific t critical values when finite Rubin degrees of freedom are
#' available and normal critical values otherwise.
#'
#' @details
#' These extractors operate on the pooled multiple-imputation result, not on
#' the individual completed-data fits. They therefore return Rubin-pooled MSM
#' coefficient summaries for the inferential target defined by
#' [qgcomp.glm.multi.mi()], rather than per-imputation coefficient tables,
#' pooled prediction objects, or pooled diagnostics. If you need to inspect the
#' stored completed-data fits directly, fit with `keep_fits = TRUE` and extract
#' the individual `"qgcompmulti"` objects from `object$fits$imputation_fits`.
#'
#' @seealso [summary.qgcompmulti_mi()], [print.qgcompmulti_mi()],
#'   [broom::tidy()], [broom::glance()], [qgcomp.glm.multi.mi()]
#' @name qgcompmulti-mi-extractors
#'
#' @export
coef.qgcompmulti_mi <- function(object, ...) {
  validate_qgcompmulti_mi(object)
  object$results$coefficients
}

#' @rdname qgcompmulti-mi-extractors
#' @export
vcov.qgcompmulti_mi <- function(object, ...) {
  validate_qgcompmulti_mi(object)
  object$results$vcov
}

#' @rdname qgcompmulti-mi-extractors
#' @export
confint.qgcompmulti_mi <- function(object, parm = NULL, level = 0.95, ...) {
  validate_qgcompmulti_mi(object)
  qgcompmulti_validate_conf_level(level)

  coefficients <- coef(object)
  coef_names <- qgcompmulti_resolve_parm(
    parm = parm,
    coef_names = names(coefficients)
  )

  coefficients <- coefficients[coef_names]
  std_error <- object$results$std_error[coef_names]
  df <- object$results$df[coef_names]

  build_qgcompmulti_mi_confint(
    coefficients = coefficients,
    std_error = std_error,
    df = df,
    level = level
  )
}
