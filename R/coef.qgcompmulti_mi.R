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
#' @param method Optional interval method. Pooled multiple-imputation
#'   coefficient reporting supports only `"wald"` in Version 0.5.0. `NULL` uses
#'   the fitted object's stored default, which is `"wald"` for pooled MI fits.
#' @param ... Not used.
#'
#' @return
#' `coef()` returns a named numeric vector of pooled MSM coefficients on the
#' fitting scale.
#'
#' `vcov()` returns the pooled covariance matrix aligned with `coef(object)`.
#'
#' `confint()` returns pooled Wald-style confidence intervals. The Wald
#' calculation uses Rubin-pooled coefficients, standard errors, and
#' term-specific degrees of freedom on the fitting scale. For odds-ratio and
#' rate-ratio estimands, returned interval limits are exponentiated for display.
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
confint.qgcompmulti_mi <- function(object, parm = NULL, level = 0.95, method = NULL, ...) {
  validate_qgcompmulti_mi(object)
  qgcompmulti_validate_conf_level(level)

  method <- qgcompmulti_resolve_interval_method(
    method = method,
    default_method = object$analysis$default_interval_method,
    context = "mi"
  )
  qgcompmulti_require_wald_interval_method(
    method = method,
    object_label = "pooled qgcompmulti multiple-imputation coefficient"
  )

  coefficients <- coef(object)
  coef_names <- qgcompmulti_resolve_parm(
    parm = parm,
    coef_names = names(coefficients)
  )

  coefficients <- coefficients[coef_names]
  std_error <- object$results$std_error[coef_names]
  df <- object$results$df[coef_names]

  ci <- build_qgcompmulti_mi_confint(
    coefficients = coefficients,
    std_error = std_error,
    df = df,
    level = level
  )

  ci[] <- qgcompmulti_transform_msm_coefficients(
    as.numeric(ci),
    estimand_scale = object$analysis$estimand_scale,
    direction = "to_display"
  )
  ci
}
