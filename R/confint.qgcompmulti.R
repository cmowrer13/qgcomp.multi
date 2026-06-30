#' Confidence intervals for qgcompmulti coefficients
#'
#' Returns confidence intervals for the marginal structural model (MSM)
#' coefficients reported by a fitted `qgcompmulti` object. Wald, percentile,
#' and basic bootstrap intervals are supported for single-fit objects. For
#' ratio estimands, intervals are computed on the MSM fitting scale and then
#' transformed to the user-facing ratio scale.
#'
#' @param object A fitted `qgcompmulti` object.
#' @param parm Optional specification of which coefficients to include. May be
#' `NULL` for all coefficients, an integer vector of coefficient positions, or
#' a character vector of coefficient names.
#' @param level Confidence level for the returned intervals. Must be strictly
#' between 0 and 1.
#' @param method Optional interval method. `NULL` uses the fitted object's
#'   stored default interval method. Supported values are `"wald"`,
#'   `"percentile"`, and `"basic"`. The basic method is the reverse-percentile
#'   bootstrap interval.
#' @param ... Not used.
#'
#' @return A numeric matrix with one row per selected coefficient and two
#' columns giving the lower and upper confidence limits on the display scale
#' recorded in `object$analysis$estimand_scale`.
#'
#' @details
#' Wald confidence intervals are based on the MSM coefficient estimates and the
#' stored bootstrap covariance matrix returned by [vcov.qgcompmulti()].
#' Percentile intervals use empirical quantiles of the stored bootstrap
#' coefficient draws. Basic intervals use the reverse-percentile construction,
#' reflecting the bootstrap quantiles around the full-sample estimate.
#'
#' All interval calculations are carried out on the MSM fitting scale. When the
#' active estimand is `"odds_ratio"` or `"rate_ratio"`, the returned limits are
#' exponentiated for display. The rows of the returned matrix align with
#' [coef.qgcompmulti()], but the numeric scale may differ for ratio estimands.
#' Pooled multiple-imputation objects intentionally remain Wald-only.
#'
#' @seealso [coef.qgcompmulti()], [vcov.qgcompmulti()],
#'   [summary.qgcompmulti()], [qgcomp.glm.multi()]
#' @export

confint.qgcompmulti <- function(object, parm = NULL, level = 0.95, method = NULL, ...) {
  validate_qgcompmulti(object)
  qgcompmulti_validate_conf_level(level)

  method <- qgcompmulti_resolve_interval_method(
    method = method,
    default_method = object$analysis$default_interval_method,
    context = "single_fit"
  )

  coefficients <- coef(object)
  coef_names <- qgcompmulti_resolve_parm(
    parm = parm,
    coef_names = names(coefficients)
  )

  coefficients <- coefficients[coef_names]
  vc <- vcov(object)[coef_names, coef_names, drop = FALSE]
  std_error <- setNames(sqrt(diag(vc)), coef_names)

  ci <- qgcompmulti_build_single_fit_confint(
    coefficients = coefficients,
    std_error = std_error,
    coef_draws = object$bootstrap$coef_draws,
    level = level,
    method = method
  )

  ci[] <- qgcompmulti_transform_msm_coefficients(
    as.numeric(ci),
    estimand_scale = object$analysis$estimand_scale,
    direction = "to_display"
  )
  ci
}
