#' Wald confidence intervals for qgcompmulti coefficients
#'
#' Returns Wald confidence intervals for the marginal structural model (MSM)
#' coefficients reported by a fitted `qgcompmulti` object.
#'
#' @param object A fitted `qgcompmulti` object.
#' @param parm Optional specification of which coefficients to include. May be
#' `NULL` for all coefficients, an integer vector of coefficient positions, or
#' a character vector of coefficient names.
#' @param level Confidence level for the returned intervals. Must be strictly
#' between 0 and 1.
#' @param ... Not used
#'
#' @return A numeric matrix with one row per selected coefficient and two
#' columns giving the lower and upper Wald confidence limits.
#'
#' @details
#' Confidence intervals are Wald intervals based on the MSM coefficient
#' estimates and the stored covariance matrix returned by [vcov.qgcompmulti()].
#' The rows of the returned matrix align with [coef.qgcompmulti()].
#'
#' @seealso [coef.qgcompmulti()], [vcov.qgcompmulti()],
#'   [summary.qgcompmulti()], [qgcomp.glm.multi()]
#' @export

confint.qgcompmulti <- function(object, parm = NULL, level = 0.95, ...) {
  validate_qgcompmulti(object)
  qgcompmulti_validate_conf_level(level)

  coefficients <- coef(object)
  coef_names <- qgcompmulti_resolve_parm(
    parm = parm,
    coef_names = names(coefficients)
  )

  coefficients <- coefficients[coef_names]
  vc <- vcov(object)[coef_names, coef_names, drop = FALSE]
  std_error <- setNames(sqrt(diag(vc)), coef_names)

  build_qgcompmulti_confint(
    coefficients = coefficients,
    std_error = std_error,
    level = level
  )
}
