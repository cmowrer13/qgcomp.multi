#' Tidy coefficient summaries for pooled qgcompmulti multiple-imputation fits
#'
#' Returns a broom-style coefficient table for the pooled marginal structural
#' model (MSM) coefficients stored in a fitted `qgcompmulti_mi` object. The core
#' `estimate` and `std.error` columns stay on the Rubin-pooled fitting scale,
#' while display columns expose the active estimand scale for reporting.
#'
#' @param x A fitted `qgcompmulti_mi` object.
#' @param conf.int Logical; if `TRUE`, add confidence interval columns.
#' @param conf.level Confidence level for interval columns when
#'   `conf.int = TRUE`.
#' @param method Optional interval method. Pooled multiple-imputation
#'   coefficient reporting supports only `"wald"` in Version 0.5.0.
#' @param ... Not used.
#'
#' @return A data frame with one row per pooled MSM coefficient. Columns
#'   `estimate`, `std.error`, `statistic`, `df`, and `p.value` are on the
#'   fitting scale. Columns `display.estimate`, `display.conf.low`, and
#'   `display.conf.high` are on the active estimand scale when present.
#'
#' @details
#' Pooled multiple-imputation inference is carried out on the fitting scale.
#' For odds-ratio and rate-ratio estimands, the display columns transform only
#' the final pooled estimates and interval limits. Rubin pooling, standard
#' errors, test statistics, and degrees of freedom remain on the fitting scale.
#'
#' @seealso [broom::glance()], [coef.qgcompmulti_mi()], [vcov.qgcompmulti_mi()],
#'   [confint.qgcompmulti_mi()], [summary.qgcompmulti_mi()],
#'   [qgcomp.glm.multi.mi()]
#' @exportS3Method broom::tidy
tidy.qgcompmulti_mi <- function(x, conf.int = FALSE, conf.level = 0.95, method = NULL, ...) {
  validate_qgcompmulti_mi(x)

  if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int)) {
    stop("`conf.int` must be either `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (isTRUE(conf.int)) {
    qgcompmulti_validate_conf_level(conf.level)
    method <- qgcompmulti_resolve_interval_method(
      method = method,
      default_method = x$analysis$default_interval_method,
      context = "mi"
    )
    qgcompmulti_require_wald_interval_method(
      method = method,
      object_label = "pooled qgcompmulti multiple-imputation tidy coefficient"
    )
  }

  coefficients <- coef(x)
  std_error <- x$results$std_error[names(coefficients)]
  df <- x$results$df[names(coefficients)]
  statistic <- coefficients / std_error
  p_value <- qgcompmulti_mi_p_values(statistic = statistic, df = df)
  display_estimate <- qgcompmulti_transform_msm_coefficients(
    coefficients,
    estimand_scale = x$analysis$estimand_scale,
    direction = "to_display"
  )

  out <- data.frame(
    term = names(coefficients),
    estimate = unname(coefficients),
    display.estimate = unname(display_estimate),
    std.error = unname(std_error),
    statistic = unname(statistic),
    df = unname(df),
    p.value = unname(p_value),
    estimand_scale = x$analysis$estimand_scale,
    msm_fitting_scale = x$analysis$msm_fitting_scale,
    estimate_scale = "fitting",
    display_scale = x$analysis$estimand_scale,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (isTRUE(conf.int)) {
    fit_ci <- build_qgcompmulti_mi_confint(
      coefficients = coefficients,
      std_error = std_error,
      df = df,
      level = conf.level
    )
    display_ci <- stats::confint(x, level = conf.level, method = method)
    out$conf.low <- unname(fit_ci[, 1])
    out$conf.high <- unname(fit_ci[, 2])
    out$display.conf.low <- unname(display_ci[, 1])
    out$display.conf.high <- unname(display_ci[, 2])
  }

  out
}
