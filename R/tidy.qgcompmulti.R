#' Tidy coefficient summaries for qgcompmulti fits
#'
#' Returns a broom-style coefficient table for the marginal structural model
#' (MSM) coefficients stored in a fitted `qgcompmulti` object. The core
#' `estimate` and `std.error` columns stay on the MSM fitting scale so
#' downstream tools such as `mice::pool()` can use them coherently. Display
#' columns expose the active estimand scale for ordinary reporting.
#'
#' @param x A fitted `qgcompmulti` object.
#' @param conf.int Logical; if `TRUE`, add confidence interval columns.
#' @param conf.level Confidence level for interval columns when `conf.int = TRUE`.
#' @param method Optional interval method for confidence intervals. `NULL` uses
#'   the fitted object's stored default interval method. Supported values are
#'   `"wald"`, `"percentile"`, and `"basic"`.
#' @param ... Not used.
#'
#' @return A data frame with one row per MSM coefficient. Columns `estimate`,
#'   `std.error`, `statistic`, and `p.value` are on the fitting scale. Columns
#'   `display.estimate`, `display.conf.low`, and `display.conf.high` are on the
#'   active estimand scale when present.
#'
#' @details
#' Keeping `estimate` and `std.error` on the fitting scale preserves machine
#' pooling behavior for multiple-imputation workflows that call
#' [mice::pool()]. For odds-ratio and rate-ratio estimands, use
#' `display.estimate` and the display confidence interval columns for
#' user-facing ratio summaries.
#'
#' @seealso [broom::glance()], [coef.qgcompmulti()], [vcov.qgcompmulti()],
#'   [confint.qgcompmulti()], [summary.qgcompmulti()], [qgcomp.glm.multi()]
#' @exportS3Method broom::tidy
tidy.qgcompmulti <- function(x, conf.int = FALSE, conf.level = 0.95, method = NULL, ...) {
  validate_qgcompmulti(x)
  if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int)) {
    stop("`conf.int` must be either `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (isTRUE(conf.int)) {
    qgcompmulti_validate_conf_level(conf.level)
    method <- qgcompmulti_resolve_interval_method(
      method = method,
      default_method = x$analysis$default_interval_method,
      context = "single_fit"
    )
  }
  coefficients <- coef(x)
  vc <- vcov(x)
  std_error <- sqrt(diag(vc))
  std_error <- setNames(as.numeric(std_error), names(coefficients))
  statistic <- coefficients / std_error
  p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
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
    p.value = unname(p_value),
    estimand_scale = x$analysis$estimand_scale,
    msm_fitting_scale = x$analysis$msm_fitting_scale,
    estimate_scale = "fitting",
    display_scale = x$analysis$estimand_scale,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (isTRUE(conf.int)) {
    fit_ci <- qgcompmulti_build_single_fit_confint(
      coefficients = coefficients,
      std_error = std_error,
      coef_draws = x$bootstrap$coef_draws,
      level = conf.level,
      method = method
    )
    display_ci <- stats::confint(x, level = conf.level, method = method)
    out$conf.low <- unname(fit_ci[, 1])
    out$conf.high <- unname(fit_ci[, 2])
    out$display.conf.low <- unname(display_ci[, 1])
    out$display.conf.high <- unname(display_ci[, 2])
  }
  out
}
