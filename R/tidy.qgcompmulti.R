#' Tidy coefficient summaries for qgcompmulti fits
#'
#' Returns a broom-style coefficient table for the marginal structural model
#' (MSM) coefficients stored in a fitted `qgcompmulti` object. The output is
#' coefficient-centric and uses the internal machine-readable MSM term names.
#'
#' @param x A fitted `qgcompmulti` object.
#' @param conf.int Logical; if `TRUE`, add Wald confidence interval columns.
#' @param conf.level Confidence level for interval columns when `conf.int = TRUE`.
#' @param ... Not used.
#'
#' @return A data frame with one row per MSM coefficient and columns `term`,
#'   `estimate`, `std.error`, `statistic`, and `p.value`. When
#'   `conf.int = TRUE`, the returned data frame also includes `conf.low` and
#'   `conf.high`.
#'
#' @details
#' This method is intentionally coefficient-centric. It does not add fit-level
#' metadata columns or presentation-oriented labels because those belong in
#' [broom::glance()] or the existing print and summary methods.
#'
#' @seealso [broom::glance()], [coef.qgcompmulti()], [vcov.qgcompmulti()],
#'   [confint.qgcompmulti()], [summary.qgcompmulti()], [qgcomp.glm.multi()]
#' @exportS3Method broom::tidy
tidy.qgcompmulti <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  validate_qgcompmulti(x)
  if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int)) {
    stop("`conf.int` must be either `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (isTRUE(conf.int)) {
    qgcompmulti_validate_conf_level(conf.level)
  }
  coefficients <- coef(x)
  vc <- vcov(x)
  std_error <- sqrt(diag(vc))
  std_error <- setNames(as.numeric(std_error), names(coefficients))
  statistic <- coefficients / std_error
  p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  out <- data.frame(
    term = names(coefficients),
    estimate = unname(coefficients),
    std.error = unname(std_error),
    statistic = unname(statistic),
    p.value = unname(p_value),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (isTRUE(conf.int)) {
    ci <- stats::confint(x, level = conf.level)
    out$conf.low <- unname(ci[, 1])
    out$conf.high <- unname(ci[, 2])
  }
  out
}
