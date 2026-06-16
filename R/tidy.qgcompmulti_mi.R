#' Tidy coefficient summaries for pooled qgcompmulti multiple-imputation fits
#'
#' Returns a broom-style coefficient table for the pooled marginal structural
#' model (MSM) coefficients stored in a fitted `qgcompmulti_mi` object.
#'
#' @param x A fitted `qgcompmulti_mi` object.
#' @param conf.int Logical; if `TRUE`, add confidence interval columns.
#' @param conf.level Confidence level for interval columns when
#'   `conf.int = TRUE`.
#' @param ... Not used.
#'
#' @return A data frame with one row per pooled MSM coefficient and columns
#'   `term`, `estimate`, `std.error`, `statistic`, `df`, and `p.value`. When
#'   `conf.int = TRUE`, the returned data frame also includes `conf.low` and
#'   `conf.high`.
#'
#' @seealso [broom::glance()], [coef.qgcompmulti_mi()], [vcov.qgcompmulti_mi()],
#'   [confint.qgcompmulti_mi()], [summary.qgcompmulti_mi()],
#'   [qgcomp.glm.multi.mi()]
#' @exportS3Method broom::tidy
tidy.qgcompmulti_mi <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  validate_qgcompmulti_mi(x)

  if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int)) {
    stop("`conf.int` must be either `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (isTRUE(conf.int)) {
    qgcompmulti_validate_conf_level(conf.level)
  }

  coefficients <- coef(x)
  std_error <- x$results$std_error[names(coefficients)]
  df <- x$results$df[names(coefficients)]
  statistic <- coefficients / std_error
  p_value <- qgcompmulti_mi_p_values(statistic = statistic, df = df)

  out <- data.frame(
    term = names(coefficients),
    estimate = unname(coefficients),
    std.error = unname(std_error),
    statistic = unname(statistic),
    df = unname(df),
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
