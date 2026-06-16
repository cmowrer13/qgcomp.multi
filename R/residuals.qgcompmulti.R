#' @rdname qgcompmulti-extractors
#' @export
residuals.qgcompmulti <- function(object, ...) {
  validate_qgcompmulti(object)

  if (is.null(object$fits$outcome_fit)) {
    return(NULL)
  }

  stats::residuals(object$fits$outcome_fit, ...)
}
