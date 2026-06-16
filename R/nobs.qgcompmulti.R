#' @importFrom stats nobs
#' @rdname qgcompmulti-extractors
#' @export

nobs.qgcompmulti <- function(object, ...) {
  validate_qgcompmulti(object)
  object$data_info$n_used
}

#' @rdname qgcompmulti-extractors
#' @export
df.residual.qgcompmulti <- function(object, ...) {
  validate_qgcompmulti(object)
  as.numeric(max(nobs(object) - length(coef(object)), 1L))
}
