#' @importFrom stats nobs
#' @rdname qgcompmulti-extractors
#' @export

nobs.qgcompmulti <- function(object, ...) {
  validate_qgcompmulti(object)
  object$data_info$n_used
}
