#' @rdname qgcompmulti-extractors
#' @export

vcov.qgcompmulti <- function(object, ...) {
  validate_qgcompmulti(object)
  object$results$vcov
}
