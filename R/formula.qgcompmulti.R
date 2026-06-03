#' @rdname qgcompmulti-extractors
#' @export

formula.qgcompmulti <- function(x, ...) {
  validate_qgcompmulti(x)
  x$formula
}
