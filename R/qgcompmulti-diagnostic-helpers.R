# ------------------------------------------------------------------------------
# Diagnostic helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
build_qgcompmulti_bootstrap <- function(coef_draws,
                                        B_requested,
                                        failure_log = NULL) {
  if (is.null(coef_draws)) {
    B_success <- 0L
  } else {
    B_success <- nrow(as.matrix(coef_draws))
  }
  B_failed <- max(0L, as.integer(B_requested) - B_success)
  list(
    coef_draws = coef_draws,
    B_requested = as.integer(B_requested),
    B_success = B_success,
    B_failed = B_failed,
    failure_log = failure_log
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_validate_bootstrap <- function(bootstrap, labels) {
  coef_names <- labels$coefficient_names
  if (!is.null(bootstrap$coef_draws)) {
    if (!is.matrix(bootstrap$coef_draws) && !is.data.frame(bootstrap$coef_draws)) {
      stop("`bootstrap$coef_draws` must be `NULL`, a matrix, or a data frame.", call. = FALSE)
    }
    coef_draws <- as.matrix(bootstrap$coef_draws)
    if (is.null(colnames(coef_draws))) {
      stop("`bootstrap$coef_draws` must have column names.", call. = FALSE)
    }
    if (!is.null(coef_names) && !identical(colnames(coef_draws), coef_names)) {
      stop(
        "Column names of `bootstrap$coef_draws` must match `labels$coefficient_names`.",
        call. = FALSE
      )
    }
  }
  scalar_fields <- c("B_requested", "B_success", "B_failed")
  for (field in scalar_fields) {
    value <- bootstrap[[field]]
    if (!is.null(value) && length(value) != 1L) {
      stop(sprintf("`bootstrap$%s` must be length 1.", field), call. = FALSE)
    }
  }
  if (!is.null(bootstrap$coef_draws) &&
      !is.null(bootstrap$B_success) &&
      nrow(as.matrix(bootstrap$coef_draws)) != bootstrap$B_success) {
    stop("`bootstrap$B_success` must match the number of retained bootstrap draws.", call. = FALSE)
  }
  if (!is.null(bootstrap$B_requested) &&
      !is.null(bootstrap$B_success) &&
      !is.null(bootstrap$B_failed) &&
      bootstrap$B_failed != (bootstrap$B_requested - bootstrap$B_success)) {
    stop("`bootstrap$B_failed` must equal `B_requested - B_success`.", call. = FALSE)
  }
  if (!is.null(bootstrap$failure_log) && !is.data.frame(bootstrap$failure_log)) {
    stop("`bootstrap$failure_log` must be `NULL` or a data frame.", call. = FALSE)
  }
  invisible(bootstrap)
}
