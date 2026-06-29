# ------------------------------------------------------------------------------
# Confidence-region helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_region_methods <- function() {
  "bootstrap_chisq"
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_region_method <- function(method) {
  allowed <- qgcompmulti_region_methods()

  if (!is.character(method) ||
      length(method) != 1L ||
      is.na(method) ||
      !method %in% allowed) {
    stop(
      sprintf(
        "`method` must be one of: %s.",
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(method)
}

#' @keywords internal
#' @noRd
qgcompmulti_resolve_region_parm <- function(parm, coef_names) {
  resolved <- qgcompmulti_resolve_parm(parm = parm, coef_names = coef_names)

  if (length(resolved) < 2L) {
    stop(
      "Confidence regions require at least two selected coefficients.",
      call. = FALSE
    )
  }

  resolved
}

#' @keywords internal
#' @noRd
new_qgcompmulti_region <- function(parm,
                                   center,
                                   covariance,
                                   level = 0.95,
                                   method = "bootstrap_chisq",
                                   threshold = NULL,
                                   plot_data = NULL) {
  object <- list(
    parm = parm,
    center = center,
    covariance = covariance,
    level = level,
    method = method,
    threshold = threshold,
    plot_data = plot_data
  )

  validate_qgcompmulti_region(object)

  structure(object, class = "qgcompmulti_region")
}

#' @keywords internal
#' @noRd
validate_qgcompmulti_region <- function(object) {
  if (!is.list(object)) {
    stop("`object` must be a list.", call. = FALSE)
  }

  required <- c("parm", "center", "covariance", "level", "method", "threshold", "plot_data")
  if (!identical(names(object), required)) {
    stop(
      sprintf(
        "`qgcompmulti_region` objects must contain fields in this order: %s.",
        paste(required, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.character(object$parm) || length(object$parm) < 2L || anyNA(object$parm)) {
    stop("`parm` must be a character vector with at least two coefficient names.", call. = FALSE)
  }

  if (!is.numeric(object$center) ||
      is.null(names(object$center)) ||
      !identical(names(object$center), object$parm)) {
    stop("`center` must be a named numeric vector aligned with `parm`.", call. = FALSE)
  }

  if (!is.matrix(object$covariance) ||
      is.null(rownames(object$covariance)) ||
      is.null(colnames(object$covariance)) ||
      !identical(rownames(object$covariance), object$parm) ||
      !identical(colnames(object$covariance), object$parm)) {
    stop("`covariance` must be a named square matrix aligned with `parm`.", call. = FALSE)
  }

  qgcompmulti_validate_conf_level(object$level)
  qgcompmulti_validate_region_method(object$method)

  if (!is.null(object$threshold) &&
      (!is.numeric(object$threshold) ||
       length(object$threshold) != 1L ||
       is.na(object$threshold) ||
       !is.finite(object$threshold) ||
       object$threshold <= 0)) {
    stop("`threshold` must be `NULL` or a single positive finite number.", call. = FALSE)
  }

  if (!is.null(object$plot_data) && !is.data.frame(object$plot_data)) {
    stop("`plot_data` must be `NULL` or a data frame.", call. = FALSE)
  }

  invisible(object)
}
