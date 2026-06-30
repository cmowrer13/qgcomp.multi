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
                                   df = NULL,
                                   coefficient_labels = NULL,
                                   estimand_scale = "mean_difference",
                                   msm_fitting_scale = "identity",
                                   plot_data = NULL) {
  if (is.null(df)) {
    df <- length(parm)
  }
  if (is.null(threshold)) {
    threshold <- stats::qchisq(level, df = df)
  }
  if (is.null(coefficient_labels)) {
    coefficient_labels <- stats::setNames(parm, parm)
  }

  object <- list(
    parm = parm,
    center = center,
    covariance = covariance,
    level = level,
    method = method,
    threshold = threshold,
    df = as.integer(df),
    coefficient_labels = coefficient_labels,
    estimand_scale = estimand_scale,
    msm_fitting_scale = msm_fitting_scale,
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

  required <- c(
    "parm",
    "center",
    "covariance",
    "level",
    "method",
    "threshold",
    "df",
    "coefficient_labels",
    "estimand_scale",
    "msm_fitting_scale",
    "plot_data"
  )
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

  qgcompmulti_validate_estimand_scale_name(object$estimand_scale)
  qgcompmulti_msm_fitting_scale_label(object$msm_fitting_scale)
  if (!is.numeric(object$threshold) ||
      length(object$threshold) != 1L ||
      is.na(object$threshold) ||
      !is.finite(object$threshold) ||
      object$threshold <= 0) {
    stop("`threshold` must be a single positive finite number.", call. = FALSE)
  }
  if (!is.integer(object$df) ||
      length(object$df) != 1L ||
      is.na(object$df) ||
      object$df < 2L ||
      object$df != length(object$parm)) {
    stop("`df` must be an integer equal to the number of selected coefficients.", call. = FALSE)
  }
  if (!is.character(object$coefficient_labels) ||
      is.null(names(object$coefficient_labels)) ||
      !identical(names(object$coefficient_labels), object$parm) ||
      anyNA(object$coefficient_labels)) {
    stop(
      "`coefficient_labels` must be a named character vector aligned with `parm`.",
      call. = FALSE
    )
  }

  if (!is.null(object$plot_data) && !is.data.frame(object$plot_data)) {
    stop("`plot_data` must be `NULL` or a data frame.", call. = FALSE)
  }

  invisible(object)
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_region_npoints <- function(npoints) {
  if (!is.numeric(npoints) ||
      length(npoints) != 1L ||
      is.na(npoints) ||
      !is.finite(npoints) ||
      npoints != as.integer(npoints) ||
      npoints < 4L) {
    stop("`npoints` must be a single integer greater than or equal to 4.", call. = FALSE)
  }
  as.integer(npoints)
}
#' @keywords internal
#' @noRd
qgcompmulti_validate_region_covariance <- function(covariance,
                                                   parm,
                                                   tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(covariance) ||
      !is.numeric(covariance) ||
      !identical(rownames(covariance), parm) ||
      !identical(colnames(covariance), parm)) {
    stop(
      "`covariance` must be a named numeric square matrix aligned with the selected coefficients.",
      call. = FALSE
    )
  }
  if (any(!is.finite(covariance))) {
    stop("The restricted bootstrap covariance matrix contains non-finite values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(covariance, t(covariance), tolerance = tol))) {
    stop("The restricted bootstrap covariance matrix must be symmetric.", call. = FALSE)
  }
  eigen_values <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigen_values <= tol)) {
    stop(
      paste(
        "The restricted bootstrap covariance matrix is not positive definite.",
        "A confidence-region ellipsoid cannot be constructed for this parameter subset."
      ),
      call. = FALSE
    )
  }
  invisible(covariance)
}
#' @keywords internal
#' @noRd
qgcompmulti_region_plot_data <- function(center,
                                         covariance,
                                         threshold,
                                         npoints = 200L) {
  if (length(center) != 2L) {
    return(NULL)
  }
  npoints <- qgcompmulti_validate_region_npoints(npoints)
  qgcompmulti_validate_region_covariance(
    covariance = covariance,
    parm = names(center)
  )
  eig <- eigen(covariance, symmetric = TRUE)
  angles <- seq(0, 2 * pi, length.out = npoints)
  unit_circle <- rbind(cos(angles), sin(angles))
  transform <- eig$vectors %*% diag(sqrt(eig$values * threshold), nrow = 2L)
  boundary <- t(transform %*% unit_circle)
  boundary <- sweep(boundary, 2L, center, FUN = "+")
  data.frame(
    x = boundary[, 1L],
    y = boundary[, 2L],
    row.names = NULL
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_region_method_label <- function(method) {
  qgcompmulti_validate_region_method(method)
  switch(
    method,
    bootstrap_chisq = "Bootstrap covariance chi-squared ellipsoid"
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_region_scale_note <- function(estimand_scale, msm_fitting_scale) {
  qgcompmulti_validate_estimand_scale_name(estimand_scale)
  qgcompmulti_msm_fitting_scale_label(msm_fitting_scale)
  if (qgcompmulti_is_ratio_estimand(estimand_scale)) {
    return(
      paste(
        "Region geometry is computed on the MSM fitting coefficient scale.",
        "For this ratio estimand, that is the log coefficient scale; the",
        "ellipsoid is not geometrically exponentiated."
      )
    )
  }
  "Region geometry is computed on the MSM fitting coefficient scale."
}
#' @keywords internal
#' @noRd
qgcompmulti_region_axis_label <- function(label, estimand_scale, msm_fitting_scale) {
  qgcompmulti_validate_estimand_scale_name(estimand_scale)
  fitting_scale <- qgcompmulti_msm_fitting_scale_label(msm_fitting_scale)
  estimand <- tolower(qgcompmulti_estimand_label(estimand_scale))
  if (qgcompmulti_is_ratio_estimand(estimand_scale)) {
    return(sprintf("%s (log %s coefficient)", label, estimand))
  }
  sprintf("%s (%s coefficient)", label, fitting_scale)
}
