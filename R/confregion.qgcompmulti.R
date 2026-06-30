#' Confidence regions for qgcompmulti coefficients
#'
#' Builds a bootstrap-covariance chi-squared confidence region for selected
#' marginal structural model (MSM) coefficients from a fitted
#' [qgcomp.glm.multi()] object.
#'
#' @param object A fitted object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A structured `qgcompmulti_region` object containing the selected
#' coefficient names, fitting-scale center, restricted bootstrap covariance
#' matrix, confidence level, chi-squared cutoff, degrees of freedom, scale
#' metadata, and two-dimensional plot data when applicable.
#' @export
confregion <- function(object, ...) {
  UseMethod("confregion")
}
#' @rdname confregion
#' @param parm Character or integer vector identifying the MSM coefficients to
#' include in the region. The default is `c("psi1", "psi2")`. Use
#' `NULL` to include all MSM coefficients.
#' @param level Confidence level for the chi-squared cutoff. Must be a single
#' number strictly between 0 and 1.
#' @param method Region construction method. Version `0.5.0` supports only
#' `"bootstrap_chisq"`, the chi-squared ellipsoid based on the bootstrap
#' covariance matrix of the selected MSM coefficients.
#' @param npoints Number of boundary points stored for two-parameter plotting.
#' Ignored for regions with more than two selected coefficients.
#'
#' @details
#' Confidence regions are computed on the MSM fitting coefficient scale. For
#' mean-difference and risk-difference estimands this is the identity scale.
#' For odds-ratio and rate-ratio estimands this is the log coefficient scale.
#' Version `0.5.0` does not geometrically transform confidence ellipsoids onto
#' nonlinear ratio display scales.
#'
#' The region is defined as
#'
#' \deqn{(\theta - \hat\theta)' \hat\Sigma^{-1} (\theta - \hat\theta) \leq \chi^2_{p, level},}
#'
#' where `p` is the number of selected coefficients and the covariance matrix is
#' the bootstrap covariance matrix restricted to those coefficients.
#'
#' Pooled multiple-imputation confidence regions are intentionally out of
#' scope for Version `0.5.0`.
#'
#' @examples
#' \dontrun{
#' dat <- sim_mixture_data(
#'   n = 400, pA = 3, pB = 3, rho_within_A = 0.3, rho_within_B = 0.3,
#'   rho_between = 0.2, psi1 = 0.5, psi2 = 0.3, psi12 = 0.2, seed = 123
#' )
#' fit <- qgcomp.glm.multi(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = dat, mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"), B = 100, seed = 13
#' )
#' region <- confregion(fit, parm = c("psi1", "psi2"))
#' region
#' plot(region)
#' }
#'
#' @export
confregion.qgcompmulti <- function(object,
                                   parm = c("psi1", "psi2"),
                                   level = 0.95,
                                   method = "bootstrap_chisq",
                                   npoints = 200L,
                                   ...) {
  validate_qgcompmulti(object)
  qgcompmulti_validate_conf_level(level)
  qgcompmulti_validate_region_method(method)
  npoints <- qgcompmulti_validate_region_npoints(npoints)
  coef_names <- names(object$results$coefficients)
  parm <- qgcompmulti_resolve_region_parm(parm = parm, coef_names = coef_names)
  center <- object$results$coefficients[parm]
  covariance <- object$results$vcov[parm, parm, drop = FALSE]
  qgcompmulti_validate_region_covariance(covariance = covariance, parm = parm)
  df <- length(parm)
  threshold <- stats::qchisq(level, df = df)
  coefficient_labels <- object$labels$coefficient_labels[parm]
  names(coefficient_labels) <- parm
  plot_data <- qgcompmulti_region_plot_data(
    center = center,
    covariance = covariance,
    threshold = threshold,
    npoints = npoints
  )
  new_qgcompmulti_region(
    parm = parm,
    center = center,
    covariance = covariance,
    level = level,
    method = method,
    threshold = threshold,
    df = df,
    coefficient_labels = coefficient_labels,
    estimand_scale = object$analysis$estimand_scale,
    msm_fitting_scale = object$analysis$msm_fitting_scale,
    plot_data = plot_data
  )
}
#' @rdname confregion
#' @export
confregion.qgcompmulti_mi <- function(object, ...) {
  stop(
    paste(
      "Confidence regions are not supported for pooled qgcompmulti",
      "multiple-imputation objects in Version 0.5.0."
    ),
    call. = FALSE
  )
}
#' @export
print.qgcompmulti_region <- function(x, ...) {
  validate_qgcompmulti_region(x)
  cat("qgcomp.multi confidence region\n")
  cat("  Method: ", qgcompmulti_region_method_label(x$method), "\n", sep = "")
  cat("  Level: ", formatC(100 * x$level, format = "f", digits = 1), "%\n", sep = "")
  cat("  Degrees of freedom: ", x$df, "\n", sep = "")
  cat("  Chi-squared cutoff: ", formatC(x$threshold, format = "f", digits = 4), "\n", sep = "")
  cat("  Estimand scale: ", qgcompmulti_estimand_label(x$estimand_scale), "\n", sep = "")
  cat("  Fitting scale: ", qgcompmulti_msm_fitting_scale_label(x$msm_fitting_scale), "\n", sep = "")
  cat("  Parameters: ", paste(unname(x$coefficient_labels), collapse = ", "), "\n", sep = "")
  cat("  Note: ", qgcompmulti_region_scale_note(x$estimand_scale, x$msm_fitting_scale), "\n", sep = "")
  cat("\nCenter (fitting scale):\n")
  print(x$center)
  cat("\nBootstrap covariance (fitting scale):\n")
  print(x$covariance)
  invisible(x)
}
