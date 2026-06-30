# ------------------------------------------------------------------------------
# Internal bootstrap-based MSM prediction helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_msm_design_matrix <- function(grid, interaction = TRUE) {
  grid <- qgcompmulti_validate_msm_grid(grid)
  if (interaction) {
    design <- cbind(
      `(Intercept)` = 1,
      psi1 = grid$psi1,
      psi2 = grid$psi2,
      `psi1:psi2` = grid$psi1 * grid$psi2
    )
  } else {
    design <- cbind(
      `(Intercept)` = 1,
      psi1 = grid$psi1,
      psi2 = grid$psi2
    )
  }
  as.matrix(design)
}
#' @keywords internal
#' @noRd
qgcompmulti_bootstrap_msm_fitted_draws <- function(object, grid) {
  validate_qgcompmulti(object)
  if (is.null(object$bootstrap$coef_draws)) {
    stop("Stored bootstrap coefficient draws are required for interval estimation.", call. = FALSE)
  }
  grid <- qgcompmulti_validate_msm_grid(grid)
  interaction <- isTRUE(object$analysis$interaction)
  coef_names <- qgcompmulti_coef_names(interaction = interaction)
  coef_draws <- as.matrix(object$bootstrap$coef_draws)[, coef_names, drop = FALSE]
  design <- qgcompmulti_msm_design_matrix(grid, interaction = interaction)
  design %*% t(coef_draws)
}
#' @keywords internal
#' @noRd
qgcompmulti_transform_msm_fitted_draws <- function(draws,
                                                   msm_fitting_scale,
                                                   direction = c("to_response", "to_fitting")) {
  direction <- match.arg(direction)
  draws <- as.matrix(draws)
  transformed <- qgcompmulti_transform_msm_surface(
    values = draws,
    msm_fitting_scale = msm_fitting_scale,
    direction = direction
  )
  dim(transformed) <- dim(draws)
  dimnames(transformed) <- dimnames(draws)
  transformed
}
#' @keywords internal
#' @noRd
qgcompmulti_bootstrap_interval <- function(draws, level = 0.95) {
  qgcompmulti_validate_conf_level(level)
  draws <- as.matrix(draws)
  alpha <- (1 - level) / 2
  lower <- apply(draws, 1, stats::quantile, probs = alpha, na.rm = TRUE)
  upper <- apply(draws, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
  data.frame(
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    row.names = NULL
  )
}
