# ------------------------------------------------------------------------------
# Internal prediction and contrast builders
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_predict_msm_surface <- function(object,
                                            grid = NULL,
                                            level = NULL) {
  validate_qgcompmulti(object)
  grid_info <- qgcompmulti_build_msm_grid(object, grid = grid)
  grid_df <- grid_info$grid
  estimates <- data.frame(
    grid_id = grid_df$grid_id,
    psi1 = grid_df$psi1,
    psi2 = grid_df$psi2,
    estimate = as.numeric(
      stats::predict(
        object$fits$msm_fit,
        newdata = grid_df[, c("psi1", "psi2"), drop = FALSE],
        type = "response"
      )
    ),
    row.names = NULL
  )
  intervals <- NULL
  interval_type <- NULL
  uncertainty_source <- NULL
  if (!is.null(level)) {
    draws <- qgcompmulti_bootstrap_msm_fitted_draws(object, grid_df)
    intervals <- cbind(
      estimates[, c("grid_id", "psi1", "psi2"), drop = FALSE],
      qgcompmulti_bootstrap_interval(draws, level = level)
    )
    interval_type <- "bootstrap_percentile"
    uncertainty_source <- "stored_bootstrap_draws"
  }
  build_qgcompmulti_prediction_result(
    prediction_type = "msm_surface",
    grid_type = grid_info$grid_type,
    grid_scale = "msm",
    estimates = estimates,
    intervals = intervals,
    interval_type = interval_type,
    uncertainty_source = uncertainty_source
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_predict_msm_point <- function(object,
                                          psi1,
                                          psi2,
                                          level = NULL) {
  point <- c(psi1 = psi1, psi2 = psi2)
  qgcompmulti_validate_regime(point, arg = "point")
  point_result <- qgcompmulti_predict_msm_surface(
    object = object,
    level = level,
    grid = data.frame(psi1 = psi1, psi2 = psi2, row.names = NULL)
  )
  point_result$prediction_type <- "msm_point"
  point_result$grid_type <- "point_regime"
  point_result
}
#' @keywords internal
#' @noRd
qgcompmulti_predict_msm_contrast <- function(object,
                                             from,
                                             to,
                                             level = NULL) {
  qgcompmulti_validate_regime(from, arg = "from")
  qgcompmulti_validate_regime(to, arg = "to")
  from_grid <- qgcompmulti_build_msm_grid(object, at = from)$grid
  to_grid <- qgcompmulti_build_msm_grid(object, at = to)$grid
  from_est <- qgcompmulti_predict_msm_surface(object, grid = from_grid, level = NULL)$estimates$estimate
  to_est <- qgcompmulti_predict_msm_surface(object, grid = to_grid, level = NULL)$estimates$estimate
  estimates <- data.frame(
    from_psi1 = unname(from[["psi1"]]),
    from_psi2 = unname(from[["psi2"]]),
    to_psi1 = unname(to[["psi1"]]),
    to_psi2 = unname(to[["psi2"]]),
    estimate = as.numeric(to_est - from_est),
    row.names = NULL
  )
  intervals <- NULL
  interval_type <- NULL
  uncertainty_source <- NULL
  if (!is.null(level)) {
    from_draws <- qgcompmulti_bootstrap_msm_fitted_draws(object, from_grid)
    to_draws <- qgcompmulti_bootstrap_msm_fitted_draws(object, to_grid)
    contrast_draws <- to_draws - from_draws
    interval_vals <- qgcompmulti_bootstrap_interval(contrast_draws, level = level)
    intervals <- cbind(estimates[, c("from_psi1", "from_psi2", "to_psi1", "to_psi2"), drop = FALSE], interval_vals)
    interval_type <- "bootstrap_percentile"
    uncertainty_source <- "stored_bootstrap_draws"
  }
  build_qgcompmulti_prediction_result(
    prediction_type = "msm_contrast",
    grid_type = "pairwise_regime",
    grid_scale = "msm",
    estimates = estimates,
    intervals = intervals,
    interval_type = interval_type,
    uncertainty_source = uncertainty_source,
    contrast = TRUE
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_exact_fit_surface <- function(object) {
  validate_qgcompmulti(object)
  build_qgcompmulti_prediction_result(
    prediction_type = "exact_fit_surface",
    grid_type = "stored_fit_grid",
    grid_scale = "intervention",
    estimates = object$prediction$counterfactual_surface,
    data_supplied = FALSE
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_exact_predict_data <- function(object,
                                           data,
                                           grid = NULL,
                                           at = NULL) {
  validate_qgcompmulti(object)
  if (missing(data) || is.null(data)) {
    stop("Exact arbitrary prediction requires explicit `data`.", call. = FALSE)
  }
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("`data` must be a non-empty data frame.", call. = FALSE)
  }
  required_vars <- unique(c(all.vars(object$formula), object$mixtures$mix1, object$mixtures$mix2))
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      sprintf(
        "The following variables are missing from `data`: %s",
        paste(missing_vars, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  grid_info <- qgcompmulti_build_intervention_grid(object, grid = grid, at = at)
  grid_df <- grid_info$grid
  msm_grid <- qgcompmulti_map_intervention_to_msm(object, grid_df)
  estimates <- qgcompmulti_counterfactual_surface(
    outcome_fit = object$fits$outcome_fit,
    nd = data,
    mix1 = object$mixtures$mix1,
    mix2 = object$mixtures$mix2,
    intervention_grid = grid_df,
    msm_grid = msm_grid
  )
  build_qgcompmulti_prediction_result(
    prediction_type = "exact_arbitrary",
    grid_type = grid_info$grid_type,
    grid_scale = "intervention",
    estimates = estimates,
    data_supplied = TRUE
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_exact_contrast_data <- function(object,
                                            data,
                                            from,
                                            to) {
  qgcompmulti_validate_regime(from, arg = "from")
  qgcompmulti_validate_regime(to, arg = "to")
  from_est <- qgcompmulti_exact_predict_data(object, data = data, at = from)$estimates$exact_mean
  to_est <- qgcompmulti_exact_predict_data(object, data = data, at = to)$estimates$exact_mean
  estimates <- data.frame(
    from_psi1 = unname(from[["psi1"]]),
    from_psi2 = unname(from[["psi2"]]),
    to_psi1 = unname(to[["psi1"]]),
    to_psi2 = unname(to[["psi2"]]),
    estimate = as.numeric(to_est - from_est),
    row.names = NULL
  )
  build_qgcompmulti_prediction_result(
    prediction_type = "exact_contrast",
    grid_type = "pairwise_regime",
    grid_scale = "intervention",
    estimates = estimates,
    data_supplied = TRUE,
    contrast = TRUE
  )
}
