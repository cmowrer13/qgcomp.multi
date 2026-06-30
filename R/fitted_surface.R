# ------------------------------------------------------------------------------
# Internal prediction and contrast builders
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_predict_msm_target <- function(object,
                                           grid) {
  validate_qgcompmulti(object)
  grid <- qgcompmulti_validate_msm_grid(grid)
  as.numeric(
    stats::predict(
      object$fits$msm_fit,
      newdata = grid[, c("psi1", "psi2"), drop = FALSE],
      type = "response"
    )
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_predict_msm_response <- function(object,
                                             grid) {
  target <- qgcompmulti_predict_msm_target(object = object, grid = grid)
  qgcompmulti_transform_msm_surface(
    values = target,
    msm_fitting_scale = object$analysis$msm_fitting_scale,
    direction = "to_response"
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_predict_msm_surface <- function(object,
                                            grid = NULL,
                                            level = NULL,
                                            method = NULL) {
  validate_qgcompmulti(object)
  grid_info <- qgcompmulti_build_msm_grid(object, grid = grid)
  grid_df <- grid_info$grid
  estimates <- data.frame(
    grid_id = grid_df$grid_id,
    psi1 = grid_df$psi1,
    psi2 = grid_df$psi2,
    estimate = qgcompmulti_predict_msm_response(object = object, grid = grid_df),
    row.names = NULL
  )
  intervals <- NULL
  interval_type <- NULL
  uncertainty_source <- NULL
  if (!is.null(level)) {
    method <- qgcompmulti_resolve_prediction_interval_method(
      method = method,
      default_method = object$analysis$default_interval_method
    )
    target_estimates <- qgcompmulti_predict_msm_target(object = object, grid = grid_df)
    target_draws <- qgcompmulti_bootstrap_msm_fitted_draws(object, grid_df)
    interval_vals <- qgcompmulti_bootstrap_interval(
      draws = target_draws,
      estimates = target_estimates,
      level = level,
      method = method
    )
    interval_vals[] <- qgcompmulti_transform_msm_surface(
      values = as.matrix(interval_vals),
      msm_fitting_scale = object$analysis$msm_fitting_scale,
      direction = "to_response"
    )
    intervals <- cbind(
      estimates[, c("grid_id", "psi1", "psi2"), drop = FALSE],
      interval_vals
    )
    interval_type <- qgcompmulti_prediction_interval_type_for_method(method)
    uncertainty_source <- "stored_bootstrap_draws"
  }
  build_qgcompmulti_prediction_result(
    prediction_type = "msm_surface",
    grid_type = grid_info$grid_type,
    grid_scale = "msm",
    estimand_scale = "response",
    estimate_scale = "response",
    fit_estimand_scale = object$analysis$estimand_scale,
    msm_fitting_scale = object$analysis$msm_fitting_scale,
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
                                          level = NULL,
                                          method = NULL) {
  point <- c(psi1 = psi1, psi2 = psi2)
  qgcompmulti_validate_regime(point, arg = "point")
  point_result <- qgcompmulti_predict_msm_surface(
    object = object,
    level = level,
    method = method,
    grid = data.frame(psi1 = psi1, psi2 = psi2, row.names = NULL)
  )
  point_result$prediction_type <- "msm_point"
  point_result$grid_type <- "point_regime"
  point_result
}
#' @keywords internal
#' @noRd
qgcompmulti_msm_contrast_estimate <- function(object,
                                              from_grid,
                                              to_grid,
                                              contrast_scale = c("response", "estimand")) {
  contrast_scale <- match.arg(contrast_scale)
  from_target <- qgcompmulti_predict_msm_target(object = object, grid = from_grid)
  to_target <- qgcompmulti_predict_msm_target(object = object, grid = to_grid)
  if (identical(contrast_scale, "response")) {
    from_response <- qgcompmulti_transform_msm_surface(
      values = from_target,
      msm_fitting_scale = object$analysis$msm_fitting_scale,
      direction = "to_response"
    )
    to_response <- qgcompmulti_transform_msm_surface(
      values = to_target,
      msm_fitting_scale = object$analysis$msm_fitting_scale,
      direction = "to_response"
    )
    return(as.numeric(to_response - from_response))
  }
  target_difference <- to_target - from_target
  as.numeric(qgcompmulti_transform_msm_coefficients(
    values = target_difference,
    estimand_scale = object$analysis$estimand_scale,
    direction = "to_display"
  ))
}
#' @keywords internal
#' @noRd
qgcompmulti_msm_contrast_draws <- function(object,
                                           from_grid,
                                           to_grid,
                                           contrast_scale = c("response", "estimand")) {
  contrast_scale <- match.arg(contrast_scale)
  from_draws <- qgcompmulti_bootstrap_msm_fitted_draws(object, from_grid)
  to_draws <- qgcompmulti_bootstrap_msm_fitted_draws(object, to_grid)
  if (identical(contrast_scale, "response")) {
    from_response_draws <- qgcompmulti_transform_msm_fitted_draws(
      draws = from_draws,
      msm_fitting_scale = object$analysis$msm_fitting_scale,
      direction = "to_response"
    )
    to_response_draws <- qgcompmulti_transform_msm_fitted_draws(
      draws = to_draws,
      msm_fitting_scale = object$analysis$msm_fitting_scale,
      direction = "to_response"
    )
    return(to_response_draws - from_response_draws)
  }
  target_difference_draws <- to_draws - from_draws
  qgcompmulti_transform_msm_coefficients(
    values = target_difference_draws,
    estimand_scale = object$analysis$estimand_scale,
    direction = "to_display"
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_predict_msm_contrast <- function(object,
                                             from,
                                             to,
                                             contrast_scale = c("response", "estimand"),
                                             level = NULL,
                                             method = NULL) {
  contrast_scale <- match.arg(contrast_scale)
  qgcompmulti_validate_regime(from, arg = "from")
  qgcompmulti_validate_regime(to, arg = "to")
  from_grid <- qgcompmulti_build_msm_grid(object, at = from)$grid
  to_grid <- qgcompmulti_build_msm_grid(object, at = to)$grid
  estimate <- qgcompmulti_msm_contrast_estimate(
    object = object,
    from_grid = from_grid,
    to_grid = to_grid,
    contrast_scale = contrast_scale
  )
  estimates <- data.frame(
    from_psi1 = unname(from[["psi1"]]),
    from_psi2 = unname(from[["psi2"]]),
    to_psi1 = unname(to[["psi1"]]),
    to_psi2 = unname(to[["psi2"]]),
    estimate = estimate,
    row.names = NULL
  )
  intervals <- NULL
  interval_type <- NULL
  uncertainty_source <- NULL
  if (!is.null(level)) {
    method <- qgcompmulti_resolve_prediction_interval_method(
      method = method,
      default_method = object$analysis$default_interval_method
    )
    if (identical(contrast_scale, "estimand")) {
      from_target <- qgcompmulti_predict_msm_target(object = object, grid = from_grid)
      to_target <- qgcompmulti_predict_msm_target(object = object, grid = to_grid)
      from_draws <- qgcompmulti_bootstrap_msm_fitted_draws(object, from_grid)
      to_draws <- qgcompmulti_bootstrap_msm_fitted_draws(object, to_grid)
      interval_vals <- qgcompmulti_bootstrap_interval(
        draws = to_draws - from_draws,
        estimates = to_target - from_target,
        level = level,
        method = method
      )
      interval_vals[] <- qgcompmulti_transform_msm_coefficients(
        values = as.matrix(interval_vals),
        estimand_scale = object$analysis$estimand_scale,
        direction = "to_display"
      )
    } else {
      contrast_draws <- qgcompmulti_msm_contrast_draws(
        object = object,
        from_grid = from_grid,
        to_grid = to_grid,
        contrast_scale = contrast_scale
      )
      interval_vals <- qgcompmulti_bootstrap_interval(
        draws = contrast_draws,
        estimates = estimate,
        level = level,
        method = method
      )
    }
    intervals <- cbind(estimates[, c("from_psi1", "from_psi2", "to_psi1", "to_psi2"), drop = FALSE], interval_vals)
    interval_type <- qgcompmulti_prediction_interval_type_for_method(method)
    uncertainty_source <- "stored_bootstrap_draws"
  }
  build_qgcompmulti_prediction_result(
    prediction_type = "msm_contrast",
    grid_type = "pairwise_regime",
    grid_scale = "msm",
    estimand_scale = if (identical(contrast_scale, "estimand")) object$analysis$estimand_scale else "response",
    estimate_scale = contrast_scale,
    fit_estimand_scale = object$analysis$estimand_scale,
    msm_fitting_scale = object$analysis$msm_fitting_scale,
    contrast_scale = contrast_scale,
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

#' Extract stored fit-time surfaces from a qgcompmulti fit
#'
#' Returns the fit-time exact counterfactual surface, the fit-time MSM fitted
#' surface, or the stored exact-versus-MSM comparison object retained in a
#' fitted [qgcomp.glm.multi()] object.
#'
#' @param object A fitted `"qgcompmulti"` object.
#' @param type Character string indicating which stored surface object to
#' return. Supported values are `"all"`, `"exact"`, `"msm"`, and
#' `"comparison"`.
#' @param ... Unused.
#'
#' @return Either a list of stored surfaces or one stored surface data frame,
#' depending on `type`.
#'
#' @export
fitted_surface <- function(object, type = c("all", "exact", "msm", "comparison"), ...) {
  UseMethod("fitted_surface")
}
#' @export
fitted_surface.qgcompmulti <- function(object,
                                       type = c("all", "exact", "msm", "comparison"),
                                       ...) {
  validate_qgcompmulti(object)
  type <- match.arg(type)
  switch(
    type,
    all = list(
      counterfactual_surface = object$prediction$counterfactual_surface,
      msm_surface = object$prediction$msm_surface,
      surface_comparison = object$prediction$surface_comparison
    ),
    exact = object$prediction$counterfactual_surface,
    msm = object$prediction$msm_surface,
    comparison = object$prediction$surface_comparison
  )
}
