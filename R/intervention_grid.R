# ------------------------------------------------------------------------------
# Internal prediction-grid helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_validate_prediction_type <- function(prediction_type) {
  allowed <- c(
    "msm_surface",
    "msm_point",
    "msm_contrast",
    "exact_fit_surface",
    "exact_arbitrary",
    "exact_contrast"
  )
  if (!is.character(prediction_type) ||
      length(prediction_type) != 1L ||
      is.na(prediction_type) ||
      !prediction_type %in% allowed) {
    stop(
      sprintf(
        "`prediction_type` must be one of: %s.",
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(prediction_type)
}
#' @keywords internal
#' @noRd
qgcompmulti_validate_interval_type <- function(interval_type) {
  if (is.null(interval_type)) {
    return(invisible(NULL))
  }
  if (!is.character(interval_type) ||
      length(interval_type) != 1L ||
      is.na(interval_type) ||
      !interval_type %in% "bootstrap_percentile") {
    stop(
      "`interval_type` must be `NULL` or \"bootstrap_percentile\".",
      call. = FALSE
    )
  }
  invisible(interval_type)
}
#' @keywords internal
#' @noRd
qgcompmulti_validate_regime <- function(at, arg = "at") {
  if (is.null(at)) {
    return(invisible(NULL))
  }
  if (is.list(at) && !is.data.frame(at)) {
    at <- unlist(at, use.names = TRUE)
  }
  if (!is.numeric(at) || length(at) != 2L || anyNA(at)) {
    stop(
      sprintf("`%s` must be a named numeric vector or list with `psi1` and `psi2`.", arg),
      call. = FALSE
    )
  }
  if (!identical(sort(names(at)), c("psi1", "psi2"))) {
    stop(
      sprintf("`%s` must have names `psi1` and `psi2`.", arg),
      call. = FALSE
    )
  }
  if (any(!is.finite(at))) {
    stop(sprintf("`%s` must contain finite values.", arg), call. = FALSE)
  }
  invisible(at)
}
#' @keywords internal
#' @noRd
qgcompmulti_standardize_grid <- function(grid) {
  grid <- as.data.frame(grid)
  if (!all(c("psi1", "psi2") %in% names(grid))) {
    stop("Prediction grids must contain `psi1` and `psi2` columns.", call. = FALSE)
  }
  if (!all(vapply(grid[c("psi1", "psi2")], is.numeric, logical(1)))) {
    stop("Prediction grid columns `psi1` and `psi2` must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(grid$psi1)) || any(!is.finite(grid$psi2))) {
    stop("Prediction grid columns `psi1` and `psi2` must be finite.", call. = FALSE)
  }
  grid <- grid[, intersect(c("grid_id", "psi1", "psi2"), names(grid)), drop = FALSE]
  if (!"grid_id" %in% names(grid)) {
    grid <- cbind(grid_id = seq_len(nrow(grid)), grid, row.names = NULL)
  } else {
    grid <- grid[, c("grid_id", "psi1", "psi2"), drop = FALSE]
  }
  row.names(grid) <- NULL
  grid
}
#' @keywords internal
#' @noRd
qgcompmulti_validate_msm_grid <- function(grid) {
  grid <- qgcompmulti_standardize_grid(grid)
  invisible(grid)
}
#' @keywords internal
#' @noRd
qgcompmulti_validate_intervention_grid <- function(grid) {
  grid <- qgcompmulti_standardize_grid(grid)
  invisible(grid)
}
#' @keywords internal
#' @noRd
qgcompmulti_build_msm_grid <- function(object, grid = NULL, at = NULL) {
  validate_qgcompmulti(object)
  if (!is.null(grid) && !is.null(at)) {
    stop("Supply only one of `grid` or `at`.", call. = FALSE)
  }
  if (!is.null(at)) {
    qgcompmulti_validate_regime(at, arg = "at")
    point_grid <- data.frame(
      grid_id = 1L,
      psi1 = unname(at[["psi1"]]),
      psi2 = unname(at[["psi2"]]),
      row.names = NULL
    )
    return(list(grid = point_grid, grid_type = "point_regime"))
  }
  if (is.null(grid)) {
    return(list(
      grid = object$prediction$msm_grid,
      grid_type = "stored_fit_grid"
    ))
  }
  list(
    grid = qgcompmulti_validate_msm_grid(grid),
    grid_type = "user_grid"
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_build_intervention_grid <- function(object, grid = NULL, at = NULL) {
  validate_qgcompmulti(object)
  if (!is.null(grid) && !is.null(at)) {
    stop("Supply only one of `grid` or `at`.", call. = FALSE)
  }
  if (!is.null(at)) {
    qgcompmulti_validate_regime(at, arg = "at")
    point_grid <- data.frame(
      grid_id = 1L,
      psi1 = unname(at[["psi1"]]),
      psi2 = unname(at[["psi2"]]),
      row.names = NULL
    )
    return(list(grid = point_grid, grid_type = "point_regime"))
  }
  if (is.null(grid)) {
    return(list(
      grid = object$prediction$intervention_grid,
      grid_type = "stored_fit_grid"
    ))
  }
  list(
    grid = qgcompmulti_validate_intervention_grid(grid),
    grid_type = "user_grid"
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_map_intervention_to_msm <- function(object, intervention_grid) {
  validate_qgcompmulti(object)
  intervention_grid <- qgcompmulti_validate_intervention_grid(intervention_grid)
  stored_intervention <- object$prediction$intervention_grid
  stored_msm <- object$prediction$msm_grid
  if (is.null(stored_intervention) || is.null(stored_msm)) {
    return(intervention_grid)
  }
  shift1 <- unique(round(stored_intervention$psi1 - stored_msm$psi1, 10))
  shift2 <- unique(round(stored_intervention$psi2 - stored_msm$psi2, 10))
  if (length(shift1) != 1L || length(shift2) != 1L) {
    stop("Stored fit-time grids do not imply a unique intervention-to-MSM mapping.", call. = FALSE)
  }
  data.frame(
    grid_id = intervention_grid$grid_id,
    psi1 = intervention_grid$psi1 - shift1,
    psi2 = intervention_grid$psi2 - shift2,
    row.names = NULL
  )
}
#' @keywords internal
#' @noRd
build_qgcompmulti_prediction_result <- function(prediction_type,
                                                grid_type,
                                                grid_scale,
                                                estimand_scale = "response",
                                                estimates,
                                                intervals = NULL,
                                                interval_type = NULL,
                                                uncertainty_source = NULL,
                                                data_supplied = FALSE,
                                                contrast = FALSE) {
  qgcompmulti_validate_prediction_type(prediction_type)
  qgcompmulti_validate_interval_type(interval_type)
  if (!is.character(grid_type) || length(grid_type) != 1L || is.na(grid_type)) {
    stop("`grid_type` must be a single character string.", call. = FALSE)
  }
  if (!is.character(grid_scale) || length(grid_scale) != 1L || is.na(grid_scale)) {
    stop("`grid_scale` must be a single character string.", call. = FALSE)
  }
  if (!is.character(estimand_scale) || length(estimand_scale) != 1L || is.na(estimand_scale)) {
    stop("`estimand_scale` must be a single character string.", call. = FALSE)
  }
  if (!is.data.frame(estimates)) {
    stop("`estimates` must be a data frame.", call. = FALSE)
  }
  if (!is.null(intervals) && !is.data.frame(intervals)) {
    stop("`intervals` must be `NULL` or a data frame.", call. = FALSE)
  }
  list(
    prediction_type = prediction_type,
    grid_type = grid_type,
    grid_scale = grid_scale,
    estimand_scale = estimand_scale,
    estimates = estimates,
    intervals = intervals,
    interval_type = interval_type,
    uncertainty_source = uncertainty_source,
    data_supplied = data_supplied,
    contrast = contrast
  )
}
