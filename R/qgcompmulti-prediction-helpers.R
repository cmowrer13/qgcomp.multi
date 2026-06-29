# ------------------------------------------------------------------------------
# Prediction helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_resample_data <- function(data, id = NULL) {
  if (is.null(id)) {
    idx <- sample(seq_len(nrow(data)), replace = TRUE)
    return(data[idx, , drop = FALSE])
  }
  ids <- data[[id]]
  unique_ids <- unique(ids)
  samp_ids <- sample(unique_ids, length(unique_ids), replace = TRUE)
  split_data <- split(data, ids)
  do.call(
    rbind,
    split_data[as.character(samp_ids)]
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(eval.parent(substitute(expr)))
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit(
    {
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  set.seed(as.integer(seed))
  eval.parent(substitute(expr))
}
#' @keywords internal
#' @noRd
build_qgcompmulti_fit_time_grids <- function(data,
                                             mix1,
                                             mix2,
                                             q,
                                             centering = "none") {
  if (is.null(q)) {
    mix1_values <- pooled_mix_quantiles(data, mix1)
    mix2_values <- pooled_mix_quantiles(data, mix2)
    intervention_grid <- expand.grid(
      psi1 = mix1_values,
      psi2 = mix2_values
    )
    msm_grid <- intervention_grid
    if (centering == "median") {
      msm_grid$psi1 <- msm_grid$psi1 - mix1_values[2]
      msm_grid$psi2 <- msm_grid$psi2 - mix2_values[2]
    }
  } else {
    intervention_grid <- expand.grid(
      psi1 = 0:(q - 1),
      psi2 = 0:(q - 1)
    )
    msm_grid <- intervention_grid
  }
  intervention_grid <- cbind(
    grid_id = seq_len(nrow(intervention_grid)),
    intervention_grid,
    row.names = NULL
  )
  msm_grid <- cbind(
    grid_id = seq_len(nrow(msm_grid)),
    msm_grid,
    row.names = NULL
  )
  list(
    intervention_grid = intervention_grid,
    msm_grid = msm_grid
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_counterfactual_surface <- function(outcome_fit,
                                               nd,
                                               mix1,
                                               mix2,
                                               intervention_grid,
                                               msm_grid) {
  pred_fun <- function(psi1, psi2, nd) {
    nd2 <- nd
    nd2[, mix1] <- psi1
    nd2[, mix2] <- psi2
    predict(outcome_fit, newdata = nd2, type = "response")
  }
  predmat <- Map(
    pred_fun,
    intervention_grid$psi1,
    intervention_grid$psi2,
    MoreArgs = list(nd = nd)
  )
  exact_mean <- vapply(predmat, mean, numeric(1))
  data.frame(
    grid_id = intervention_grid$grid_id,
    intervention_psi1 = intervention_grid$psi1,
    intervention_psi2 = intervention_grid$psi2,
    msm_psi1 = msm_grid$psi1,
    msm_psi2 = msm_grid$psi2,
    exact_mean = exact_mean,
    row.names = NULL,
    check.names = FALSE
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_counterfactual_surface_target <- function(counterfactual_surface,
                                                      msm_fitting_scale) {
  exact_target <- qgcompmulti_transform_msm_surface(
    values = counterfactual_surface$exact_mean,
    msm_fitting_scale = msm_fitting_scale,
    direction = "to_fitting"
  )

  data.frame(
    grid_id = counterfactual_surface$grid_id,
    intervention_psi1 = counterfactual_surface$intervention_psi1,
    intervention_psi2 = counterfactual_surface$intervention_psi2,
    msm_psi1 = counterfactual_surface$msm_psi1,
    msm_psi2 = counterfactual_surface$msm_psi2,
    exact_target = exact_target,
    row.names = NULL,
    check.names = FALSE
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_msm_surface_target <- function(msm_fit,
                                           intervention_grid,
                                           msm_grid) {
  msm_target <- as.numeric(
    predict(msm_fit, newdata = msm_grid[, c("psi1", "psi2"), drop = FALSE], type = "response")
  )
  data.frame(
    grid_id = intervention_grid$grid_id,
    intervention_psi1 = intervention_grid$psi1,
    intervention_psi2 = intervention_grid$psi2,
    msm_psi1 = msm_grid$psi1,
    msm_psi2 = msm_grid$psi2,
    msm_target = msm_target,
    row.names = NULL,
    check.names = FALSE
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_msm_surface <- function(msm_surface_target,
                                    msm_fitting_scale) {
  msm_mean <- qgcompmulti_transform_msm_surface(
    values = msm_surface_target$msm_target,
    msm_fitting_scale = msm_fitting_scale,
    direction = "to_response"
  )

  data.frame(
    grid_id = msm_surface_target$grid_id,
    intervention_psi1 = msm_surface_target$intervention_psi1,
    intervention_psi2 = msm_surface_target$intervention_psi2,
    msm_psi1 = msm_surface_target$msm_psi1,
    msm_psi2 = msm_surface_target$msm_psi2,
    msm_mean = msm_mean,
    row.names = NULL,
    check.names = FALSE
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_surface_comparison <- function(counterfactual_surface,
                                           msm_surface) {
  if (!identical(counterfactual_surface$grid_id, msm_surface$grid_id)) {
    stop("Fit-time surfaces must align on `grid_id`.", call. = FALSE)
  }
  data.frame(
    grid_id = counterfactual_surface$grid_id,
    intervention_psi1 = counterfactual_surface$intervention_psi1,
    intervention_psi2 = counterfactual_surface$intervention_psi2,
    msm_psi1 = counterfactual_surface$msm_psi1,
    msm_psi2 = counterfactual_surface$msm_psi2,
    exact_mean = counterfactual_surface$exact_mean,
    msm_mean = msm_surface$msm_mean,
    residual = counterfactual_surface$exact_mean - msm_surface$msm_mean,
    row.names = NULL,
    check.names = FALSE
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_surface_comparison_target <- function(counterfactual_surface_target,
                                                  msm_surface_target) {
  if (!identical(counterfactual_surface_target$grid_id, msm_surface_target$grid_id)) {
    stop("Fit-time target surfaces must align on `grid_id`.", call. = FALSE)
  }

  data.frame(
    grid_id = counterfactual_surface_target$grid_id,
    intervention_psi1 = counterfactual_surface_target$intervention_psi1,
    intervention_psi2 = counterfactual_surface_target$intervention_psi2,
    msm_psi1 = counterfactual_surface_target$msm_psi1,
    msm_psi2 = counterfactual_surface_target$msm_psi2,
    exact_target = counterfactual_surface_target$exact_target,
    msm_target = msm_surface_target$msm_target,
    residual_target = counterfactual_surface_target$exact_target - msm_surface_target$msm_target,
    row.names = NULL,
    check.names = FALSE
  )
}
#' @keywords internal
#' @noRd
build_qgcompmulti_prediction <- function(intervention_grid = NULL,
                                         msm_grid = NULL,
                                         counterfactual_surface = NULL,
                                         msm_surface = NULL,
                                         surface_comparison = NULL,
                                         counterfactual_surface_target = NULL,
                                         msm_surface_target = NULL,
                                         surface_comparison_target = NULL) {
  list(
    intervention_grid = intervention_grid,
    msm_grid = msm_grid,
    counterfactual_surface = counterfactual_surface,
    msm_surface = msm_surface,
    surface_comparison = surface_comparison,
    counterfactual_surface_target = counterfactual_surface_target,
    msm_surface_target = msm_surface_target,
    surface_comparison_target = surface_comparison_target
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_validate_prediction <- function(prediction) {
  field_types <- list(
    intervention_grid = c("grid_id", "psi1", "psi2"),
    msm_grid = c("grid_id", "psi1", "psi2"),
    counterfactual_surface = c(
      "grid_id",
      "intervention_psi1",
      "intervention_psi2",
      "msm_psi1",
      "msm_psi2",
      "exact_mean"
    ),
    msm_surface = c(
      "grid_id",
      "intervention_psi1",
      "intervention_psi2",
      "msm_psi1",
      "msm_psi2",
      "msm_mean"
    ),
    surface_comparison = c(
      "grid_id",
      "intervention_psi1",
      "intervention_psi2",
      "msm_psi1",
      "msm_psi2",
      "exact_mean",
      "msm_mean",
      "residual"
    ),
    counterfactual_surface_target = c(
      "grid_id",
      "intervention_psi1",
      "intervention_psi2",
      "msm_psi1",
      "msm_psi2",
      "exact_target"
    ),
    msm_surface_target = c(
      "grid_id",
      "intervention_psi1",
      "intervention_psi2",
      "msm_psi1",
      "msm_psi2",
      "msm_target"
    ),
    surface_comparison_target = c(
      "grid_id",
      "intervention_psi1",
      "intervention_psi2",
      "msm_psi1",
      "msm_psi2",
      "exact_target",
      "msm_target",
      "residual_target"
    )
  )
  for (field in names(field_types)) {
    value <- prediction[[field]]
    if (!is.null(value) && !is.data.frame(value)) {
      stop(
        sprintf("`prediction$%s` must be `NULL` or a data frame.", field),
        call. = FALSE
      )
    }
    if (!is.null(value)) {
      missing_cols <- setdiff(field_types[[field]], names(value))
      if (length(missing_cols) > 0L) {
        stop(
          sprintf(
            "`prediction$%s` is missing required columns: %s",
            field,
            paste(missing_cols, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
  }
  if (!is.null(prediction$intervention_grid) &&
      !is.null(prediction$msm_grid) &&
      nrow(prediction$intervention_grid) != nrow(prediction$msm_grid)) {
    stop("`prediction$intervention_grid` and `prediction$msm_grid` must have the same number of rows.", call. = FALSE)
  }
  if (!is.null(prediction$counterfactual_surface) &&
      !is.null(prediction$msm_surface) &&
      nrow(prediction$counterfactual_surface) != nrow(prediction$msm_surface)) {
    stop("`prediction$counterfactual_surface` and `prediction$msm_surface` must have the same number of rows.", call. = FALSE)
  }
  if (!is.null(prediction$surface_comparison) &&
      !is.null(prediction$counterfactual_surface) &&
      nrow(prediction$surface_comparison) != nrow(prediction$counterfactual_surface)) {
    stop("`prediction$surface_comparison` must align with the stored fit-time surfaces.", call. = FALSE)
  }
  if (!is.null(prediction$counterfactual_surface_target) &&
      !is.null(prediction$msm_surface_target) &&
      nrow(prediction$counterfactual_surface_target) != nrow(prediction$msm_surface_target)) {
    stop(
      "`prediction$counterfactual_surface_target` and `prediction$msm_surface_target` must have the same number of rows.",
      call. = FALSE
    )
  }
  if (!is.null(prediction$surface_comparison_target) &&
      !is.null(prediction$counterfactual_surface_target) &&
      nrow(prediction$surface_comparison_target) != nrow(prediction$counterfactual_surface_target)) {
    stop(
      "`prediction$surface_comparison_target` must align with the stored target-scale surfaces.",
      call. = FALSE
    )
  }
  invisible(prediction)
}
