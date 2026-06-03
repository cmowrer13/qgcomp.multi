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
build_qgcompmulti_prediction <- function(intervention_grid = NULL,
                                         msm_grid = NULL,
                                         counterfactual_surface = NULL,
                                         msm_surface = NULL,
                                         surface_comparison = NULL) {
  list(
    intervention_grid = intervention_grid,
    msm_grid = msm_grid,
    counterfactual_surface = counterfactual_surface,
    msm_surface = msm_surface,
    surface_comparison = surface_comparison
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_validate_prediction <- function(prediction) {
  data_fields <- c(
    "intervention_grid",
    "msm_grid",
    "counterfactual_surface",
    "msm_surface",
    "surface_comparison"
  )
  for (field in data_fields) {
    value <- prediction[[field]]
    if (!is.null(value) && !is.data.frame(value)) {
      stop(
        sprintf("`prediction$%s` must be `NULL` or a data frame.", field),
        call. = FALSE
      )
    }
  }
  invisible(prediction)
}
