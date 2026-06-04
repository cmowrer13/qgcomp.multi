#' Predict from a qgcompmulti fit
#'
#' Generates prediction objects from a fitted [qgcomp.glm.multi()] model.
#' The default prediction target is the fitted marginal structural model
#' (MSM) surface. Exact fit-time surface extraction and exact arbitrary
#' prediction on user-supplied data are also supported through explicit
#' `type` values.
#'
#' @param object A fitted `"qgcompmulti"` object.
#' @param type Character string specifying the prediction target. Supported
#' values are `"msm"`, `"msm_point"`, `"msm_contrast"`, `"exact"`, and
#' `"exact_contrast"`.
#' @param grid Optional data frame with columns `psi1` and `psi2` giving a
#' user-specified prediction grid.
#' @param at Optional named numeric vector or list with entries `psi1` and
#' `psi2` giving a single intervention regime.
#' @param from,to Optional named numeric vectors or lists with entries `psi1`
#' and `psi2` defining the source and target regimes for a direct contrast.
#' @param data Optional data frame used for exact arbitrary prediction or exact
#' arbitrary contrasts. Exact arbitrary prediction requires explicit `data`.
#' @param interval Logical; if `TRUE`, returns bootstrap percentile intervals
#' when supported for the requested prediction type.
#' @param level Confidence level used when `interval = TRUE`.
#' @param ... Unused.
#'
#' @return A structured list describing the requested prediction target. The
#' returned object includes metadata identifying the prediction type, grid type,
#' grid scale, estimand scale, the prediction estimates, and any interval
#' information.
#'
#' @export
predict.qgcompmulti <- function(object,
                                type = c("msm", "msm_point", "msm_contrast", "exact", "exact_contrast"),
                                grid = NULL,
                                at = NULL,
                                from = NULL,
                                to = NULL,
                                data = NULL,
                                interval = FALSE,
                                level = 0.95,
                                ...) {
  validate_qgcompmulti(object)
  type <- match.arg(type)
  validate_predict_qgcompmulti_inputs(
    type = type,
    grid = grid,
    at = at,
    from = from,
    to = to,
    data = data,
    interval = interval,
    level = level
  )
  level_arg <- if (isTRUE(interval)) level else NULL
  switch(
    type,
    msm = qgcompmulti_predict_msm_surface(
      object = object,
      grid = grid,
      level = level_arg
    ),
    msm_point = qgcompmulti_predict_msm_point(
      object = object,
      psi1 = at[["psi1"]],
      psi2 = at[["psi2"]],
      level = level_arg
    ),
    msm_contrast = qgcompmulti_predict_msm_contrast(
      object = object,
      from = from,
      to = to,
      level = level_arg
    ),
    exact = if (is.null(data) && is.null(grid) && is.null(at)) {
      qgcompmulti_exact_fit_surface(object)
    } else {
      qgcompmulti_exact_predict_data(
        object = object,
        data = data,
        grid = grid,
        at = at
      )
    },
    exact_contrast = qgcompmulti_exact_contrast_data(
      object = object,
      data = data,
      from = from,
      to = to
    )
  )
}
