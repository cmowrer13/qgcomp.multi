#' Predict from a qgcompmulti fit
#'
#' Generates prediction objects from a fitted [qgcomp.glm.multi()] model. The
#' default target is the fitted marginal structural model (MSM) surface. Exact
#' fit-time surface extraction and exact arbitrary prediction on user-supplied
#' data are also supported through explicit `type` values.
#'
#' In `qgcomp.multi`, it helps to keep two prediction targets separate:
#'
#' \itemize{
#'   \item \strong{MSM-based predictions}, which evaluate the fitted marginal
#'   structural model on the MSM coding scale and support bootstrap percentile
#'   and basic bootstrap intervals; and
#'   \item \strong{exact counterfactual predictions}, which come directly from
#'   the fitted outcome model under specified interventions.
#' }
#'
#' @param object A fitted `"qgcompmulti"` object.
#' @param type Character string specifying the prediction target. Supported
#' values are `"msm"`, `"msm_point"`, `"msm_contrast"`, `"exact"`, and
#' `"exact_contrast"`.
#' @param grid Optional data frame with columns `psi1` and `psi2` giving a
#' user-specified prediction grid on the MSM coding scale.
#' @param at Optional named numeric vector or list with entries `psi1` and
#' `psi2` giving a single intervention regime on the MSM coding scale.
#' @param from,to Optional named numeric vectors or lists with entries `psi1`
#' and `psi2` defining the source and target regimes for a direct contrast on
#' the MSM coding scale.
#' @param contrast_scale Character string specifying the output scale for
#' `type = "msm_contrast"`. `"response"` returns the response-scale difference
#' between the `to` and `from` regimes. `"estimand"` returns the contrast on
#' the active fitted estimand scale; for odds-ratio and rate-ratio fits this is
#' computed by differencing the MSM linear predictor and exponentiating.
#' @param data Optional data frame used for exact arbitrary prediction or exact
#' arbitrary contrasts. Exact arbitrary prediction requires explicit `data`.
#' @param interval Logical; if `TRUE`, returns bootstrap intervals when
#' supported for the requested prediction type.
#' @param level Confidence level used when `interval = TRUE`.
#' @param method Optional interval method used when `interval = TRUE` for
#'   MSM-based prediction outputs. Supported values are `"percentile"` and
#'   `"basic"`. `NULL` uses the fitted object's stored default when it is a
#'   bootstrap method, and otherwise falls back to `"percentile"` because Wald
#'   intervals are coefficient-only.
#' @param ... Unused.
#'
#' @return A structured list describing the requested prediction target. The
#' returned object includes metadata identifying the prediction type, grid type,
#' grid scale, estimand scale, the prediction estimates, and any interval
#' information.
#'
#' @details
#' The default `type = "msm"` returns fitted MSM predictions on the response
#' scale. These predictions are the natural targets for surface plotting,
#' interval construction, and direct regime contrasts because the fitted object
#' retains bootstrap coefficient draws for the MSM.
#'
#' `type = "exact"` has two different behaviors:
#' \itemize{
#'   \item If `data`, `grid`, and `at` are all omitted, the function returns the
#'   exact fit-time counterfactual surface that was stored when the model was
#'   fitted.
#'   \item If `data` is supplied together with `grid` or `at`, the function
#'   computes exact counterfactual means over the supplied covariate
#'   distribution. Explicit `data` is required because the exact counterfactual
#'   mean is defined by averaging predicted outcomes over a concrete covariate
#'   distribution.
#' }
#'
#' For MSM-based prediction, user-supplied inputs such as `grid`, `at`, `from`,
#' and `to` are interpreted on the MSM coding scale. This distinction matters
#' most when `q = NULL` and `centering = "median"`, because the plotted heatmap
#' axes are labeled on the intervention-value scale even though MSM prediction
#' inputs use centered coordinates. In that setting, `(0, 0)` corresponds to the
#' pooled median intervention for both mixtures.
#'
#' User-specified MSM grids are restricted to the stored fit-time support. The
#' function allows interpolation within that support, but it does not allow
#' arbitrary extrapolation beyond the intervention range used to fit the model.
#'
#' Public interval support is limited to MSM-based predictions and direct MSM
#' contrasts. Those intervals can use `method = "percentile"` or
#' `method = "basic"`. Exact public prediction targets do not currently return
#' intervals because they would require a separate uncertainty calculation over
#' a user-supplied covariate distribution.
#'
#' Prediction surfaces and plots are response-scale summaries even when the
#' fitted coefficient estimand is an odds ratio or rate ratio. To compare two
#' regimes on the fitted ratio scale, use `type = "msm_contrast"` with
#' `contrast_scale = "estimand"`.
#'
#' @examples
#' \dontrun{
#' dat <- sim_mixture_data(
#'   n = 400,
#'   pA = 3,
#'   pB = 3,
#'   rho_within_A = 0.3,
#'   rho_within_B = 0.3,
#'   rho_between = 0.2,
#'   psi1 = 0.5,
#'   psi2 = 0.3,
#'   psi12 = 0.2,
#'   seed = 123
#' )
#'
#' fit <- qgcomp.glm.multi(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   q = 4,
#'   B = 100,
#'   seed = 13
#' )
#'
#' # Default MSM surface prediction on the stored fit-time grid
#' pred_msm <- predict(fit)
#'
#' # MSM prediction at one regime, with bootstrap interval
#' predict(
#'   fit,
#'   type = "msm_point",
#'   at = c(psi1 = 1, psi2 = 2),
#'   interval = TRUE
#' )
#'
#' # Direct MSM contrast between two regimes
#' predict(
#'   fit,
#'   type = "msm_contrast",
#'   from = c(psi1 = 0, psi2 = 0),
#'   to = c(psi1 = 3, psi2 = 3),
#'   interval = TRUE
#' )
#'
#' # Exact fit-time surface extraction
#' predict(fit, type = "exact")
#'
#' # Exact arbitrary prediction requires explicit data because the
#' # counterfactual mean must be averaged over a concrete covariate distribution
#' predict(
#'   fit,
#'   type = "exact",
#'   data = dat,
#'   at = c(psi1 = 1, psi2 = 2)
#' )
#' }
#'
#' @export
predict.qgcompmulti <- function(object,
                                type = c("msm", "msm_point", "msm_contrast", "exact", "exact_contrast"),
                                grid = NULL,
                                at = NULL,
                                from = NULL,
                                to = NULL,
                                contrast_scale = c("response", "estimand"),
                                data = NULL,
                                interval = FALSE,
                                level = 0.95,
                                method = NULL,
                                ...) {
  validate_qgcompmulti(object)
  type <- match.arg(type)
  contrast_scale_supplied <- !missing(contrast_scale)
  if (identical(type, "msm_contrast")) {
    contrast_scale <- match.arg(contrast_scale)
  } else {
    if (contrast_scale_supplied) {
      stop(
        "`contrast_scale` is only supported when `type = \"msm_contrast\"`.",
        call. = FALSE
      )
    }
    contrast_scale <- NULL
  }
  validate_predict_qgcompmulti_inputs(
    type = type,
    grid = grid,
    at = at,
    from = from,
    to = to,
    contrast_scale = contrast_scale,
    data = data,
    interval = interval,
    level = level,
    method = method
  )
  level_arg <- if (isTRUE(interval)) level else NULL
  method_arg <- if (isTRUE(interval)) method else NULL
  switch(
    type,
    msm = qgcompmulti_predict_msm_surface(
      object = object,
      grid = grid,
      level = level_arg,
      method = method_arg
    ),
    msm_point = qgcompmulti_predict_msm_point(
      object = object,
      psi1 = at[["psi1"]],
      psi2 = at[["psi2"]],
      level = level_arg,
      method = method_arg
    ),
    msm_contrast = qgcompmulti_predict_msm_contrast(
      object = object,
      from = from,
      to = to,
      contrast_scale = contrast_scale,
      level = level_arg,
      method = method_arg
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
