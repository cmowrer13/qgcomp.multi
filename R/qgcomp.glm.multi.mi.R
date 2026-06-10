#' Fit a pooled multiple-imputation quantile g-computation model
#'
#' Fits one ordinary `qgcomp.glm.multi()` model per completed imputation and
#' pools the marginal structural model (MSM) coefficient inference internally
#' using Rubin's rules. This is the primary native multiple-imputation workflow
#' for `qgcomp.multi` in Version `0.4.0`.
#'
#' @param f A model formula for the outcome regression.
#' @param data Either a `mice::mids` object or a non-empty list of completed
#'   data frames.
#' @param mix1 A character vector giving the names of the variables in the
#'   first exposure mixture.
#' @param mix2 A character vector giving the names of the variables in the
#'   second exposure mixture.
#' @param interaction Logical; if `TRUE`, include the mixture interaction in
#'   each per-imputation fit and in the pooled MSM.
#' @param family A GLM family object for the per-imputation outcome model.
#' @param q Integer greater than or equal to `2` giving the number of quantiles
#'   used to discretize the exposure variables, or `NULL` to fit on the original
#'   exposure scale.
#' @param centering Character string controlling how the MSM intervention
#'   variables are coded when `q = NULL`. Must be one of `"none"` or
#'   `"median"`.
#' @param B Integer greater than or equal to `2` giving the number of bootstrap
#'   replications used within each completed-data fit.
#' @param id Optional character string giving the name of a cluster identifier
#'   variable for clustered bootstrap resampling.
#' @param MCsize Optional integer controlling the Monte Carlo sample size used
#'   within each completed-data fit. If `NULL`, the wrapper uses the completed
#'   dataset row count.
#' @param seed Optional integer master seed. When supplied, the wrapper
#'   deterministically derives one distinct fit-specific seed per imputation and
#'   stores both the master seed and the fit-specific seeds in the pooled
#'   result.
#' @param keep_fits Logical; if `TRUE`, retain the full per-imputation
#'   `qgcompmulti` fitted objects inside the pooled result. Defaults to `FALSE`
#'   to keep the pooled object lean.
#' @param progress Logical; if `TRUE`, pass the existing serial bootstrap
#'   progress display through to each per-imputation fit.
#'
#' @return An object of class `"qgcompmulti_mi"` containing pooled MSM
#'   coefficients, pooled covariance and standard errors, Rubin pooling
#'   components, MI metadata, and optionally the retained per-imputation fitted
#'   objects.
#'
#' @details
#' This wrapper is an orchestrator over repeated ordinary fits rather than a
#' second modeling engine. After normalizing the imputations into a completed
#' data-frame list, the function validates that the completed datasets are
#' structurally compatible, fits one `qgcomp.glm.multi()` model per imputation,
#' and pools the resulting MSM coefficient vectors and covariance matrices using
#' the internal Rubin-pooling helpers defined for the `"qgcompmulti_mi"` class.
#'
#' Pooled results in `0.4.0` are intentionally focused on MSM inference and a
#' compact summary layer. The pooled object does not yet try to replicate the
#' full prediction and diagnostics surface of ordinary single-fit
#' `qgcompmulti` objects.
#'
#' When `data` is a `mids` object, the optional `mice` package must be
#' installed. Inputs supplied as a plain list must already be fully completed:
#' required analysis variables cannot still contain missing values.
#'
#' @examples
#' dat <- sim_mixture_data(
#'   n = 120,
#'   pA = 3,
#'   pB = 3,
#'   rho_within_A = 0.3,
#'   rho_within_B = 0.3,
#'   rho_between = 0.2,
#'   psi1 = 0.5,
#'   psi2 = 0.3,
#'   psi12 = 0.2,
#'   return_quantized = FALSE,
#'   seed = 123
#' )
#'
#' imp_list <- lapply(
#'   seq_len(3),
#'   function(i) {
#'     dat_i <- dat
#'     dat_i$X1 <- dat_i$X1 + rnorm(nrow(dat_i), sd = 0.02 * i)
#'     dat_i$W2 <- dat_i$W2 + rnorm(nrow(dat_i), sd = 0.02 * i)
#'     dat_i
#'   }
#' )
#'
#' fit_mi <- qgcomp.glm.multi.mi(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = imp_list,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   interaction = TRUE,
#'   q = 4,
#'   B = 5,
#'   seed = 13
#' )
#'
#' fit_mi$results$coef_table
#'
#' @export
qgcomp.glm.multi.mi <- function(f,
                                data,
                                mix1,
                                mix2,
                                interaction = TRUE,
                                family = gaussian(),
                                q = 4,
                                centering = "none",
                                B = 200,
                                id = NULL,
                                MCsize = NULL,
                                seed = NULL,
                                keep_fits = FALSE,
                                progress = FALSE) {
  call <- match.call()

  validate_qgcomp_multi_mi_inputs(
    f = f,
    data = data,
    mix1 = mix1,
    mix2 = mix2,
    interaction = interaction,
    family = family,
    q = q,
    centering = centering,
    id = id,
    MCsize = MCsize,
    B = B,
    seed = seed,
    keep_fits = keep_fits,
    progress = progress
  )

  input_type <- qgcompmulti_mi_input_type(data)
  completed <- qgcompmulti_mi_completed_data(data, input_type = input_type)
  mcsize_value <- if (is.null(MCsize)) nrow(completed[[1]]) else as.integer(MCsize)

  qgcompmulti_mi_validate_completed_data(
    completed = completed,
    f = f,
    mix1 = mix1,
    mix2 = mix2,
    interaction = interaction,
    family = family,
    q = q,
    centering = centering,
    id = id,
    MCsize = mcsize_value,
    B = B,
    progress = progress
  )

  fit_seeds <- qgcompmulti_mi_fit_seeds(length(completed), seed = seed)
  fits <- qgcompmulti_mi_fit_models(
    completed = completed,
    f = f,
    mix1 = mix1,
    mix2 = mix2,
    interaction = interaction,
    family = family,
    q = q,
    centering = centering,
    B = B,
    id = id,
    MCsize = mcsize_value,
    fit_seeds = fit_seeds,
    progress = progress
  )

  qgcompmulti_pool_mi_fits(
    fits = fits,
    input_type = input_type,
    keep_fits = keep_fits,
    seed = seed,
    fit_seeds = if (is.null(seed)) NULL else as.integer(unlist(fit_seeds, use.names = FALSE)),
    call = call,
    formula = f
  )
}
