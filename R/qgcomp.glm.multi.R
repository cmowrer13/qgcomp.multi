#' Fit a quantile g-computation model for two exposure mixtures
#'
#' Fits a two-mixture g-computation model with an optional interaction between
#' the mixture intervention variables. The function fits an outcome regression,
#' computes predicted potential outcomes under joint interventions on the two
#' mixtures, and then fits a marginal structural model (MSM) to summarize the
#' resulting intervention-response surface. The current interface is designed for
#' analyses with exactly two mixtures, supplied through `mix1` and `mix2`.
#' Uncertainty is estimated using a nonparametric bootstrap.
#'
#' The fitted MSM has the form
#'
#' \deqn{
#' E[Y^{x(q1), w(q2)}] = \psi_0 + \psi_1 q_1 + \psi_2 q_2 + \psi_{12} q_1 q_2
#' }
#'
#' when `interaction = TRUE`. When `q` is an integer, `q1` and `q2` index
#' quantized joint intervention levels for mixtures 1 and 2. When `q = NULL`,
#' the same linear MSM is fit over original-scale intervention values using a
#' `3 x 3` grid formed from the pooled 25th, 50th, and 75th percentile values
#' within each mixture. In that setting, every component in a mixture is set to
#' the same pooled mixture-specific value under a given intervention.
#'
#' @param f A model formula for the outcome regression. The formula should
#' include the outcome and any baseline covariates. Mixture variables listed
#' in `mix1` and `mix2` should also appear in the formula if they are to be
#' included in the outcome model. The `mix1` and `mix2` arguments define how
#' exposures are grouped into mixtures, but they do not add variables to the
#' formula automatically.
#' @param data A data frame containing the outcome, exposure variables, and
#' any covariates in the model.
#' @param mix1 A character vector giving the names of the variables in the
#' first exposure mixture.
#' @param mix2 A character vector giving the names of the variables in the
#' second exposure mixture.
#' @param interaction Logical; if `TRUE`, the package includes an interaction
#' term in both the outcome regression and the fitted MSM. In the current
#' implementation, the outcome model is augmented with a product between
#' the sums of the components in `mix1` and `mix2`, and the MSM includes the
#' `psi1 * psi2` interaction term. If `FALSE`, both models are fit without that
#' interaction.
#' @param family A GLM family object (e.g., `gaussian()`, `binomial()`,
#' `poisson()`) specifying the outcome model.
#' @param q Integer greater than or equal to 2 giving the number of quantiles
#' used to discretize the exposure variables, or `NULL` to skip quantization and
#' fit the outcome model on the original exposure scale. When `q = NULL`, the
#' fit-time intervention grid is defined by the pooled 25th, 50th, and 75th
#' percentile values within each mixture, and under each intervention every
#' component in a mixture is set to the same pooled mixture-specific value.
#' @param centering Character string controlling how the marginal structural
#' model intervention variables are coded when `q = NULL`. Must be one of
#' `"none"` or `"median"`. With `"none"`, the MSM uses the raw intervention
#' values. With `"median"`, the MSM uses intervention values centered at the
#' pooled median within each mixture. Centering affects only the MSM fit, not
#' the outcome regression. This argument is ignored when `q` is numeric.
#' @param B Integer greater than or equal to 2 giving the number of bootstrap
#' replications used for standard error estimation.
#' @param id Optional character string giving the name of a cluster identifier
#' variable. If supplied, bootstrap resampling is performed at the cluster
#' level rather than the observation level.
#' @param MCsize Optional integer controlling the Monte Carlo sample size used
#' in the g-computation step. If equal to `nrow(data)`, the empirical
#' covariate distribution is fully enumerated. Smaller values approximate the
#' marginalization step using a random subsample, which can reduce computation
#' time in large datasets. When `id` is supplied and `MCsize < nrow(data)`,
#' the approximation is implemented by sampling `MCsize` clusters with
#' replacement.
#' @param seed Optional integer random seed used to make bootstrap resampling
#' and any Monte Carlo subsampling reproducible. If `NULL`, the current RNG
#' state is used and not modified by `qgcomp.glm.multi()`. When supplied, the
#' same inputs and the same seed reproduce the same fitted object and
#' inferential results.
#'
#' @return An object of class `"qgcompmulti"` representing the fitted
#' two-mixture quantile g-computation model. Major components include:
#' \describe{
#'   \item{`data_info`}{Outcome name, sample-size metadata, and indicators for
#'   quantization and clustered fitting.}
#'   \item{`mixtures`}{The mixture definitions, quantization setting `q`, and
#'   original-scale centering choice when `q = NULL`.}
#'   \item{`analysis`}{Model settings such as the GLM family, interaction
#'   status, bootstrap count, cluster identifier, Monte Carlo size, and any
#'   supplied random seed.}
#'   \item{`fits`}{The fitted outcome regression and marginal structural model
#'   (MSM) objects.}
#'   \item{`prediction`}{Stored fit-time prediction objects for later
#'   prediction, plotting, and diagnostic methods, including the intervention
#'   grid, the MSM-coded grid, the exact counterfactual surface, the
#'   corresponding MSM fitted surface, and a comparison object on the common
#'   fit-time grid.}
#'   \item{`bootstrap`}{Retained bootstrap coefficient draws, bootstrap
#'   replication counts and lightweight failure data.}
#'   \item{`results`}{The MSM coefficient vector, standard errors, covariance
#'   matrix, and coefficient table.}
#'   \item{`labels`}{Internal coefficient names and human-readable labels used by
#'   the print and summary methods.}
#' }
#'
#' @details
#' This function extends quantile g-computation to two exposure mixtures by
#' evaluating predicted outcomes over a two-dimensional intervention grid.
#' For each bootstrap replication, the observed data are resampled, exposures
#' are either quantized or left on their original scale depending on `q`, the
#' outcome model is fit, and predicted potential outcomes are computed under
#' joint interventions on the two mixtures. A marginal structural model is then
#' fit to those predicted counterfactual means to obtain the reported mixture
#' effect estimates.
#'
#' When `q` is an integer, the intervention grid is `0, 1, ..., q - 1` for each
#' mixture. In that setting, `psi1` and `psi2` are interpreted as the change in
#' the marginal mean outcome associated with simultaneously increasing every
#' component in the corresponding mixture by one quantile, holding the other
#' mixture intervention level fixed at the lowest quantile.
#'
#' When `q = NULL`, the exposure variables are left on their original analysis
#' scales. The intervention grid is then defined by assigning every component in
#' a mixture to a common pooled percentile value from that mixture. In this
#' setting, the MSM coefficients are defined with respect to that original-scale
#' intervention coding, so their units depend on the measurement scale of the
#' underlying exposures. If `centering = "median"`, the MSM is fit on centered
#' intervention values and the intercept corresponds to the pooled-median
#' intervention for both mixtures.
#'
#' If `interaction = TRUE`, the current implementation adds an interaction term
#' to both the outcome regression and the MSM. This means the MSM is summarizing
#' a counterfactual surface implied by an interacting outcome model rather than
#' simply adding an interaction at the final summary step.
#'
#' The outcome model is fit using `glm()`, so this function can be used with
#' Gaussian, binomial, Poisson, and other generalized linear models supported by
#' the supplied formula and family specification. Predicted potential outcomes
#' are computed on the response scale and the MSM is fit using an identity link,
#' regardless of the outcome-model link function.
#'
#' Interpretation therefore depends on the outcome type:
#' \itemize{
#'   \item For continuous outcomes, parameters represent mean differences.
#'   \item For binary outcomes, parameters represent differences in predicted
#'   probabilities, that is, risk differences.
#'   \item For count outcomes, parameters represent differences in expected counts.
#' }
#'
#' For causal interpretation, the usual identifying conditions for
#' g-computation still apply: consistency, conditional exchangeability,
#' positivity for the interventions under study, and adequate specification of
#' the outcome model. The support diagnostic can help users inspect the
#' intervention grid, but it should not be read as a full positivity proof.
#'
#' @examples
#' dat <- sim_mixture_data(
#'   n = 500,
#'   pA = 4,
#'   pB = 4,
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
#' fit <- qgcomp.glm.multi(
#'   f = Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3", "X4"),
#'   mix2 = c("W1", "W2", "W3", "W4"),
#'   interaction = TRUE,
#'   q = 4,
#'   B = 100,
#'   MCsize = nrow(dat),
#'   seed = 13
#' )
#'
#' fit
#' summary(fit)
#' coef(fit)
#' confint(fit)
#'
#' # Public prediction and plotting workflow
#' \dontrun{
#' predict(fit)
#' predict(
#'   fit,
#'   type = "msm_contrast",
#'   from = c(psi1 = 0, psi2 = 0),
#'   to = c(psi1 = 3, psi2 = 3),
#'   interval = TRUE
#' )
#' plot(fit)
#' plot(fit, style = "contour")
#'
#' # Diagnostics
#' support(fit)
#' diagnostics(fit, type = "bootstrap")
#' adequacy(fit)
#'
#' # Sensitivity helpers
#' mcsize_sensitivity(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   MCsize_values = c(250, 500),
#'   q = 4,
#'   B = 100,
#'   seed = 13
#' )
#' }
#'
#' dat_cont <- sim_mixture_data(
#'   n = 500,
#'   pA = 4,
#'   pB = 4,
#'   rho_within_A = 0.3,
#'   rho_within_B = 0.3,
#'   rho_between = 0.2,
#'   psi1 = 0.5,
#'   psi2 = 0.3,
#'   psi12 = 0.2,
#'   return_quantized = FALSE,
#'   seed = 321
#' )
#'
#' fit_cont <- qgcomp.glm.multi(
#'   f = Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + C,
#'   data = dat_cont,
#'   mix1 = c("X1", "X2", "X3", "X4"),
#'   mix2 = c("W1", "W2", "W3", "W4"),
#'   interaction = TRUE,
#'   q = NULL,
#'   centering = "median",
#'   B = 100,
#'   MCsize = nrow(dat_cont),
#'   seed = 13
#' )
#'
#' fit_cont
#' summary(fit_cont)
#' coef(fit_cont)
#' confint(fit_cont)
#'
#' @importFrom stats as.formula coef gaussian glm predict quantile rnorm setNames terms vcov
#'
#' @export

qgcomp.glm.multi <- function(f,
                             data,
                             mix1,
                             mix2,
                             interaction = TRUE,
                             family = gaussian(),
                             q = 4,
                             centering = "none",
                             B = 200,
                             id = NULL,
                             MCsize = nrow(data),
                             seed = NULL) {
  call <- match.call()

  validate_qgcomp_multi_inputs(
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
    seed = seed
  )

  qgcompmulti_with_seed(
    seed,
    {
      data_q <- quantize_mixtures(data, mix1, mix2, q)
      full_fit <- qgcompmulti_msm_fit(
        f = f,
        data = data_q,
        mix1 = mix1,
        mix2 = mix2,
        interaction = interaction,
        family = family,
        q = q,
        centering = centering,
        id = id,
        MCsize = MCsize
      )
      labels <- build_qgcompmulti_labels(interaction = interaction)
      coefs <- full_fit$coefficients
      psi_hat <- vector("list", B)
      failure_log <- vector("list", 0L)
      for (b in seq_len(B)) {
        data_b <- qgcompmulti_resample_data(
          data = data_q,
          id = id
        )
        boot_fit <- tryCatch(
          qgcompmulti_msm_fit(
            f = f,
            data = data_b,
            mix1 = mix1,
            mix2 = mix2,
            interaction = interaction,
            family = family,
            q = q,
            centering = centering,
            id = id,
            MCsize = MCsize
          ),
          error = function(e) e
        )
        if (inherits(boot_fit, "error")) {
          failure_log[[length(failure_log) + 1L]] <- data.frame(
            replicate = b,
            message = conditionMessage(boot_fit),
            row.names = NULL,
            check.names = FALSE
          )
          next
        }
        psi_hat[[b]] <- boot_fit$coefficients
      }
      success_idx <- which(!vapply(psi_hat, is.null, logical(1)))
      if (length(success_idx) == 0L) {
        stop("All bootstrap replications failed.", call. = FALSE)
      }
      psi_hat <- do.call(
        rbind,
        lapply(
          psi_hat[success_idx],
          function(x) matrix(x, nrow = 1L, dimnames = list(NULL, names(x)))
        )
      )
      std.err <- build_qgcompmulti_std_error(
        coef_draws = psi_hat,
        coef_names = names(coefs)
      )
      varcov <- build_qgcompmulti_vcov(
        coef_draws = psi_hat,
        coef_names = names(coefs)
      )
      data_info <- build_qgcompmulti_data_info(
        data = data,
        formula = f,
        q = q,
        id = id,
        n_used = full_fit$n_used
      )
      mixtures <- build_qgcompmulti_mixtures(
        mix1 = mix1,
        mix2 = mix2,
        q = q,
        centering = centering
      )
      analysis <- build_qgcompmulti_analysis(
        interaction = interaction,
        family = family,
        B = B,
        id = id,
        MCsize = MCsize,
        seed = seed
      )
      fits <- build_qgcompmulti_fits(full_fit)
      prediction <- build_qgcompmulti_prediction(
        intervention_grid = full_fit$intervention_grid,
        msm_grid = full_fit$msm_grid,
        counterfactual_surface = full_fit$counterfactual_surface,
        msm_surface = full_fit$msm_surface,
        surface_comparison = full_fit$surface_comparison
      )
      failure_log_df <- if (length(failure_log) == 0L) {
        NULL
      } else {
        do.call(rbind, failure_log)
      }
      bootstrap <- build_qgcompmulti_bootstrap(
        coef_draws = psi_hat,
        B_requested = B,
        failure_log = failure_log_df
      )
      results <- build_qgcompmulti_results(coefs, std.err, varcov)
      fit <- new_qgcompmulti(
        call = call,
        formula = f,
        data_info = data_info,
        mixtures = mixtures,
        analysis = analysis,
        fits = fits,
        prediction = prediction,
        bootstrap = bootstrap,
        results = results,
        labels = labels
      )
      fit
    }
  )
}
