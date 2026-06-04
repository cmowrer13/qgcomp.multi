#' @keywords internal
#' @noRd
pooled_mix_quantiles <- function(data, vars, probs = c(0.25, 0.50, 0.75)){

  as.numeric(quantile(
    unlist(data[, vars, drop = FALSE], use.names = FALSE),
    probs = probs,
    na.rm = TRUE
  ))
}

#' Fit the marginal structural model for two-mixture quantile g-computation
#'
#' Fits the outcome regression, computes predicted potential outcomes under
#' joint interventions on two exposure mixtures, and estimates the marginal
#' structural model (MSM) coefficients that summarize the resulting
#' two-dimensional dose-response surface. This function is primarily an
#' internal helper used by [qgcomp.glm.multi()] and will not typically
#' need to be called directly by users.
#'
#' Given a fitted outcome model, the function evaluates predicted outcomes over
#' either a grid of uniform mixture quantiles or, when `q = NULL`, a 3 x 3 grid
#' formed from the pooled 25th, 50th, and 75th percentiles of each mixture, and
#' regresses these predictions on the mixture intervention levels to obtain MSM
#' parameters. When `interaction = TRUE`, the fitted MSM has the form
#'
#' \deqn{
#' E[Y^{x(q1), w(q2)}] = \psi_1 q_1 + \psi_2 q_2 + \psi_{12} q_1 q_2
#' }
#'
#' where `psi1` and `psi2` represent the main effects of increasing the
#' intervention level for mixtures 1 and 2 by one unit on the MSM scale, and
#' `psi12` represents their interaction on that same scale. For quantized
#' analyses this is a one-quantile increase; for `q = NULL` it is a one-unit
#' increase in the original-scale intervention coding.
#'
#' When `q = NULL` and `centering = "median"`, the MSM is fit using intervention
#' values centered at the pooled median of each mixture.
#'
#' @param f A model formula for the outcome regression. The formula should
#' include the outcome and any baseline covariates. Mixture variables listed
#' in `mix1` and `mix2` should also appear in the formula if they are to be
#' included in the outcome model.
#' @param data A data frame containing the outcome, exposure variables, and
#' any covariates in the model.
#' @param mix1 A character vector giving the names of the variables in the
#' first exposure mixture.
#' @param mix2 A character vector giving the names of the variables in the
#' second exposure mixture.
#' @param interaction Logical; if `TRUE`, includes an interaction term between
#' the two mixture indices in the marginal structural model. If `FALSE`, only
#' main effects are estimated.
#' @param family A GLM family object (e.g., `gaussian()`, `binomial()`,
#' `poisson()`) specifying the outcome model.
#' @param q Integer greater than or equal to 2 giving the number of quantiles
#' used to discretize the exposure variables and define the intervention grid,
#' or `NULL` to skip quantization and define a 3 x 3 intervention grid using
#' pooled 25th, 50th, and 75th percentiles for each mixture.
#' @param centering Character string controlling how the marginal structural
#' model intervention variables are coded when `q = NULL`. Must be one of `"none"`
#' or `"median"`. Centering affects only the MSM predictors and does not change
#' the outcome regression fit. This argument is ignored when `q` is numeric.
#' @param id Optional character string giving the name of a cluster identifier
#' variable. If supplied, Monte Carlo subsampling is performed at the cluster
#' level rather than the observation level.
#' @param MCsize Optional integer controlling the Monte Carlo sample size used
#' to approximate the marginalization step in g-computation. If equal to
#' `nrow(data)`, all observations are used. Smaller values compute predicted
#' outcomes over a random subsample drawn from the empirical distribution. When
#' `id` is supplied and `MCsize < nrow(data)`, the approximation is implemented
#' by sampling `MCsize` clusters with replacement.
#' @param seed Optional integer random seed used to make Monte Carlo
#' subsampling reproducible when `MCsize < nrow(data)`. If `NULL`, the current
#' RNG state is used and not modified by `qgcompmulti_msm_fit()`.
#'
#' @return A list with components:
#' \describe{
#' \item{outcome_fit}{The fitted outcome regression model object.}
#' \item{msm_fit}{The fitted marginal structural model object.}
#' \item{coefficients}{A named vector of MSM coefficients.}
#' \item{n_used}{The number of observations used in the MSM pseudo-dataset
#' prediction step.}
#' }
#'
#' @details
#' This function carries out the core g-computation step for the two-mixture
#' extension. It first fits a generalized linear model for the outcome using
#' `glm()`. It then constructs a grid of intervention levels over the two
#' mixtures, replaces the observed mixture values with the intervention values,
#' and computes predicted outcomes for each observation (or for a Monte Carlo
#' subsample if `MCsize < nrow(data)`). These predicted outcomes are stacked
#' into a pseudo-dataset and used to fit a marginal structural model that
#' summarizes the dose-response surface. Note that predicted potential
#' outcomes are computed on the response scale and the marginal
#' structural model is fit using an identity link, regardless of the outcome
#' model.
#'
#' When `q` is an integer, the intervention grid is `0, 1, ..., q - 1` for each
#' mixture, corresponding to simultaneous quantile increases in all components
#' of that mixture.
#'
#' When `q = NULL`, each mixture is instead set to common pooled percentile
#' values on the original exposure scale. The resulting MSM coefficients are
#' therefore scale-dependent and should be interpreted in the units of the
#' underlying exposures. If `centering = "median"`, the intercept corresponds
#' to the pooled-median intervention for both mixtures.
#'
#' Because the number of intervention combinations grows as `q^2`, this step
#' can become computationally expensive for large datasets. The `MCsize`
#' argument reduces this burden by approximating the empirical covariate
#' distribution using a random subset while leaving the outcome model fit
#' unchanged.
#'
#' @examples
#' dat <- sim_mixture_data(
#'   n = 500,
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
#' qgcompmulti_msm_fit(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   interaction = TRUE,
#'   q = 4,
#'   MCsize = nrow(dat),
#'   seed = 13
#' )
#'
#' @export

qgcompmulti_msm_fit <- function(f,
                                data,
                                mix1,
                                mix2,
                                interaction = TRUE,
                                family = gaussian(),
                                q = 4,
                                centering = "none",
                                id = NULL,
                                MCsize = nrow(data),
                                seed = NULL) {

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
    seed = seed
  )

  qgcompmulti_with_seed(
    seed,
    {
      response <- deparse(f[[2]])
      rhs_terms <- attr(terms(f), "term.labels")
      rhs <- paste(rhs_terms, collapse = " + ")
      if (interaction) {
        m1 <- paste(mix1, collapse = " + ")
        m2 <- paste(mix2, collapse = " + ")
        interaction_term <- paste0("I((", m1, ") * (", m2, "))")
        rhs <- paste(rhs, interaction_term, sep = " + ")
      }
      outcome_formula <- as.formula(paste(response, "~", rhs))
      outcome_fit <- glm(outcome_formula, data = data, family = family)
      grids <- build_qgcompmulti_fit_time_grids(
        data = data,
        mix1 = mix1,
        mix2 = mix2,
        q = q,
        centering = centering
      )
      intervention_grid <- grids$intervention_grid
      msm_grid <- grids$msm_grid
      if (MCsize == nrow(data)) {
        nd <- data
      } else if (is.null(id)) {
        samp_idx <- sample(seq_len(nrow(data)), MCsize, replace = TRUE)
        nd <- data[samp_idx, , drop = FALSE]
      } else {
        ids <- data[[id]]
        unique_ids <- unique(ids)
        samp_ids <- sample(unique_ids, MCsize, replace = TRUE)
        split_data <- split(data, ids)
        nd <- do.call(
          rbind,
          split_data[as.character(samp_ids)]
        )
      }
      n_used <- nrow(nd)
      counterfactual_surface <- qgcompmulti_counterfactual_surface(
        outcome_fit = outcome_fit,
        nd = nd,
        mix1 = mix1,
        mix2 = mix2,
        intervention_grid = intervention_grid,
        msm_grid = msm_grid
      )
      msmdat <- data.frame(
        Ya   = rep(counterfactual_surface$exact_mean, each = n_used),
        psi1 = rep(msm_grid$psi1, each = n_used),
        psi2 = rep(msm_grid$psi2, each = n_used)
      )
      if (interaction) {
        msm_formula <- Ya ~ psi1 * psi2
      } else {
        msm_formula <- Ya ~ psi1 + psi2
      }
      msmfit <- glm(
        msm_formula,
        data = msmdat
      )
      coefficients <- qgcompmulti_extract_msm_coefficients(
        msm_fit = msmfit,
        interaction = interaction
      )
      msm_surface <- qgcompmulti_msm_surface(
        msm_fit = msmfit,
        intervention_grid = intervention_grid,
        msm_grid = msm_grid
      )
      surface_comparison <- qgcompmulti_surface_comparison(
        counterfactual_surface = counterfactual_surface,
        msm_surface = msm_surface
      )
      list(
        outcome_fit = outcome_fit,
        msm_fit = msmfit,
        coefficients = coefficients,
        n_used = n_used,
        intervention_grid = intervention_grid,
        msm_grid = msm_grid,
        counterfactual_surface = counterfactual_surface,
        msm_surface = msm_surface,
        surface_comparison = surface_comparison
      )
    }
  )
}
