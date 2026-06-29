#' @keywords internal
#' @noRd
pooled_mix_quantiles <- function(data, vars, probs = c(0.25, 0.50, 0.75)){

  as.numeric(quantile(
    unlist(data[, vars, drop = FALSE], use.names = FALSE),
    probs = probs,
    na.rm = TRUE
  ))
}

#' Fit the core outcome-model and MSM components for two-mixture g-computation
#'
#' Fits the outcome regression, computes predicted potential outcomes under joint
#' interventions on two exposure mixtures, and estimates the marginal structural
#' model (MSM) coefficients that summarize the resulting two-dimensional
#' intervention-response surface.
#'
#' This is a lower-level fitting helper used internally by
#' [qgcomp.glm.multi()]. It is exported because it can be useful for method
#' development, testing, and direct inspection of the fitted outcome model, MSM,
#' and stored fit-time surfaces. Most users will want [qgcomp.glm.multi()]
#' instead.
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
#' @param interaction Logical; if `TRUE`, the package includes an interaction
#' term in both the outcome regression and the fitted MSM. In the current
#' implementation, the outcome model is augmented with a cross-product between
#' the sums of the components in `mix1` and `mix2`, and the MSM includes the
#' `psi1 * psi2` interaction term. If `FALSE`, both models are fit without that
#' interaction.
#' @param family A GLM family object (e.g., `gaussian()`, `binomial()`,
#' `poisson()`) specifying the outcome model.
#' @param estimand_scale Optional character string naming the fitted MSM
#' estimand scale. Supported values depend on `family` and its link. If `NULL`,
#' the Version `0.5.0` defaults are used.
#' @param q Integer greater than or equal to 2 giving the number of quantiles
#' used to discretize the exposure variables, or `NULL` to skip quantization and
#' fit the outcome model on the original exposure scale. When `q = NULL`, the
#' fit-time intervention grid is defined by the pooled 25th, 50th, and 75th
#' percentile values within each mixture, and under each intervention every
#' component in a mixture is set to the same pooled mixture-specific value.
#' @param centering Character string controlling how the marginal structural
#' model intervention variables are coded when `q = NULL`. Must be one of `"none"`
#' or `"median"`. Centering affects only the MSM predictors and does not change
#' the outcome regression fit. This argument is ignored when `q` is numeric.
#' @param id Optional character string giving the name of a cluster identifier
#' variable. If supplied, Monte Carlo subsampling is performed at the cluster
#' level rather than the observation level.
#' @param MCsize Optional integer controlling the Monte Carlo sample size used
#' to approximate the marginalization step in g-computation. If `MCsize` is
#' greater than or equal to the current analysis sample size, all observations
#' are used. Smaller values compute predicted outcomes over a random subsample
#' drawn from the empirical distribution. When `id` is supplied and
#' `MCsize < nrow(data)`, the approximation is implemented by sampling `MCsize`
#' clusters with replacement.
#' @param seed Optional integer random seed used to make Monte Carlo
#' subsampling reproducible when `MCsize < nrow(data)`. If `NULL`, the current
#' RNG state is used and not modified by `qgcompmulti_msm_fit()`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`outcome_fit`}{The fitted outcome regression model object.}
#'   \item{`msm_fit`}{The fitted marginal structural model object.}
#'   \item{`coefficients`}{A named vector of MSM coefficients.}
#'   \item{`n_used`}{The number of observations used in the g-computation
#'   prediction step.}
#'   \item{`intervention_grid`}{The fit-time intervention grid on the
#'   intervention-value scale.}
#'   \item{`msm_grid`}{The corresponding fit-time grid on the MSM coding scale.}
#'   \item{`counterfactual_surface`}{The exact fit-time counterfactual mean
#'   surface implied by the fitted outcome model.}
#'   \item{`msm_surface`}{The fitted MSM surface evaluated on the common
#'   fit-time grid.}
#'   \item{`surface_comparison`}{A direct exact-versus-MSM comparison object on
#'   the common fit-time grid.}
#'   \item{`counterfactual_surface_target`}{The transformed fit-time MSM target
#'   surface used in the MSM fit.}
#'   \item{`msm_surface_target`}{The fitted MSM target-scale surface evaluated
#'   on the common fit-time grid.}
#'   \item{`surface_comparison_target`}{A direct exact-versus-MSM comparison
#'   object on the target scale.}
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
#' outcomes are computed on the response scale. The marginal structural model
#' is then fit either on that response-scale surface or on a transformed target
#' surface implied by `estimand_scale`.
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
                                estimand_scale = NULL,
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
    estimand_scale = estimand_scale,
    q = q,
    centering = centering,
    id = id,
    MCsize = MCsize,
    seed = seed
  )

  qgcompmulti_with_seed(
    seed,
    {
      estimand <- qgcompmulti_resolve_estimand_spec(
        family = family,
        estimand_scale = estimand_scale,
        mode = "planned"
      )
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
      if (MCsize >= nrow(data)) {
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
      counterfactual_surface_target <- qgcompmulti_counterfactual_surface_target(
        counterfactual_surface = counterfactual_surface,
        msm_fitting_scale = estimand$msm_fitting_scale
      )
      msmdat <- data.frame(
        Ya   = rep(counterfactual_surface_target$exact_target, each = n_used),
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
      msm_surface_target <- qgcompmulti_msm_surface_target(
        msm_fit = msmfit,
        intervention_grid = intervention_grid,
        msm_grid = msm_grid
      )
      msm_surface <- qgcompmulti_msm_surface(
        msm_surface_target = msm_surface_target,
        msm_fitting_scale = estimand$msm_fitting_scale
      )
      surface_comparison <- qgcompmulti_surface_comparison(
        counterfactual_surface = counterfactual_surface,
        msm_surface = msm_surface
      )
      surface_comparison_target <- qgcompmulti_surface_comparison_target(
        counterfactual_surface_target = counterfactual_surface_target,
        msm_surface_target = msm_surface_target
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
        surface_comparison = surface_comparison,
        counterfactual_surface_target = counterfactual_surface_target,
        msm_surface_target = msm_surface_target,
        surface_comparison_target = surface_comparison_target
      )
    }
  )
}
