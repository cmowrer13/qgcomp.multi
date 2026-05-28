#' Fit a quantile g-computation model for two exposure mixtures
#'
#' Fits an extension of quantile g-computation for settings with two exposure
#' mixtures and an optional interaction term. The function first quantizes the
#' specified exposure variables, fits an outcome regression model, computes
#' predicted potential outcomes under joint interventions on the two mixtures,
#' and then fits a marginal structural model (MSM) to summarize the resulting
#' dose-response surface. Uncertainty is estimated using a nonparametric
#' bootstrap.
#'
#' The fitted MSM has the form
#'
#' \deqn{
#' E[Y^{x(q1), w(q2)}] = \psi_1 q_1 + \psi_2 q_2 + \psi_{12} q_1 q_2
#' }
#'
#' when `interaction = TRUE`, where `psi1` and `psi2` represent the main
#' effects of one-quantile increases in mixtures 1 and 2, respectively, and
#' `psi12` represents their interaction on the MSM scale.
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
#' @param q Integer giving the number of quantiles used to discretize the
#' exposure variables, or `NULL` to skip quantization and fit the outcome
#' model on the original exposure scale.
#' @param centering Character string controlling how the marginal structural
#' model intervention variables are coded when `q = NULL`. Must be one of
#' `"none"` or `"median"`. With `"none"`, the MSM uses the raw intervention
#' values. With `"median"`, the MSM uses intervention values centered at the
#' pooled median within each mixure. Centering affects only the MSM fit, not
#' the outcome regression.
#' @param B Integer giving the number of bootstrap replications used for
#' standard error estimation.
#' @param id Optional character string giving the name of a cluster identifier
#' variable. If supplied, bootstrap resampling is performed at the cluster
#' level rather than the observation level.
#' @param MCsize Optional integer controlling the Monte Carlo sample size used
#' in the g-computation step. If equal to `nrow(data)`, the empirical
#' covariate distribution is fully enumerated. Smaller values approximate the
#' marginalization step using a random subsample, which can reduce computation
#' time in large datasets.
#'
#' @return A list with components:
#' \describe{
#'   \item{coefs}{Estimated coefficients from the marginal structural model.}
#'   \item{std.err}{Bootstrap standard errors for the estimated coefficients.}
#'   \item{varcov}{Estimated covariance matrix for `psi1` and `psi2`.}
#'   \item{coef_table}{A coefficient table containing estimates, standard
#'   errors, z-statistics, and p-values.}
#' }
#'
#' @details
#' This function extends quantile g-computation to two exposure mixtures by
#' evaluating predicted outcomes over a two-dimensional intervention grid.
#' For each bootstrap replication, the observed data are resampled, exposures
#' are either quantized or left on their original scale depending on `q`, the
#' outcome model is fit, and predicted potential outcomesare computed under
#' either uniform quantile interventions or a 3 x 3 grid of pooled 25th/50th/75th
#' percentile interventions for the two mixtures. A marginal structural model is
#' then fit to these predicted outcomes to obtain the mixture effect estimates.
#'
#' The `MCsize` argument provides a Monte Carlo approximation to the
#' marginalization step in g-computation. Rather than averaging predictions
#' over all observations, the MSM can be fit using predictions from a random
#' subsample of size `MCsize`. This can substantially reduce computation time
#' when the sample size is large or when the intervention grid is large.
#'
#' The outcome model is fit using `glm()`, so this function can be used with
#' Gaussian, binomial, Poisson, and other generalized linear models supported
#' by the supplied formula and family specification. Note that predicted
#' potential outcomes are computed on the response scale and the marginal
#' structural model is fit using an identity link, regardless of the outcome
#' model.
#'
#' Interpretation of parameters therefore depends on the outcome type:
#' \itemize{
#'   \item For continuous outcomes, parameters represent mean differences.
#'   \item For binary outcomes, parameters represent differences in predicted
#'   probabilities (risk differences).
#'   \item For count outcomes, parameters represent differences in expected counts.
#' }
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
#' fit <- qgcomp.glm.multi.boot(
#'   f = Y ~ X1 + X2 + X3 + X4 + W1 + W2 + W3 + W4 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3", "X4"),
#'   mix2 = c("W1", "W2", "W3", "W4"),
#'   interaction = TRUE,
#'   q = 4,
#'   B = 100,
#'   MCsize = nrow(dat)
#' )
#'
#' fit$coef_table
#'
#' @export

qgcomp.glm.multi.boot <- function(f,
                                  data,
                                  mix1,
                                  mix2,
                                  interaction = TRUE,
                                  family = gaussian(),
                                  q = 4,
                                  centering = "none",
                                  B = 200,
                                  id = NULL,
                                  MCsize = nrow(data)){

  data_q <- quantize_mixtures(data, mix1, mix2, q)

  coefs <- msm_fit(f, data_q, mix1, mix2, interaction, family, q, centering, id, MCsize)

  if (interaction){
    psi_hat <- matrix(NA, B, 4)
  } else {
    psi_hat <- matrix(NA, B, 3)
  }

  for (b in 1:B){

    if (is.null(id)) {

      idx <- sample(1:nrow(data_q), replace = T)
      data_b <- data_q[idx, ]

    } else {

      ids <- data_q[[id]]
      unique_ids <- unique(ids)

      samp_ids <- sample(unique_ids, length(unique_ids), replace = TRUE)

      split_data <- split(data_q, ids)

      data_b <- do.call(
        rbind,
        split_data[as.character(samp_ids)]
      )
    }

    psi_hat[b, ] <- msm_fit(f, data_b, mix1, mix2, interaction, family, q, centering, id, MCsize)
  }

  int <- var(psi_hat[,1])
  var1 <- var(psi_hat[,2])
  var2 <- var(psi_hat[,3])

  if (interaction){
    var12 <- var(psi_hat[,4])
  }

  cov12 <- cov(psi_hat[,2], psi_hat[,3])

  if (interaction) {
    std.err <- c(sqrt(int), sqrt(var1), sqrt(var2), sqrt(var12))
    names(std.err) <- c("intercept", "psi1", "psi2", "psi1:psi2")
  } else {
    std.err <- c(sqrt(int), sqrt(var1), sqrt(var2))
    names(std.err) <- c("intercept", "psi1", "psi2")
  }

  varcov <- matrix(c(var1, cov12, cov12, var2), nrow = 2, byrow = T)

  fit <- (list(coefs = coefs,
               std.err = std.err,
               varcov = varcov))

  make_coef_table <- function(coefs, se) {
    z_val <- coefs / se
    pval <- 2 * pnorm(abs(z_val), lower.tail = FALSE)

    cbind(
      Estimate = coefs,
      `Std. Error` = se,
      `z value` = z_val,
      `Pr(>|z|)` = pval
    )
  }

  fit$coef_table <- make_coef_table(coefs, std.err)

  fit
}

