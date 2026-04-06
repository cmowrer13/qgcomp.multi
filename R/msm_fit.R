#' Fit the marginal structural model for two-mixture quantile g-computation
#'
#' Fits the outcome regression, computes predicted potential outcomes under
#' joint interventions on two exposure mixtures, and estimates the marginal
#' structural model (MSM) coefficients that summarize the resulting
#' two-dimensional dose-response surface. This function is primarily an
#' internal helper used by [qgcomp.glm.multi.boot()] and will not typically
#' need to be called directly by users.
#'
#' Given a fitted outcome model, the function evaluates predicted outcomes over
#' a grid of uniform mixture quantiles and regresses these predictions on the
#' mixture intervention levels to obtain MSM parameters. When
#' `interaction = TRUE`, the fitted MSM has the form
#'
#' \deqn{
#' E[Y^{x(q1), w(q2)}] = \psi_1 q_1 + \psi_2 q_2 + \psi_{12} q_1 q_2
#' }
#'
#' where `psi1` and `psi2` represent the main effects of one-quantile
#' increases in mixtures 1 and 2, respectively, and `psi12` represents their
#' interaction on the MSM scale.
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
#' exposure variables and define the intervention grid.
#' @param id Optional character string giving the name of a cluster identifier
#' variable. If supplied, Monte Carlo subsampling is performed at the cluster
#' level rather than the observation level.
#' @param MCsize Optional integer controlling the Monte Carlo sample size used
#' to approximate the marginalization step in g-computation. If equal to
#' `nrow(data)`, all observations are used. Smaller values compute predicted
#' outcomes over a random subsample drawn from the empirical distribution.
#'
#' @return A numeric vector of estimated marginal structural model coefficients.
#' When `interaction = TRUE`, the returned coefficients are `(Intercept)`,
#' `psi1`, `psi2`, and `psi1:psi2`. When `interaction = FALSE`, only the
#' intercept and main-effect coefficients are returned.
#'
#' @details
#' This function carries out the core g-computation step for the two-mixture
#' extension. It first fits a generalized linear model for the outcome using
#' `glm()`. It then constructs a grid of intervention levels over the two
#' mixtures, replaces the observed mixture values with the intervention values,
#' and computes predicted outcomes for each observation (or for a Monte Carlo
#' subsample if `MCsize < nrow(data)`). These predicted outcomes are stacked
#' into a pseudo-dataset and used to fit a marginal structural model that
#' summarizes the dose-response surface. Note that predicted  potential
#' outcomes are computed on the response scale and the marginal
#' structural model is fit using an identity link, regardless of the outcome
#' model.
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
#' msm_fit(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   interaction = TRUE,
#'   q = 4,
#'   MCsize = nrow(dat)
#' )
#'
#' @export

msm_fit <- function(f,
                    data,
                    mix1,
                    mix2,
                    interaction = TRUE,
                    family = gaussian(),
                    q = 4,
                    id = NULL,
                    MCsize = nrow(data)){

  response <- deparse(f[[2]])

  rhs_terms <- attr(terms(f), "term.labels")

  rhs <- paste(rhs_terms, collapse = " + ")

  if (interaction) {
    m1 <- paste(mix1, collapse = " + ")
    m2 <- paste(mix2, collapse = " + ")

    interaction_term <- paste0("I((", m1, ") * (", m2, "))")
    rhs <- paste(rhs, interaction_term, sep = " + ")
  }

  newf <- as.formula(paste(response, "~", rhs))

  fit <- glm(newf, data = data, family = family)

  intgrid <- expand.grid(
    psi1 = 0:(q-1),
    psi2 = 0:(q-1)
  )

  if (MCsize == nrow(data)) {

    nd <- data

  } else {

    if (is.null(id)) {

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
  }

  nobs <- nrow(nd)

  pred_fun <- function(psi1, psi2, nd) {
    nd2 <- nd
    nd2[, mix1] <- psi1
    nd2[, mix2] <- psi2

    predict(fit, newdata = nd2, type = "response")
  }

  predmat <- Map(
    pred_fun,
    intgrid$psi1,
    intgrid$psi2,
    MoreArgs = list(nd = nd)
  )

  msmdat <- data.frame(
    Ya   = unlist(predmat),
    psi1 = rep(intgrid$psi1, each = nobs),
    psi2 = rep(intgrid$psi2, each = nobs)
  )

  if (interaction) {
    msm_form <- Ya ~ psi1 * psi2
  }  else {
    msm_form <- Ya ~ psi1 + psi2
  }

  msmfit <- glm(
    msm_form,
    data = msmdat
  )

  return(msmfit$coef)

}
