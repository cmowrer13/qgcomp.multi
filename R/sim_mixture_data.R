# ------------------------------------------------------------------------------
# Internal simulation helpers
# ------------------------------------------------------------------------------
qgcompmulti_validate_sim_scalar <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("`%s` must be a single finite numeric value.", arg), call. = FALSE)
  }
}

qgcompmulti_validate_sim_probability <- function(x, arg, open = FALSE) {
  qgcompmulti_validate_sim_scalar(x, arg)

  if (isTRUE(open)) {
    if (x <= 0 || x >= 1) {
      stop(sprintf("`%s` must lie strictly between 0 and 1.", arg), call. = FALSE)
    }
  } else if (x < 0 || x > 1) {
    stop(sprintf("`%s` must lie between 0 and 1.", arg), call. = FALSE)
  }
}

qgcompmulti_validate_sim_count <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || x < 1L) {
    stop(sprintf("`%s` must be a single positive integer.", arg), call. = FALSE)
  }
}

qgcompmulti_validate_sim_correlation <- function(x, arg) {
  qgcompmulti_validate_sim_scalar(x, arg)
  if (x < -1 || x > 1) {
    stop(sprintf("`%s` must lie between -1 and 1.", arg), call. = FALSE)
  }
}

qgcompmulti_validate_sim_correlation_matrix <- function(Sigma) {
  chol_ok <- tryCatch(
    {
      chol(Sigma)
      TRUE
    },
    error = function(...) FALSE
  )

  if (!chol_ok) {
    stop(
      paste(
        "The requested within- and between-mixture correlations do not produce",
        "a positive-definite exposure correlation matrix."
      ),
      call. = FALSE
    )
  }
}

qgcompmulti_quantize_exposures <- function(mat, q) {
  out <- apply(mat, 2, function(v) dplyr::ntile(v, q) - 1L)

  if (is.null(dim(out))) {
    out <- matrix(out, ncol = 1L)
  }

  storage.mode(out) <- "numeric"
  out
}

#' Generate simulated dataset with two exposure mixtures
#'
#' Generates a dataset for simulation studies involving two mixtures of
#' correlated exposures and their interaction. Exposures are generated from
#' a multivariate normal distribution with user-specified within-mixture and
#' between-mixture correlations, then discretized into quantiles. The outcome
#' is generated from family- and estimand-scale-specific models where the
#' parameters correspond to simultaneous one-quantile increases in each
#' mixture.
#'
#' The exposure summaries are:
#'
#' \deqn{S_A = \textrm{rowMeans}(X_q), \qquad S_B = \textrm{rowMeans}(W_q)}
#'
#' where `X_q` and `W_q` are the quantized mixture components. The outcome is
#' then generated from one of the following models:
#'
#' \itemize{
#'   \item Gaussian mean-difference generation:
#'   \deqn{Y = \mu_0 + \psi_1 S_A + \psi_2 S_B + \psi_{12} S_A S_B + \beta_C C + \epsilon}
#'   \item Binomial odds-ratio generation:
#'   \deqn{\textrm{logit}\{P(Y = 1)\} = \textrm{logit}(p_0) + \log(\psi_1) S_A + \log(\psi_2) S_B + \log(\psi_{12}) S_A S_B + \beta_C C}
#'   \item Binomial risk-difference generation:
#'   \deqn{P(Y = 1) = p_0 + \psi_1 S_A + \psi_2 S_B + \psi_{12} S_A S_B + \beta_C C}
#'   \item Poisson rate-ratio generation:
#'   \deqn{\log\{E(Y)\} = \log(\lambda_0) + \log(\psi_1) S_A + \log(\psi_2) S_B + \log(\psi_{12}) S_A S_B + \beta_C C}
#'   \item Poisson mean-difference generation:
#'   \deqn{E(Y) = \lambda_0 + \psi_1 S_A + \psi_2 S_B + \psi_{12} S_A S_B + \beta_C C}
#' }
#'
#' Ratio-valued inputs are accepted on their natural scale. For example,
#' `psi1 = 1.5` with `family = poisson()` and `estimand_scale = "rate_ratio"`
#' means that a one-quantile increase in mixture A corresponds to a rate ratio
#' of `1.5` before any interaction modification by mixture B.
#'
#' @param n Number of observations.
#' @param pA Number of components in mixture A.
#' @param pB Number of components in mixture B.
#' @param rho_within_A Correlation between components within mixture A.
#' @param rho_within_B Correlation between components within mixture B.
#' @param rho_between Correlation between components across mixtures.
#' @param psi1 Effect of a one-quantile increase in mixture A.
#' @param psi2 Effect of a one-quantile increase in mixture B.
#' @param psi12 Interaction effect representing the change in the effect of
#' mixture A when mixture B increases by one quantile.
#' @param sigma_eps Standard deviation of the outcome error term for Gaussian
#' generation.
#' @param q Number of quantiles used to discretize exposures.
#' @param return_quantized Logical; if TRUE returns quantized exposures,
#' otherwise returns the original continuous exposures.
#' @param seed Optional random seed for reproducibility.
#' @param family Outcome family. Supported families are `gaussian()`,
#' `binomial()`, and `poisson()`.
#' @param estimand_scale Optional outcome-generating estimand scale. Supported
#' values mirror the fitting API: `"mean_difference"`, `"risk_difference"`,
#' `"odds_ratio"`, and `"rate_ratio"`. If omitted, defaults follow the `0.5.0`
#' fitting conventions.
#' @param baseline_mean Baseline mean used for Gaussian generation.
#' @param baseline_risk Baseline risk used for binomial generation.
#' @param baseline_rate Baseline rate used for Poisson generation.
#' @param beta_C Effect of covariate `C` in the outcome-generating model.
#'
#' @return A data frame containing the simulated outcome `Y`, mixture
#' exposures (`X1...XpA` and `W1...WpB`), and covariate `C`.
#'
#' @details
#' Continuous exposures are generated from a multivariate normal distribution
#' with a block correlation structure. Exposures are then discretized into
#' quantile categories ranging from `0` to `q - 1`. The outcome is always
#' generated from the quantized mixture summaries so that the requested
#' parameters `psi1`, `psi2`, and `psi12` correspond directly to the intended
#' one-quantile mixture effects on the chosen generating scale.
#'
#' Additive binary-risk and additive Poisson-mean generation are supported, but
#' the requested baseline and effect values must imply valid probabilities or
#' strictly positive means for all simulated observations.
#'
#' @examples
#' dat_gaussian <- sim_mixture_data(
#'   n = 1000,
#'   pA = 4,
#'   pB = 4,
#'   rho_within_A = 0.5,
#'   rho_within_B = 0.5,
#'   rho_between = 0.3,
#'   psi1 = 0.5,
#'   psi2 = 0.3,
#'   psi12 = 0.2,
#'   return_quantized = TRUE,
#'   seed = 123
#' )
#'
#' dat_binomial <- sim_mixture_data(
#'   n = 1000,
#'   pA = 3,
#'   pB = 3,
#'   rho_within_A = 0.4,
#'   rho_within_B = 0.4,
#'   rho_between = 0.2,
#'   psi1 = 1.5,
#'   psi2 = 1.2,
#'   psi12 = 1.1,
#'   family = binomial(),
#'   baseline_risk = 0.2,
#'   beta_C = 0.25,
#'   seed = 123
#' )
#'
#' dat_poisson <- sim_mixture_data(
#'   n = 1000,
#'   pA = 3,
#'   pB = 3,
#'   rho_within_A = 0.4,
#'   rho_within_B = 0.4,
#'   rho_between = 0.2,
#'   psi1 = 1.5,
#'   psi2 = 1.2,
#'   psi12 = 1.1,
#'   family = poisson(),
#'   baseline_rate = 1,
#'   beta_C = 0.15,
#'   seed = 123
#' )
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom dplyr ntile
sim_mixture_data <- function(
    n,
    pA,
    pB,
    rho_within_A,
    rho_within_B,
    rho_between,
    psi1,
    psi2,
    psi12,
    sigma_eps = 1,
    q = 4,
    return_quantized = FALSE,
    seed = NULL,
    family = stats::gaussian(),
    estimand_scale = NULL,
    baseline_mean = 0,
    baseline_risk = 0.2,
    baseline_rate = 1,
    beta_C = 1
) {

  qgcompmulti_validate_sim_count(n, "n")
  qgcompmulti_validate_sim_count(pA, "pA")
  qgcompmulti_validate_sim_count(pB, "pB")
  qgcompmulti_validate_sim_correlation(rho_within_A, "rho_within_A")
  qgcompmulti_validate_sim_correlation(rho_within_B, "rho_within_B")
  qgcompmulti_validate_sim_correlation(rho_between, "rho_between")
  qgcompmulti_validate_sim_scalar(psi1, "psi1")
  qgcompmulti_validate_sim_scalar(psi2, "psi2")
  qgcompmulti_validate_sim_scalar(psi12, "psi12")
  qgcompmulti_validate_sim_scalar(sigma_eps, "sigma_eps")
  qgcompmulti_validate_sim_scalar(baseline_mean, "baseline_mean")
  qgcompmulti_validate_sim_probability(baseline_risk, "baseline_risk", open = TRUE)
  qgcompmulti_validate_sim_scalar(baseline_rate, "baseline_rate")
  qgcompmulti_validate_sim_scalar(beta_C, "beta_C")

  if (sigma_eps <= 0) {
    stop("`sigma_eps` must be strictly positive.", call. = FALSE)
  }

  if (baseline_rate <= 0) {
    stop("`baseline_rate` must be strictly positive.", call. = FALSE)
  }

  validate_q_argument(q)

  if (!is.logical(return_quantized) || length(return_quantized) != 1L ||
      is.na(return_quantized)) {
    stop("`return_quantized` must be either `TRUE` or `FALSE`.", call. = FALSE)
  }

  if (!is.null(seed) && !is_scalar_whole_number(seed)) {
    stop("`seed` must be `NULL` or a single integer.", call. = FALSE)
  }

  if (!inherits(family, "family")) {
    stop(
      "`family` must be a valid GLM family object, such as gaussian(), binomial(), or poisson().",
      call. = FALSE
    )
  }

  family_name <- family$family
  link <- family$link

  if (!family_name %in% c("gaussian", "binomial", "poisson")) {
    stop(
      "`sim_mixture_data()` currently supports only gaussian(), binomial(), and poisson() families.",
      call. = FALSE
    )
  }

  qgcompmulti_validate_estimand_scale(
    family_name = family_name,
    link = link,
    estimand_scale = estimand_scale,
    allow_null = TRUE
  )

  estimand_scale <- if (is.null(estimand_scale)) {
    qgcompmulti_default_estimand_scale(
      family_name = family_name,
      link = link,
      mode = "planned"
    )
  } else {
    estimand_scale
  }

  if (estimand_scale %in% c("odds_ratio", "rate_ratio")) {
    if (any(c(psi1, psi2, psi12) <= 0)) {
      stop(
        paste(
          "Ratio-scale simulation requires `psi1`, `psi2`, and `psi12`",
          "to be strictly positive."
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  p_total <- pA + pB

  Sigma <- matrix(rho_between, nrow = p_total, ncol = p_total)
  Sigma[1:pA, 1:pA] <- rho_within_A
  Sigma[(pA + 1):p_total, (pA + 1):p_total] <- rho_within_B
  diag(Sigma) <- 1

  qgcompmulti_validate_sim_correlation_matrix(Sigma)

  Z <- MASS::mvrnorm(n, mu = rep(0, p_total), Sigma = Sigma)

  X <- Z[, 1:pA, drop = FALSE]
  W <- Z[, (pA + 1):p_total, drop = FALSE]

  colnames(X) <- paste0("X", seq_len(pA))
  colnames(W) <- paste0("W", seq_len(pB))

  Xq <- qgcompmulti_quantize_exposures(X, q)
  Wq <- qgcompmulti_quantize_exposures(W, q)

  colnames(Xq) <- colnames(X)
  colnames(Wq) <- colnames(W)

  SA <- rowMeans(Xq)
  SB <- rowMeans(Wq)
  C <- stats::rnorm(n)

  if (identical(family_name, "gaussian")) {
    mean_y <- baseline_mean +
      psi1 * SA +
      psi2 * SB +
      psi12 * SA * SB +
      beta_C * C

    Y <- as.numeric(mean_y + stats::rnorm(n, mean = 0, sd = sigma_eps))
  } else if (identical(family_name, "binomial")) {
    if (identical(estimand_scale, "odds_ratio")) {
      lp <- stats::qlogis(baseline_risk) +
        log(psi1) * SA +
        log(psi2) * SB +
        log(psi12) * SA * SB +
        beta_C * C

      prob <- stats::plogis(lp)
    } else {
      prob <- baseline_risk +
        psi1 * SA +
        psi2 * SB +
        psi12 * SA * SB +
        beta_C * C

      if (any(prob < 0 | prob > 1)) {
        stop(
          paste(
            "The requested binomial additive simulation produced probabilities",
            "outside [0, 1]. Adjust `baseline_risk`, `beta_C`, or the `psi` values."
          ),
          call. = FALSE
        )
      }
    }

    Y <- stats::rbinom(n = n, size = 1L, prob = prob)
  } else {
    if (identical(estimand_scale, "rate_ratio")) {
      eta <- log(baseline_rate) +
        log(psi1) * SA +
        log(psi2) * SB +
        log(psi12) * SA * SB +
        beta_C * C

      mu <- exp(eta)
    } else {
      mu <- baseline_rate +
        psi1 * SA +
        psi2 * SB +
        psi12 * SA * SB +
        beta_C * C

      if (any(mu <= 0)) {
        stop(
          paste(
            "The requested Poisson additive simulation produced nonpositive means.",
            "Adjust `baseline_rate`, `beta_C`, or the `psi` values."
          ),
          call. = FALSE
        )
      }
    }

    Y <- stats::rpois(n = n, lambda = mu)
  }

  data_cont <- data.frame(Y = Y, X, W, C = C)
  data_quant <- data.frame(Y = Y, Xq, Wq, C = C)

  if (return_quantized) {
    data_quant
  } else {
    data_cont
  }
}
