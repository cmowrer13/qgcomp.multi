#' Generate simulated dataset with two exposure mixtures
#'
#' Generates a dataset for simulation studies involving two mixtures of
#' correlated exposures and their interaction. Exposures are generated from
#' a multivariate normal distribution with user-specified within-mixture and
#' between-mixture correlations, then discretized into quantiles. The outcome
#' is generated from a model where the effects correspond to simultaneous
#' one-quantile increases in each mixture.
#'
#' The data-generating model is:
#'
#' Y = psi1 * S_A + psi2 * S_B + psi12 * S_A * S_B + C + error
#'
#' where S_A and S_B are the average quantile levels of the exposures in
#' mixture A and mixture B, respectively.
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
#' @param sigma_eps Standard deviation of the outcome error term.
#' @param q Number of quantiles used to discretize exposures.
#' @param return_quantized Logical; if TRUE returns quantized exposures,
#' otherwise returns the original continuous exposures.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A data frame containing the simulated outcome `Y`, mixture
#' exposures (`X1...XpA` and `W1...WpB`), and covariate `C`.
#'
#' @details
#' Continuous exposures are generated from a multivariate normal distribution
#' with a block correlation structure. Exposures are then discretized into
#' quantile categories ranging from 0 to q-1. The outcome is generated using
#' the quantized exposures so that the parameters psi1, psi2, and psi12
#' correspond directly to the estimands in two-mixture quantile
#' g-computation models.
#'
#' @examples
#' dat <- sim_mixture_data(
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
  seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  p_total <- pA + pB

  Sigma <- matrix(rho_between, nrow = p_total, ncol = p_total)

  Sigma[1:pA, 1:pA] <- rho_within_A
  Sigma[(pA + 1):(p_total), (pA + 1):(p_total)] <- rho_within_B

  diag(Sigma) <- 1

  Z <- MASS::mvrnorm(n, mu = rep(0, p_total), Sigma = Sigma)

  X <- Z[, 1:pA, drop = FALSE]
  W <- Z[, (pA + 1):(p_total), drop = FALSE]

  colnames(X) <- paste0("X", 1:pA)
  colnames(W) <- paste0("W", 1:pB)

  quantize <- function(mat, q) {
    apply(mat, 2, function(v) {
      dplyr::ntile(v, q) - 1
    })
  }

  Xq <- quantize(X, q)
  Wq <- quantize(W, q)

  colnames(Xq) <- colnames(X)
  colnames(Wq) <- colnames(W)

  SA <- rowMeans(Xq)
  SB <- rowMeans(Wq)

  C <- rnorm(n)

  Y <- as.numeric(psi1 * SA +
                    psi2 * SB +
                    psi12 * SA * SB +
                    C +
                    rnorm(n, 0, sigma_eps))

  data_cont <- data.frame(Y, X, W, C)
  data_quant <- data.frame(Y, Xq, Wq, C)

  if (return_quantized) {
    return(data_quant)
  } else {
    return(data_cont)
  }
}
