#' Fit quantile g-computation model for two mixtures
#'
#' @param f Model formula
#' @param data Data frame
#' @param mix1 Names of mixture 1 variables
#' @param mix2 Names of mixture 2 variables
#' @param interaction Include interaction term
#' @param q Number of quantiles
#' @param B Number of bootstrap iterations
#' @param id Cluster id
#' @param MCsize Monte Carlo sample size
#'
#' @return Fitted qgcomp model
#' @export

qgcomp.glm.multi.boot <- function(f,
                                  data,
                                  mix1,
                                  mix2,
                                  interaction = TRUE,
                                  q = 4,
                                  B = 200,
                                  id = NULL,
                                  MCsize = nrow(data)){

  data_q <- quantize_mixtures(data, mix1, mix2, q)

  coefs <- msm_fit(f, data_q, mix1, mix2, interaction, q, id, MCsize)

  if (interaction){
    psi_hat <- matrix(NA, B, 3)
  } else {
    psi_hat <- matrix(NA, B, 2)
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

    psi_hat[b, ] <- msm_fit(f, data_b, mix1, mix2, interaction, q, id, MCsize)
  }

  var1 <- var(psi_hat[,1])
  var2 <- var(psi_hat[,2])

  if (interaction){
    var12 <- var(psi_hat[,3])
  }

  cov12 <- cov(psi_hat[,1], psi_hat[,2])

  if (interaction) {
    std.err <- c(sqrt(var1), sqrt(var2), sqrt(var12))
    names(std.err) <- c("psi1", "psi2", "psi1:psi2")
  } else {
    std.err <- c(sqrt(var1), sqrt(var2))
    names(std.err) <- c("psi1", "psi2")
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
