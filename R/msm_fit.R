#' Fit marginal structural model for two mixtures
#'
#' @param f Model formula
#' @param data Data frame
#' @param mix1 Names of mixture 1 variables
#' @param mix2 Names of mixture 2 variables
#' @param interaction Include interaction term
#' @param q Number of quantiles
#' @param id Cluster id
#' @param MCsize Monte Carlo sample size
#'
#' @return Estimated MSM coefficients
#' @export

msm_fit <- function(f,
                    data,
                    mix1,
                    mix2,
                    interaction = TRUE,
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

  fit <- glm(newf, data = data)

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
