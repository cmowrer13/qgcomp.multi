#' Quantize exposure variables for two-mixture quantile g-computation
#'
#' Converts continuous exposure variables belonging to two mixtures into
#' discrete quantile categories. Each exposure is partitioned into `q`
#' approximately equal-sized groups using empirical quantiles, and replaced
#' with integer values ranging from 0 to `q - 1`. This transformation places
#' mixture components on a common scale, enabling interpretation of model
#' parameters as effects of simultaneous one-quantile increases in each
#' component.
#'
#' @param data A data frame containing the exposure variables and any
#' additional variables (e.g., outcome or covariates).
#' @param mix1 A character vector giving the names of the variables in the
#' first exposure mixture.
#' @param mix2 A character vector giving the names of the variables in the
#' second exposure mixture.
#' @param q Integer giving the number of quantiles used to discretize each
#' exposure variable, or `NULL` to skip quantization and return the mixture
#' variables unchanged.
#'
#' @return A data frame of the same dimensions as `data`, When `q` is numeric,
#' the variables listed in `mix1` and `mix2` have been replaced by their
#' quantized versions. Quantized variables take integer values in
#' `{0, 1, ..., q - 1}`. When `q = NULL`, the returned data are unchanged.
#'
#' @details
#' Quantization is a key step in quantile g-computation, as it ensures that
#' exposure variables measured on different scales are made comparable. By
#' transforming each exposure into quantile categories, a one-unit increase
#' corresponds to a shift of one quantile, allowing mixture effects to be
#' interpreted as the joint effect of increasing all components by one
#' quantile.
#'
#' Quantiles assign observations
#' to groups of approximately equal size. The resulting categories are then
#' shifted to start at 0 (rather than 1) to align with the intervention levels
#' used in the marginal structural model.
#'
#' @examples
#' dat <- sim_mixture_data(
#'   n = 100,
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
#' dat_q <- quantize_mixtures(
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   q = 4
#' )
#'
#' head(dat_q)
#'
#' @export

quantize_mixtures <- function(data, mix1, mix2, q = 4){

  if (is.null(q)) {
    return(data)
  }

  vars <- c(mix1, mix2)

  for (v in vars){

    bins <- cut(
      data[[v]],
      quantile(data[[v]],
               probs = seq(0, 1, length.out = q + 1),
               na.rm = TRUE),
      include.lowest = TRUE,
      labels = FALSE
    )

    data[[v]] <- as.numeric(bins - 1)
  }

  data
}
