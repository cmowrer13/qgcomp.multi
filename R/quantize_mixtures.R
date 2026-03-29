#' Quantize mixture variables
#'
#' Converts dataset with continuous mixture components into quantized version.
#'
#' @param data Data frame
#' @param mix1 Names of mixture 1 variables
#' @param mix2 Names of mixture 2 variables
#' @param q Number of quantiles
#'
#' @return Data frame with quantized variables
#' @export

quantize_mixtures <- function(data, mix1, mix2, q = 4){

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
