#' Glance summaries for pooled qgcompmulti multiple-imputation fits
#'
#' Returns a one-row metadata summary for a fitted `qgcompmulti_mi` object.
#' The output is intentionally compact and inference-focused.
#'
#' @param x A fitted `qgcompmulti_mi` object.
#' @param ... Not used.
#'
#' @return A one-row data frame containing compact pooled-fit metadata for the
#'   `qgcompmulti_mi` object, including sample-size information, the number of
#'   imputations, input-type metadata, mixture coding metadata, Monte Carlo
#'   size, interaction status, retained-fit status, and clustering metadata.
#'
#' @seealso [broom::tidy()], [summary.qgcompmulti_mi()], [print.qgcompmulti_mi()],
#'   [qgcomp.glm.multi.mi()]
#' @exportS3Method broom::glance
glance.qgcompmulti_mi <- function(x, ...) {
  validate_qgcompmulti_mi(x)

  data.frame(
    n_input = x$data_info$n_input,
    n_used = x$data_info$n_used,
    m = x$mi_info$m,
    input_type = x$mi_info$input_type,
    family = x$analysis$family_name,
    link = x$analysis$link,
    quantized = x$data_info$quantized,
    q = if (is.null(x$mixtures$q)) NA_integer_ else as.integer(x$mixtures$q),
    centering = x$mixtures$centering,
    MCsize = x$analysis$MCsize,
    interaction = x$analysis$interaction,
    keep_fits = x$mi_info$keep_fits,
    has_clusters = x$data_info$has_clusters,
    cluster_var = if (isTRUE(x$data_info$has_clusters)) x$data_info$cluster_var else NA_character_,
    n_clusters = if (isTRUE(x$data_info$has_clusters)) x$data_info$n_clusters else NA_integer_,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
