#' Glance summaries for qgcompmulti fits
#'
#' Returns a one-row metadata summary for a fitted `qgcompmulti` object. The
#' output stays compact and fit-oriented rather than repeating coefficient-level
#' results or diagnostic summaries.
#'
#' @param x A fitted `qgcompmulti` object.
#' @param ... Not used.
#'
#' @return A one-row data frame containing compact fit metadata for the
#'   `qgcompmulti` object, including sample-size information, outcome-family
#'   metadata, mixture coding metadata, bootstrap replication counts, Monte
#'   Carlo size, interaction status, and clustering metadata.
#'
#' @details
#' For original-scale fits (`q = NULL`), `q` is returned as `NA_integer_` and
#' the `centering` column records how the MSM intervention variables were coded.
#' For quantized fits, `centering` still reports the stored setting from the fit
#' object, but it is not inferentially active.
#'
#' @seealso [broom::tidy()], [summary.qgcompmulti()], [print.qgcompmulti()],
#'   [qgcomp.glm.multi()]
#' @exportS3Method broom::glance
glance.qgcompmulti <- function(x, ...) {
  validate_qgcompmulti(x)
  data.frame(
    n_input = x$data_info$n_input,
    n_used = nobs(x),
    family = x$analysis$family_name,
    link = x$analysis$link,
    quantized = x$data_info$quantized,
    q = if (is.null(x$mixtures$q)) NA_integer_ else as.integer(x$mixtures$q),
    centering = x$mixtures$centering,
    B_requested = x$bootstrap$B_requested,
    B_success = x$bootstrap$B_success,
    B_failed = x$bootstrap$B_failed,
    MCsize = x$analysis$MCsize,
    interaction = x$analysis$interaction,
    has_clusters = x$data_info$has_clusters,
    cluster_var = if (isTRUE(x$data_info$has_clusters)) x$data_info$cluster_var else NA_character_,
    n_clusters = if (isTRUE(x$data_info$has_clusters)) x$data_info$n_clusters else NA_integer_,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
