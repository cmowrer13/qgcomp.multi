# ------------------------------------------------------------------------------
# Plot and presentation helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_deparse_one_line <- function(x) {
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}
#' @keywords internal
#' @noRd
qgcompmulti_family_label <- function(analysis) {
  paste0(analysis$family_name, " (", analysis$link, ")")
}
#' @keywords internal
#' @noRd
qgcompmulti_exposure_mode_label <- function(mixtures) {
  if (is.null(mixtures$q)) {
    return(
      paste0(
        "Original-scale exposures (centering = \"",
        mixtures$centering,
        "\")"
      )
    )
  }
  paste0("Quantized exposures (q = ", mixtures$q, ")")
}
#' @keywords internal
#' @noRd
qgcompmulti_interaction_label <- function(analysis) {
  if (isTRUE(analysis$interaction)) {
    return("MSM interaction: included")
  }
  "MSM interaction: not included"
}
#' @keywords internal
#' @noRd
qgcompmulti_mixture_definition <- function(label, vars) {
  paste0(label, ": ", paste(vars, collapse = ", "))
}
#' @keywords internal
#' @noRd
qgcompmulti_labeled_coef_table <- function(results, labels) {
  coef_table <- results$coef_table
  if (is.null(coef_table)) {
    return(NULL)
  }
  if (!is.data.frame(coef_table)) {
    stop("`results$coef_table` must be a data frame.", call. = FALSE)
  }
  coef_names <- rownames(coef_table)
  coef_labels <- labels$coefficient_labels[coef_names]
  if (anyNA(coef_labels)) {
    stop(
      "Missing coefficient labels for one or more MSM coefficients.",
      call. = FALSE
    )
  }
  out <- coef_table
  rownames(out) <- unname(coef_labels)
  out
}
#' @keywords internal
#' @noRd
qgcompmulti_outcome_model_info <- function(outcome_fit) {
  if (is.null(outcome_fit)) {
    return(
      list(
        model_class = NULL,
        n_parameters = NULL,
        aic = NULL,
        null_deviance = NULL,
        residual_deviance = NULL
      )
    )
  }
  list(
    model_class = class(outcome_fit)[1],
    n_parameters = length(stats::coef(outcome_fit)),
    aic = stats::AIC(outcome_fit),
    null_deviance = outcome_fit$null.deviance,
    residual_deviance = outcome_fit$deviance
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_format_metric <- function(x, digits = 3) {
  if (is.null(x) || is.na(x) || !is.finite(x)) {
    return("NA")
  }
  formatC(x, digits = digits, format = "f")
}
#' @keywords internal
#' @noRd
qgcompmulti_print_call <- function(call) {
  cat("Call:\n")
  print(call)
}
#' @keywords internal
#' @noRd
qgcompmulti_print_mixtures <- function(mixtures) {
  cat("Mixtures:\n")
  cat(
    "  ",
    qgcompmulti_mixture_definition(mixtures$labels[["mix1"]], mixtures$mix1),
    "\n",
    sep = ""
  )
  cat(
    "  ",
    qgcompmulti_mixture_definition(mixtures$labels[["mix2"]], mixtures$mix2),
    "\n",
    sep = ""
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_print_msm_table <- function(msm_table) {
  if (is.null(msm_table) || nrow(msm_table) == 0L) {
    cat("MSM coefficients: none available\n")
    return(invisible(NULL))
  }
  cat("MSM coefficients:\n")
  stats::printCoefmat(
    as.matrix(msm_table),
    P.values = TRUE,
    has.Pvalue = TRUE,
    signif.stars = getOption("show.signif.stars")
  )
  invisible(msm_table)
}
