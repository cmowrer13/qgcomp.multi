#' Sensitivity to Monte Carlo size in qgcompmulti fits
#'
#' Re-fits [qgcomp.glm.multi()] across multiple `MCsize` values while
#' preserving the rest of the analysis specification.
#'
#' @param f,data,mix1,mix2,interaction,family,q,centering,B,id,seed Arguments
#' passed through to [qgcomp.glm.multi()].
#' @param MCsize_values Integer vector of Monte Carlo sizes to compare.
#' @param keep_fits Logical; if `TRUE`, retain the full fitted objects.
#'
#' @return An object of class `"qgcompmulti_mcsize_sensitivity"`.
#' @export
mcsize_sensitivity <- function(f,
                               data,
                               mix1,
                               mix2,
                               MCsize_values,
                               interaction = TRUE,
                               family = gaussian(),
                               q = 4,
                               centering = "none",
                               B = 200,
                               id = NULL,
                               seed = NULL,
                               keep_fits = TRUE) {
  MCsize_values <- qgcompmulti_mcsize_values(MCsize_values, nrow(data))
  fits <- lapply(MCsize_values, function(mcsize) {
    qgcomp.glm.multi(
      f = f,
      data = data,
      mix1 = mix1,
      mix2 = mix2,
      interaction = interaction,
      family = family,
      q = q,
      centering = centering,
      B = B,
      id = id,
      MCsize = mcsize,
      seed = seed
    )
  })
  results_table <- qgcompmulti_build_sensitivity_table(
    values = as.list(MCsize_values),
    parameter_name = "MCsize",
    fits = fits
  )
  structure(
    list(
      workflow_type = "mcsize_sensitivity",
      call_template = match.call(),
      settings_fixed = list(
        interaction = interaction,
        family = family,
        q = q,
        centering = centering,
        B = B,
        id = id,
        seed = seed
      ),
      MCsize_values = MCsize_values,
      results_table = results_table,
      fits = if (isTRUE(keep_fits)) fits else NULL,
      notes = "Only MCsize varies across refits; all other analysis settings are held fixed."
    ),
    class = "qgcompmulti_mcsize_sensitivity"
  )
}
#' Sensitivity to quantization choice in qgcompmulti fits
#'
#' Re-fits [qgcomp.glm.multi()] across multiple integer `q` values while
#' preserving the rest of the analysis specification.
#'
#' @param f,data,mix1,mix2,interaction,family,centering,B,id,MCsize,seed
#' Arguments passed through to [qgcomp.glm.multi()].
#' @param q_values Integer vector of quantization choices to compare.
#' @param keep_fits Logical; if `TRUE`, retain the full fitted objects.
#'
#' @return An object of class `"qgcompmulti_q_sensitivity"`.
#' @export
q_sensitivity <- function(f,
                          data,
                          mix1,
                          mix2,
                          q_values,
                          interaction = TRUE,
                          family = gaussian(),
                          centering = "none",
                          B = 200,
                          id = NULL,
                          MCsize = nrow(data),
                          seed = NULL,
                          keep_fits = TRUE) {
  q_values <- qgcompmulti_q_values(q_values)
  fits <- lapply(q_values, function(q) {
    qgcomp.glm.multi(
      f = f,
      data = data,
      mix1 = mix1,
      mix2 = mix2,
      interaction = interaction,
      family = family,
      q = q,
      centering = centering,
      B = B,
      id = id,
      MCsize = MCsize,
      seed = seed
    )
  })
  results_table <- qgcompmulti_build_sensitivity_table(
    values = as.list(q_values),
    parameter_name = "q",
    fits = fits
  )
  structure(
    list(
      workflow_type = "q_sensitivity",
      call_template = match.call(),
      settings_fixed = list(
        interaction = interaction,
        family = family,
        centering = centering,
        B = B,
        id = id,
        MCsize = MCsize,
        seed = seed
      ),
      q_values = q_values,
      results_table = results_table,
      fits = if (isTRUE(keep_fits)) fits else NULL,
      comparability_note = qgcompmulti_q_comparability_note(),
      notes = "Only q varies across refits; all other analysis settings are held fixed."
    ),
    class = "qgcompmulti_q_sensitivity"
  )
}
#' @export
print.qgcompmulti_mcsize_sensitivity <- function(x, ...) {
  cat("qgcompmulti MCsize sensitivity\n\n")
  cat(x$notes, "\n\n", sep = "")
  qgcompmulti_print_sensitivity_table(x$results_table)
  invisible(x)
}
#' @export
print.qgcompmulti_q_sensitivity <- function(x, ...) {
  cat("qgcompmulti q sensitivity\n\n")
  cat(x$notes, "\n\n", sep = "")
  cat("Comparability note:\n")
  cat("  ", x$comparability_note, "\n\n", sep = "")
  qgcompmulti_print_sensitivity_table(x$results_table)
  invisible(x)
}
