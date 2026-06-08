#' Sensitivity to Monte Carlo size in qgcompmulti fits
#'
#' Re-fits [qgcomp.glm.multi()] across multiple `MCsize` values while
#' preserving the rest of the analysis specification.
#'
#' This helper is intended to answer a practical computational question: are
#' the fitted results reasonably stable as the Monte Carlo approximation size
#' changes?
#'
#' @param f,data,mix1,mix2,interaction,family,q,centering,B,id,seed Arguments
#' passed through to [qgcomp.glm.multi()].
#' @param MCsize_values Integer vector of Monte Carlo sizes to compare.
#' @param keep_fits Logical; if `TRUE`, retain the full fitted objects.
#'
#' @return An object of class `"qgcompmulti_mcsize_sensitivity"`.
#'
#' @details
#' `MCsize` sensitivity is implemented as a repeated-fit workflow rather than a
#' diagnostic of one existing fit. The helper keeps the model formula, outcome
#' family, mixture definitions, interaction setting, quantization choice,
#' bootstrap count, clustering, and seed fixed while varying only `MCsize`.
#'
#' The resulting object is designed for stability assessment rather than
#' automatic tuning. Users should look for whether the fitted coefficients,
#' adequacy summaries, and bootstrap behavior are reasonably consistent across
#' the requested `MCsize` values.
#'
#' @examples
#' \dontrun{
#' dat <- sim_mixture_data(
#'   n = 400,
#'   pA = 3,
#'   pB = 3,
#'   rho_within_A = 0.3,
#'   rho_within_B = 0.3,
#'   rho_between = 0.2,
#'   psi1 = 0.5,
#'   psi2 = 0.3,
#'   psi12 = 0.2,
#'   seed = 123
#' )
#'
#' mcsize_sensitivity(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   MCsize_values = c(100, 200, 400),
#'   q = 4,
#'   B = 100,
#'   seed = 13
#' )
#' }
#'
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
#' Sensitivity to bootstrap iteration count in qgcompmulti fits
#'
#' Re-fits [qgcomp.glm.multi()] across multiple `B` values while preserving
#' the rest of the analysis specification.
#'
#' This helper is intended to answer a practical inferential question: are the
#' reported bootstrap-based standard errors and related summaries reasonably
#' stable as the number of bootstrap replications changes?
#'
#' @param f,data,mix1,mix2,interaction,family,q,centering,id,MCsize,seed
#' Arguments passed through to [qgcomp.glm.multi()].
#' @param B_values Integer vector of bootstrap iteration counts to compare.
#' @param keep_fits Logical; if `TRUE`, retain the full fitted objects.
#'
#' @return An object of class `"qgcompmulti_b_sensitivity"`.
#'
#' @details
#' `B` sensitivity is implemented as a repeated-fit workflow in which all
#' settings other than `B` are held fixed.
#'
#' When `seed` is supplied, this helper treats it as a master seed and
#' deterministically derives one distinct fit-specific seed for each requested
#' value of `B`. This preserves reproducibility of the overall sensitivity
#' workflow while avoiding reuse of the same bootstrap resamples across the
#' different refits.
#'
#' The resulting object is designed for stability assessment rather than
#' automatic tuning. Users should focus especially on whether bootstrap-based
#' standard errors, bootstrap retention counts, and broad qualitative
#' conclusions are reasonably consistent across the requested `B` values.
#'
#' @examples
#' \dontrun{
#' dat <- sim_mixture_data(
#'   n = 400,
#'   pA = 3,
#'   pB = 3,
#'   rho_within_A = 0.3,
#'   rho_within_B = 0.3,
#'   rho_between = 0.2,
#'   psi1 = 0.5,
#'   psi2 = 0.3,
#'   psi12 = 0.2,
#'   seed = 123
#' )
#'
#' b_sensitivity(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   B_values = c(50, 100, 200),
#'   q = 4,
#'   MCsize = nrow(dat),
#'   seed = 13
#' )
#' }
#'
#' @export
b_sensitivity <- function(f,
                          data,
                          mix1,
                          mix2,
                          B_values,
                          interaction = TRUE,
                          family = gaussian(),
                          q = 4,
                          centering = "none",
                          id = NULL,
                          MCsize = nrow(data),
                          seed = NULL,
                          keep_fits = TRUE) {
  B_values <- qgcompmulti_b_values(B_values)
  fit_seeds <- qgcompmulti_sensitivity_refit_seeds(length(B_values), seed)
  fits <- lapply(seq_along(B_values), function(i) {
    qgcomp.glm.multi(
      f = f,
      data = data,
      mix1 = mix1,
      mix2 = mix2,
      interaction = interaction,
      family = family,
      q = q,
      centering = centering,
      B = B_values[[i]],
      id = id,
      MCsize = MCsize,
      seed = fit_seeds[[i]]
    )
  })
  results_table <- qgcompmulti_build_sensitivity_table(
    values = as.list(B_values),
    parameter_name = "B",
    fits = fits
  )
  structure(
    list(
      workflow_type = "b_sensitivity",
      call_template = match.call(),
      settings_fixed = list(
        interaction = interaction,
        family = family,
        q = q,
        centering = centering,
        id = id,
        MCsize = MCsize,
        seed = seed
      ),
      B_values = B_values,
      results_table = results_table,
      fits = if (isTRUE(keep_fits)) fits else NULL,
      notes = paste(
        "Only B varies across refits; all other analysis settings are held fixed.",
        "When `seed` is supplied, it is treated as a master seed and deterministically expanded into distinct fit-specific seeds so bootstrap draws are not reused across different values of B."
      )
    ),
    class = "qgcompmulti_b_sensitivity"
  )
}
#' Sensitivity to quantization choice in qgcompmulti fits
#'
#' Re-fits [qgcomp.glm.multi()] across multiple integer `q` values while
#' preserving the rest of the analysis specification.
#'
#' This helper is intended for robustness assessment, not for coefficient
#' ranking across different quantization choices.
#'
#' @param f,data,mix1,mix2,interaction,family,centering,B,id,MCsize,seed
#' Arguments passed through to [qgcomp.glm.multi()].
#' @param q_values Integer vector of quantization choices to compare.
#' @param keep_fits Logical; if `TRUE`, retain the full fitted objects.
#'
#' @return An object of class `"qgcompmulti_q_sensitivity"`.
#'
#' @details
#' `q` sensitivity is implemented as a repeated-fit workflow in which all
#' settings other than `q` are held fixed.
#'
#' Users should be cautious when comparing raw coefficient magnitudes across
#' different values of `q`. A larger `q` implies a smaller one-quantile
#' intervention step, so smaller coefficients may be expected mechanically even
#' when the broader qualitative pattern of the fitted surface is stable. For
#' that reason, the printed sensitivity object includes an explicit
#' comparability note.
#'
#' The helper therefore supports sensitivity assessment, not a claim that one
#' choice of `q` produces “stronger” or “weaker” effects based only on raw
#' one-step coefficient magnitudes.
#'
#' @examples
#' \dontrun{
#' dat <- sim_mixture_data(
#'   n = 400,
#'   pA = 3,
#'   pB = 3,
#'   rho_within_A = 0.3,
#'   rho_within_B = 0.3,
#'   rho_between = 0.2,
#'   psi1 = 0.5,
#'   psi2 = 0.3,
#'   psi12 = 0.2,
#'   seed = 123
#' )
#'
#' q_sensitivity(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   q_values = c(3, 4, 5),
#'   B = 100,
#'   seed = 13
#' )
#' }
#'
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
print.qgcompmulti_b_sensitivity <- function(x, ...) {
  cat("qgcompmulti B sensitivity\n\n")
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
