# ------------------------------------------------------------------------------
# Object contract helpers
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
qgcompmulti_required_top_components <- function() {
  c(
    "call",
    "formula",
    "data_info",
    "mixtures",
    "analysis",
    "fits",
    "bootstrap",
    "results",
    "labels"
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_component_fields <- function() {
  list(
    data_info = c(
      "n_input",
      "n_used",
      "outcome",
      "has_clusters",
      "cluster_var",
      "n_clusters",
      "quantized"
    ),
    mixtures = c(
      "mix1",
      "mix2",
      "q",
      "centering"
    ),
    analysis = c(
      "interaction",
      "family",
      "family_name",
      "link",
      "B",
      "id",
      "MCsize",
      "seed"
    ),
    fits = c(
      "outcome_fit",
      "msm_fit"
    ),
    bootstrap = c(
      "coef_draws",
      "B_requested",
      "B_success"
    ),
    results = c(
      "coefficients",
      "std_error",
      "vcov",
      "coef_table"
    ),
    labels = c(
      "mixture_labels",
      "coefficient_names",
      "coefficient_labels"
    )
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_init_component <- function(component_name) {
  fields <- qgcompmulti_component_fields()[[component_name]]

  if (is.null(fields)) {
    stop(sprintf("Unknown qgcompmulti component: `%s`.", component_name), call. = FALSE)
  }

  setNames(vector("list", length(fields)), fields)
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_component <- function(x, component_name, required_fields) {
  if (!is.list(x)) {
    stop(sprintf("`%s` must be a list.", component_name), call. = FALSE)
  }

  actual_fields <- names(x)

  if (is.null(actual_fields) || !identical(actual_fields, required_fields)) {
    stop(
      sprintf(
        "`%s` must contain fields in this order: %s.",
        component_name,
        paste(required_fields, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(x)
}

# ------------------------------------------------------------------------------
# Coefficient naming and label helpers
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
qgcompmulti_coef_names <- function(interaction = TRUE) {
  if (!is.logical(interaction) || length(interaction) != 1L || is.na(interaction)) {
    stop("`interaction` must be either `TRUE` or `FALSE`.", call. = FALSE)
  }

  if (interaction) {
    c("(Intercept)", "psi1", "psi2", "psi1:psi2")
  } else {
    c("(Intercept)", "psi1", "psi2")
  }
}

#' @keywords internal
#' @noRd
qgcompmulti_coef_labels <- function(interaction = TRUE) {
  if (interaction) {
    c(
      "Intercept",
      "Mixture 1 main effect",
      "Mixture 2 main effect",
      "Mixture interaction"
    )
  } else {
    c(
      "Intercept",
      "Mixture 1 main effect",
      "Mixture 2 main effect"
    )
  }
}

#' @keywords internal
#' @noRd
qgcompmulti_default_mixture_labels <- function() {
  c(
    mix1 = "Mixture 1",
    mix2 = "Mixture 2"
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_labels <- function(interaction = TRUE,
                                     mixture_labels = qgcompmulti_default_mixture_labels()) {
  coef_names <- qgcompmulti_coef_names(interaction = interaction)
  coef_labels <- qgcompmulti_coef_labels(interaction = interaction)

  if (!is.character(mixture_labels) || !identical(names(mixture_labels), c("mix1", "mix2"))) {
    stop(
      "`mixture_labels` must be a named character vector with names `mix1` and `mix2`.",
      call. = FALSE
    )
  }

  list(
    mixture_labels = mixture_labels,
    coefficient_names = coef_names,
    coefficient_labels = setNames(coef_labels, coef_names)
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_standardize_coefficients <- function(coefficients, interaction = TRUE) {
  expected_names <- qgcompmulti_coef_names(interaction = interaction)

  if (!is.numeric(coefficients)) {
    stop("`coefficients` must be numeric.", call. = FALSE)
  }

  missing_names <- setdiff(expected_names, names(coefficients))

  if (length(missing_names) > 0L) {
    stop(
      sprintf(
        "Missing expected coefficient names: %s",
        paste(missing_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  coefficients <- coefficients[expected_names]

  setNames(as.numeric(coefficients), expected_names)
}

#' @keywords internal
#' @noRd
qgcompmulti_extract_msm_coefficients <- function(msm_fit, interaction = TRUE) {
  if (is.null(msm_fit) || !inherits(msm_fit, "glm")) {
    stop("`msm_fit` must be a fitted `glm` object.", call. = FALSE)
  }

  qgcompmulti_standardize_coefficients(
    coefficients = stats::coef(msm_fit),
    interaction = interaction
  )
}

# ------------------------------------------------------------------------------
# Bootstrap and core statistical result helpers
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
qgcompmulti_initialize_boot_draws <- function(B, interaction = TRUE) {
  coef_names <- qgcompmulti_coef_names(interaction = interaction)

  matrix(
    NA_real_,
    nrow = as.integer(B),
    ncol = length(coef_names),
    dimnames = list(NULL, coef_names)
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_resample_data <- function(data, id = NULL) {
  if (is.null(id)) {
    idx <- sample(seq_len(nrow(data)), replace = TRUE)
    return(data[idx, , drop = FALSE])
  }

  ids <- data[[id]]
  unique_ids <- unique(ids)
  samp_ids <- sample(unique_ids, length(unique_ids), replace = TRUE)
  split_data <- split(data, ids)

  do.call(
    rbind,
    split_data[as.character(samp_ids)]
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(eval.parent(substitute(expr)))
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit(
    {
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  set.seed(as.integer(seed))
  eval.parent(substitute(expr))
}

#' @keywords internal
#' @noRd
build_qgcompmulti_coef_table <- function(coefficients, std_error) {
  if (is.null(coefficients) || is.null(std_error)) {
    return(NULL)
  }

  if (!is.numeric(coefficients) || !is.numeric(std_error)) {
    stop("`coefficients` and `std_error` must be numeric vectors.", call. = FALSE)
  }

  if (length(coefficients) != length(std_error)) {
    stop("`coefficients` and `std_error` must have the same length.", call. = FALSE)
  }

  if (is.null(names(coefficients)) || is.null(names(std_error))) {
    stop("`coefficients` and `std_error` must both be named.", call. = FALSE)
  }

  if (!identical(names(coefficients), names(std_error))) {
    stop("`coefficients` and `std_error` must have identical names.", call. = FALSE)
  }

  z_value <- coefficients / std_error
  p_value <- 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)

  out <- data.frame(
    Estimate = unname(coefficients),
    `Std. Error` = unname(std_error),
    `z value` = unname(z_value),
    `Pr(>|z|)` = unname(p_value),
    row.names = names(coefficients),
    check.names = FALSE
  )

  out
}

#' @keywords internal
#' @noRd
build_qgcompmulti_vcov <- function(coef_draws, coef_names = colnames(coef_draws)) {
  if (is.null(coef_draws)) {
    return(NULL)
  }

  if (!is.matrix(coef_draws) && !is.data.frame(coef_draws)) {
    stop("`coef_draws` must be a matrix or data frame.", call. = FALSE)
  }

  coef_draws <- as.matrix(coef_draws)

  if (is.null(colnames(coef_draws))) {
    stop("`coef_draws` must have column names.", call. = FALSE)
  }

  if (is.null(coef_names)) {
    stop("`coef_names` must not be `NULL`.", call. = FALSE)
  }

  missing_names <- setdiff(coef_names, colnames(coef_draws))

  if (length(missing_names) > 0L) {
    stop(
      sprintf(
        "Missing expected bootstrap coefficient names: %s",
        paste(missing_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  coef_draws <- coef_draws[, coef_names, drop = FALSE]

  stats::cov(coef_draws)
}

#' @keywords internal
#' @noRd
build_qgcompmulti_std_error <- function(coef_draws, coef_names = colnames(coef_draws)) {
  vc <- build_qgcompmulti_vcov(
    coef_draws = coef_draws,
    coef_names = coef_names
  )

  if (is.null(vc)) {
    return(NULL)
  }

  se <- sqrt(diag(vc))
  setNames(as.numeric(se), names(se))
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_conf_level <- function(level) {
  if (!is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1) {
    stop(
      "`level` must be a single number strictly between 0 and 1.",
      call. = FALSE
    )
  }

  invisible(level)
}

#' @keywords internal
#' @noRd
qgcompmulti_resolve_parm <- function(parm, coef_names) {
  if (!is.character(coef_names) || length(coef_names) == 0L || anyNA(coef_names)) {
    stop("`coef_names` must be a non-empty character vector.", call. = FALSE)
  }

  if (is.null(parm)) {
    return(coef_names)
  }

  if (is.character(parm)) {
    if (anyNA(parm)) {
      stop("`parm` must not contain missing values.", call. = FALSE)
    }

    missing_names <- setdiff(parm, coef_names)

    if (length(missing_names) > 0L) {
      stop(
        sprintf(
          "Unknown coefficient names in `parm`: %s",
          paste(missing_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    return(parm)
  }

  if (is.numeric(parm)) {
    if (anyNA(parm) ||
        any(!is.finite(parm)) ||
        any(parm != as.integer(parm)) ||
        any(parm < 1L) ||
        any(parm > length(coef_names))) {
      stop(
        "`parm` must contain valid coefficient positions.",
        call. = FALSE
      )
    }

    return(coef_names[parm])
  }

  stop(
    "`parm` must be `NULL`, a character vector of coefficient names, or an integer vector of coefficient positions.",
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_confint_colnames <- function(level) {
  probs <- 100 * c((1 - level) / 2, 1 - (1 - level) / 2)
  labels <- trimws(format(probs, scientific = FALSE, trim = TRUE))
  labels <- sub("\\.?0+$", "", labels)

  paste0(labels, " %")
}

#' @keywords internal
#' @noRd
build_qgcompmulti_confint <- function(coefficients, std_error, level = 0.95) {
  if (!is.numeric(coefficients) || is.null(names(coefficients))) {
    stop(
      "`coefficients` must be a named numeric vector.",
      call. = FALSE
    )
  }

  if (!is.numeric(std_error) || is.null(names(std_error))) {
    stop(
      "`std_error` must be a named numeric vector.",
      call. = FALSE
    )
  }

  if (!identical(names(coefficients), names(std_error))) {
    stop(
      "`coefficients` and `std_error` must have identical names.",
      call. = FALSE
    )
  }

  qgcompmulti_validate_conf_level(level)

  zcrit <- stats::qnorm((1 + level) / 2)
  out <- cbind(
    unname(coefficients - zcrit * std_error),
    unname(coefficients + zcrit * std_error)
  )

  rownames(out) <- names(coefficients)
  colnames(out) <- qgcompmulti_confint_colnames(level)

  out
}

# ------------------------------------------------------------------------------
# Object assembly helpers
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
qgcompmulti_response_name <- function(formula) {
  deparse(formula[[2]])
}

#' @keywords internal
#' @noRd
build_qgcompmulti_data_info <- function(data, formula, q, id, n_used) {
  list(
    n_input = nrow(data),
    n_used = n_used,
    outcome = qgcompmulti_response_name(formula),
    has_clusters = !is.null(id),
    cluster_var = id,
    n_clusters = if (is.null(id)) NULL else length(unique(data[[id]])),
    quantized = !is.null(q)
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_mixtures <- function(mix1, mix2, q, centering) {
  list(
    mix1 = mix1,
    mix2 = mix2,
    q = q,
    centering = centering
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_analysis <- function(interaction, family, B, id, MCsize, seed) {
  list(
    interaction = interaction,
    family = family,
    family_name = family$family,
    link = family$link,
    B = B,
    id = id,
    MCsize = MCsize,
    seed = if (is.null(seed)) NULL else as.integer(seed)
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_fits <- function(full_fit) {
  list(
    outcome_fit = full_fit$outcome_fit,
    msm_fit = full_fit$msm_fit
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_bootstrap <- function(coef_draws, B_requested) {
  list(
    coef_draws = coef_draws,
    B_requested = as.integer(B_requested),
    B_success = nrow(coef_draws)
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_results <- function(coefficients, std_error, vcov) {
  list(
    coefficients = coefficients,
    std_error = std_error,
    vcov = vcov,
    coef_table = build_qgcompmulti_coef_table(coefficients, std_error)
  )
}

# ------------------------------------------------------------------------------
# Print/summary formatting helpers
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
qgcompmulti_summary_components <- function() {
  c(
    "call",
    "formula",
    "fit_overview",
    "mixtures",
    "msm_table",
    "outcome_model_info",
    "labels"
  )
}

#' @keywords internal
#' @noRd
validate_summary_qgcompmulti <- function(x) {
  if (!is.list(x)) {
    stop("`summary.qgcompmulti` objects must be lists.", call. = FALSE)
  }

  if (!identical(names(x), qgcompmulti_summary_components())) {
    stop(
      sprintf(
        "`summary.qgcompmulti` objects must contain components in this order: %s.",
        paste(qgcompmulti_summary_components(), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(x$call) && !is.call(x$call)) {
    stop("`call` must be `NULL` or a matched call.", call. = FALSE)
  }

  if (!is.null(x$formula) && !inherits(x$formula, "formula")) {
    stop("`formula` must be `NULL` or a formula.", call. = FALSE)
  }

  if (!is.list(x$fit_overview)) {
    stop("`fit_overview` must be a list.", call. = FALSE)
  }

  if (!is.list(x$mixtures)) {
    stop("`mixtures` must be a list.", call. = FALSE)
  }

  if (!is.data.frame(x$msm_table)) {
    stop("`msm_table` must be a data frame.", call. = FALSE)
  }

  if (!is.list(x$outcome_model_info)) {
    stop("`outcome_model_info` must be a list.", call. = FALSE)
  }

  qgcompmulti_validate_labels(x$labels)

  invisible(x)
}

#' @keywords internal
#' @noRd
build_qgcompmulti_summary <- function(object) {
  validate_qgcompmulti(object)

  summary_obj <- list(
    call = object$call,
    formula = object$formula,
    fit_overview = list(
      outcome = object$data_info$outcome,
      family = qgcompmulti_family_label(object$analysis),
      n_input = object$data_info$n_input,
      n_used = object$data_info$n_used,
      exposure_mode = qgcompmulti_exposure_mode_label(object$mixtures),
      interaction = qgcompmulti_interaction_label(object$analysis),
      bootstrap_requested = object$bootstrap$B_requested,
      bootstrap_success = object$bootstrap$B_success,
      MCsize = object$analysis$MCsize,
      seed = object$analysis$seed,
      has_clusters = object$data_info$has_clusters,
      cluster_var = object$data_info$cluster_var,
      n_clusters = object$data_info$n_clusters
    ),
    mixtures = list(
      labels = object$labels$mixture_labels,
      mix1 = object$mixtures$mix1,
      mix2 = object$mixtures$mix2
    ),
    msm_table = qgcompmulti_labeled_coef_table(object$results, object$labels),
    outcome_model_info = qgcompmulti_outcome_model_info(object$fits$outcome_fit),
    labels = object$labels
  )

  validate_summary_qgcompmulti(summary_obj)

  summary_obj
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

# ------------------------------------------------------------------------------
# Validation helpers for assembled objects
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_validate_labels <- function(labels) {
  coef_names <- labels$coefficient_names
  coef_labels <- labels$coefficient_labels
  mix_labels <- labels$mixture_labels

  if (!is.null(mix_labels)) {
    if (!is.character(mix_labels) || !identical(names(mix_labels), c("mix1", "mix2"))) {
      stop(
        "`labels$mixture_labels` must be a named character vector with names `mix1` and `mix2`.",
        call. = FALSE
      )
    }
  }

  if (!is.null(coef_names) && !is.character(coef_names)) {
    stop("`labels$coefficient_names` must be `NULL` or a character vector.", call. = FALSE)
  }

  if (!is.null(coef_labels)) {
    if (!is.character(coef_labels)) {
      stop("`labels$coefficient_labels` must be `NULL` or a character vector.", call. = FALSE)
    }

    if (is.null(names(coef_labels))) {
      stop("`labels$coefficient_labels` must be named.", call. = FALSE)
    }

    if (!is.null(coef_names) && !identical(names(coef_labels), coef_names)) {
      stop(
        "Names of `labels$coefficient_labels` must match `labels$coefficient_names`.",
        call. = FALSE
      )
    }
  }

  invisible(labels)
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_bootstrap <- function(bootstrap, labels) {
  coef_names <- labels$coefficient_names

  if (!is.null(bootstrap$coef_draws)) {
    if (!is.matrix(bootstrap$coef_draws) && !is.data.frame(bootstrap$coef_draws)) {
      stop("`bootstrap$coef_draws` must be `NULL`, a matrix, or a data frame.", call. = FALSE)
    }

    coef_draws <- as.matrix(bootstrap$coef_draws)

    if (is.null(colnames(coef_draws))) {
      stop("`bootstrap$coef_draws` must have column names.", call. = FALSE)
    }

    if (!is.null(coef_names) && !identical(colnames(coef_draws), coef_names)) {
      stop(
        "Column names of `bootstrap$coef_draws` must match `labels$coefficient_names`.",
        call. = FALSE
      )
    }
  }

  if (!is.null(bootstrap$B_requested) && length(bootstrap$B_requested) != 1L) {
    stop("`bootstrap$B_requested` must be length 1.", call. = FALSE)
  }

  if (!is.null(bootstrap$B_success) && length(bootstrap$B_success) != 1L) {
    stop("`bootstrap$B_success` must be length 1.", call. = FALSE)
  }

  if (!is.null(bootstrap$coef_draws) &&
      !is.null(bootstrap$B_success) &&
      nrow(as.matrix(bootstrap$coef_draws)) != bootstrap$B_success) {
    stop("`bootstrap$B_success` must match the number of retained bootstrap draws.", call. = FALSE)
  }

  invisible(bootstrap)
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_fits <- function(fits) {
  if (!is.null(fits$outcome_fit) && !inherits(fits$outcome_fit, "glm")) {
    stop("`fits$outcome_fit` must be `NULL` or a `glm` object.", call. = FALSE)
  }

  if (!is.null(fits$msm_fit) && !inherits(fits$msm_fit, "glm")) {
    stop("`fits$msm_fit` must be `NULL` or a `glm` object.", call. = FALSE)
  }

  invisible(fits)
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_results <- function(results, labels) {
  coef_names <- labels$coefficient_names

  if (!is.null(results$coefficients)) {
    if (!is.numeric(results$coefficients)) {
      stop("`results$coefficients` must be `NULL` or numeric.", call. = FALSE)
    }

    if (is.null(names(results$coefficients))) {
      stop("`results$coefficients` must be named.", call. = FALSE)
    }

    if (!is.null(coef_names) && !identical(names(results$coefficients), coef_names)) {
      stop(
        "Names of `results$coefficients` must match `labels$coefficient_names`.",
        call. = FALSE
      )
    }
  }

  if (!is.null(results$std_error)) {
    if (!is.numeric(results$std_error)) {
      stop("`results$std_error` must be `NULL` or numeric.", call. = FALSE)
    }

    if (is.null(names(results$std_error))) {
      stop("`results$std_error` must be named.", call. = FALSE)
    }

    if (!is.null(coef_names) && !identical(names(results$std_error), coef_names)) {
      stop(
        "Names of `results$std_error` must match `labels$coefficient_names`.",
        call. = FALSE
      )
    }
  }

  if (!is.null(results$vcov)) {
    if (!is.matrix(results$vcov)) {
      stop("`results$vcov` must be `NULL` or a matrix.", call. = FALSE)
    }

    rn <- rownames(results$vcov)
    cn <- colnames(results$vcov)

    if (is.null(rn) || is.null(cn) || !identical(rn, cn)) {
      stop("`results$vcov` must have identical row and column names.", call. = FALSE)
    }

    if (!is.null(coef_names) && !identical(rn, coef_names)) {
      stop(
        "Row and column names of `results$vcov` must match `labels$coefficient_names`.",
        call. = FALSE
      )
    }
  }

  if (!is.null(results$coef_table)) {
    if (!is.data.frame(results$coef_table)) {
      stop("`results$coef_table` must be `NULL` or a data frame.", call. = FALSE)
    }

    if (!is.null(coef_names) && !identical(rownames(results$coef_table), coef_names)) {
      stop(
        "Row names of `results$coef_table` must match `labels$coefficient_names`.",
        call. = FALSE
      )
    }
  }

  invisible(results)
}
