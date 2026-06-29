#' @keywords internal
#' @noRd
qgcompmulti_mi_extract_coefficients <- function(fit) {
  validate_qgcompmulti(fit)
  coef(fit)
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_extract_vcov <- function(fit) {
  validate_qgcompmulti(fit)
  vcov(fit)
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_validate_wald_interval_method <- function(method) {
  if (!identical(method, "wald")) {
    stop(
      paste(
        "Pooled multiple-imputation fits support only Wald interval methods in Version 0.5.0.",
        "Per-imputation fits must use `analysis$default_interval_method = \"wald\"`."
      ),
      call. = FALSE
    )
  }

    qgcompmulti_validate_interval_method(method = method, context = "mi")
    invisible(method)
  }

#' @keywords internal
#' @noRd
qgcompmulti_mi_validate_fit_list <- function(fits) {
  if (!is.list(fits) || length(fits) == 0L) {
    stop("`fits` must be a non-empty list of `qgcompmulti` objects.", call. = FALSE)
  }

  invalid_idx <- which(!vapply(fits, inherits, logical(1), what = "qgcompmulti"))
  if (length(invalid_idx) > 0L) {
    stop(
      sprintf(
        "All elements of `fits` must inherit from `qgcompmulti`. Problem at position(s): %s",
        paste(invalid_idx, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  reference <- fits[[1]]
  coef_names <- names(coef(reference))

  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    validate_qgcompmulti(fit)
    qgcompmulti_mi_validate_wald_interval_method(fit$analysis$default_interval_method)

    if (!identical(
      qgcompmulti_deparse_one_line(fit$formula),
      qgcompmulti_deparse_one_line(reference$formula)
    )) {
      stop("All per-imputation fits must use the same formula.", call. = FALSE)
    }

    if (!identical(fit$mixtures, reference$mixtures)) {
      stop("All per-imputation fits must use the same mixture definitions and coding settings.", call. = FALSE)
    }

    for (field_name in c(
      "interaction",
      "family_name",
      "link",
      "estimand_scale",
      "estimand_scale_defaulted",
      "msm_fitting_scale",
      "default_interval_method",
      "B",
      "id",
      "MCsize"
    )) {
      if (!identical(fit$analysis[[field_name]], reference$analysis[[field_name]])) {
        stop(
          sprintf("All per-imputation fits must agree on `analysis$%s`.", field_name),
          call. = FALSE
        )
      }
    }

    for (field_name in c("n_input", "n_used", "outcome", "has_clusters", "cluster_var", "n_clusters", "quantized")) {
      if (!identical(fit$data_info[[field_name]], reference$data_info[[field_name]])) {
        stop(
          sprintf("All per-imputation fits must agree on `data_info$%s`.", field_name),
          call. = FALSE
        )
      }
    }

    if (!identical(names(coef(fit)), coef_names)) {
      stop("All per-imputation fits must target the same MSM coefficient names.", call. = FALSE)
    }

    if (!identical(fit$labels, reference$labels)) {
      stop("All per-imputation fits must use the same coefficient and mixture labels.", call. = FALSE)
    }
  }

  invisible(fits)
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_named_matrix_mean <- function(matrices) {
  Reduce(`+`, matrices) / length(matrices)
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_diag_metric <- function(within_diag, between_diag, total_diag, m) {
  riv <- ifelse(
    between_diag <= 0,
    0,
    ifelse(within_diag <= 0, Inf, ((1 + 1 / m) * between_diag) / within_diag)
  )

  df <- ifelse(
    between_diag <= 0,
    Inf,
    (m - 1) * (1 + 1 / riv)^2
  )

  fmi <- ifelse(total_diag <= 0, 0, ((1 + 1 / m) * between_diag) / total_diag)

  list(riv = riv, df = df, fmi = fmi)
}

#' @keywords internal
#' @noRd
qgcompmulti_pool_rubin <- function(coefficients, vcovs) {
  if (!is.list(coefficients) || length(coefficients) == 0L) {
    stop("`coefficients` must be a non-empty list of named numeric vectors.", call. = FALSE)
  }
  if (!is.list(vcovs) || length(vcovs) != length(coefficients)) {
    stop("`vcovs` must be a list of matrices with the same length as `coefficients`.", call. = FALSE)
  }

  m <- length(coefficients)
  coef_names <- names(coefficients[[1]])

  if (!is.numeric(coefficients[[1]]) || is.null(coef_names)) {
    stop("Each element of `coefficients` must be a named numeric vector.", call. = FALSE)
  }

  for (i in seq_along(coefficients)) {
    if (!is.numeric(coefficients[[i]]) || !identical(names(coefficients[[i]]), coef_names)) {
      stop("All coefficient vectors must be numeric and share identical names.", call. = FALSE)
    }
    if (!is.matrix(vcovs[[i]]) ||
        !identical(rownames(vcovs[[i]]), coef_names) ||
        !identical(colnames(vcovs[[i]]), coef_names)) {
      stop("All covariance matrices must align with the pooled coefficient names.", call. = FALSE)
    }
  }

  coef_mat <- do.call(
    rbind,
    lapply(coefficients, function(x) matrix(unname(x), nrow = 1L, dimnames = list(NULL, coef_names)))
  )

  qbar <- setNames(colMeans(coef_mat), coef_names)
  within_var <- qgcompmulti_mi_named_matrix_mean(vcovs)
  between_var <- if (m == 1L) {
    matrix(0, nrow = length(coef_names), ncol = length(coef_names), dimnames = list(coef_names, coef_names))
  } else {
    stats::cov(coef_mat)
  }
  total_var <- within_var + (1 + 1 / m) * between_var
  std_error <- setNames(sqrt(diag(total_var)), coef_names)

  metrics <- qgcompmulti_mi_diag_metric(
    within_diag = diag(within_var),
    between_diag = diag(between_var),
    total_diag = diag(total_var),
    m = m
  )

  list(
    coefficients = qbar,
    std_error = std_error,
    vcov = total_var,
    within_var = within_var,
    between_var = between_var,
    total_var = total_var,
    df = setNames(as.numeric(metrics$df), coef_names),
    riv = setNames(as.numeric(metrics$riv), coef_names),
    fmi = setNames(as.numeric(metrics$fmi), coef_names)
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_mi_coef_table <- function(coefficients, std_error, df) {
  if (!identical(names(coefficients), names(std_error)) || !identical(names(coefficients), names(df))) {
    stop("`coefficients`, `std_error`, and `df` must have identical names.", call. = FALSE)
  }

  statistic <- coefficients / std_error
  p_value <- vapply(
    seq_along(statistic),
    function(i) {
      if (is.finite(df[[i]])) {
        2 * stats::pt(abs(statistic[[i]]), df = df[[i]], lower.tail = FALSE)
      } else {
        2 * stats::pnorm(abs(statistic[[i]]), lower.tail = FALSE)
      }
    },
    numeric(1)
  )

  data.frame(
    Estimate = unname(coefficients),
    `Std. Error` = unname(std_error),
    `t value` = unname(statistic),
    df = unname(df),
    `Pr(>|t|)` = unname(p_value),
    row.names = names(coefficients),
    check.names = FALSE
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_mi_data_info <- function(fits) {
  reference <- fits[[1]]
  list(
    n_input = reference$data_info$n_input,
    n_used = reference$data_info$n_used,
    outcome = reference$data_info$outcome,
    has_clusters = reference$data_info$has_clusters,
    cluster_var = reference$data_info$cluster_var,
    n_clusters = reference$data_info$n_clusters,
    quantized = reference$data_info$quantized
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_mi_mixtures <- function(fits) {
  fits[[1]]$mixtures
}

#' @keywords internal
#' @noRd
build_qgcompmulti_mi_analysis <- function(fits) {
  reference <- fits[[1]]
  list(
    interaction = reference$analysis$interaction,
    family = reference$analysis$family,
    family_name = reference$analysis$family_name,
    link = reference$analysis$link,
    estimand_scale = reference$analysis$estimand_scale,
    estimand_scale_defaulted = reference$analysis$estimand_scale_defaulted,
    msm_fitting_scale = reference$analysis$msm_fitting_scale,
    default_interval_method = "wald",
    B = reference$analysis$B,
    id = reference$analysis$id,
    MCsize = reference$analysis$MCsize
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_mi_info <- function(fits,
                                      input_type = "completed_list",
                                      keep_fits = FALSE,
                                      seed = NULL,
                                      fit_seeds = NULL) {
  list(
    m = length(fits),
    input_type = input_type,
    keep_fits = keep_fits,
    seed = if (is.null(seed)) NULL else as.integer(seed),
    fit_seeds = if (is.null(fit_seeds)) NULL else as.integer(fit_seeds),
    n_input_per_imputation = as.integer(vapply(fits, function(x) x$data_info$n_input, numeric(1))),
    n_used_per_imputation = as.integer(vapply(fits, function(x) x$data_info$n_used, numeric(1)))
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_mi_results <- function(pooled) {
  list(
    coefficients = pooled$coefficients,
    std_error = pooled$std_error,
    vcov = pooled$vcov,
    coef_table = build_qgcompmulti_mi_coef_table(
      coefficients = pooled$coefficients,
      std_error = pooled$std_error,
      df = pooled$df
    ),
    within_var = pooled$within_var,
    between_var = pooled$between_var,
    total_var = pooled$total_var,
    df = pooled$df,
    riv = pooled$riv,
    fmi = pooled$fmi
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_pool_mi_fits <- function(fits,
                                     input_type = "completed_list",
                                     keep_fits = FALSE,
                                     seed = NULL,
                                     fit_seeds = NULL,
                                     call = NULL,
                                     formula = NULL) {
  qgcompmulti_mi_validate_fit_list(fits)

  pooled <- qgcompmulti_pool_rubin(
    coefficients = lapply(fits, qgcompmulti_mi_extract_coefficients),
    vcovs = lapply(fits, qgcompmulti_mi_extract_vcov)
  )

  reference <- fits[[1]]

  new_qgcompmulti_mi(
    call = call,
    formula = if (is.null(formula)) reference$formula else formula,
    data_info = build_qgcompmulti_mi_data_info(fits),
    mixtures = build_qgcompmulti_mi_mixtures(fits),
    mi_info = build_qgcompmulti_mi_info(
      fits = fits,
      input_type = input_type,
      keep_fits = keep_fits,
      seed = seed,
      fit_seeds = fit_seeds
    ),
    analysis = build_qgcompmulti_mi_analysis(fits),
    fits = list(
      imputation_fits = if (isTRUE(keep_fits)) fits else NULL
    ),
    results = build_qgcompmulti_mi_results(pooled),
    labels = reference$labels
  )
}

qgcompmulti_mi_input_label <- function(input_type) {
  switch(
    input_type,
    completed_list = "Completed data list",
    mids = "mice::mids",
    input_type
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_p_values <- function(statistic, df) {
  if (!is.numeric(statistic) || !is.numeric(df) || !identical(names(statistic), names(df))) {
    stop("`statistic` and `df` must be named numeric vectors with identical names.", call. = FALSE)
  }

  setNames(
    vapply(
      seq_along(statistic),
      function(i) {
        if (is.finite(df[[i]]) && df[[i]] > 0) {
          2 * stats::pt(abs(statistic[[i]]), df = df[[i]], lower.tail = FALSE)
        } else {
          2 * stats::pnorm(abs(statistic[[i]]), lower.tail = FALSE)
        }
      },
      numeric(1)
    ),
    names(statistic)
  )
}

#' @keywords internal
#' @noRd
build_qgcompmulti_mi_confint <- function(coefficients, std_error, df, level = 0.95) {
  if (!is.numeric(coefficients) || is.null(names(coefficients))) {
    stop("`coefficients` must be a named numeric vector.", call. = FALSE)
  }
  if (!is.numeric(std_error) || is.null(names(std_error))) {
    stop("`std_error` must be a named numeric vector.", call. = FALSE)
  }
  if (!is.numeric(df) || is.null(names(df))) {
    stop("`df` must be a named numeric vector.", call. = FALSE)
  }
  if (!identical(names(coefficients), names(std_error)) ||
      !identical(names(coefficients), names(df))) {
    stop("`coefficients`, `std_error`, and `df` must have identical names.", call. = FALSE)
  }

  qgcompmulti_validate_conf_level(level)

  crit <- vapply(
    seq_along(df),
    function(i) {
      if (is.finite(df[[i]]) && df[[i]] > 0) {
        stats::qt((1 + level) / 2, df = df[[i]])
      } else {
        stats::qnorm((1 + level) / 2)
      }
    },
    numeric(1)
  )

  out <- cbind(
    unname(coefficients - crit * std_error),
    unname(coefficients + crit * std_error)
  )
  rownames(out) <- names(coefficients)
  colnames(out) <- qgcompmulti_confint_colnames(level)
  out
}

#' @keywords internal
#' @noRd
build_qgcompmulti_mi_pooling_table <- function(object) {
  validate_qgcompmulti_mi(object)

  coef_names <- object$labels$coefficient_names
  out <- data.frame(
    df = unname(object$results$df[coef_names]),
    RIV = unname(object$results$riv[coef_names]),
    FMI = unname(object$results$fmi[coef_names]),
    row.names = coef_names,
    check.names = FALSE
  )
  rownames(out) <- unname(object$labels$coefficient_labels[coef_names])
  out
}

#' @keywords internal
#' @noRd
qgcompmulti_print_pooling_table <- function(pooling_table) {
  if (is.null(pooling_table) || nrow(pooling_table) == 0L) {
    cat("Pooling diagnostics: none available\n")
    return(invisible(NULL))
  }

  cat("Pooling diagnostics:\n")
  stats::printCoefmat(
    as.matrix(pooling_table),
    P.values = FALSE,
    has.Pvalue = FALSE,
    signif.stars = FALSE
  )
  invisible(pooling_table)
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_summary_components <- function() {
  c(
    "call",
    "formula",
    "fit_overview",
    "mi_overview",
    "mixtures",
    "msm_table",
    "pooling_table",
    "labels"
  )
}

#' @keywords internal
#' @noRd
validate_summary_qgcompmulti_mi <- function(x) {
  if (!is.list(x)) {
    stop("`summary.qgcompmulti_mi` objects must be lists.", call. = FALSE)
  }
  if (!identical(names(x), qgcompmulti_mi_summary_components())) {
    stop(
      sprintf(
        "`summary.qgcompmulti_mi` objects must contain components in this order: %s.",
        paste(qgcompmulti_mi_summary_components(), collapse = ", ")
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
  if (!is.list(x$mi_overview)) {
    stop("`mi_overview` must be a list.", call. = FALSE)
  }
  if (!is.list(x$mixtures)) {
    stop("`mixtures` must be a list.", call. = FALSE)
  }
  if (!is.data.frame(x$msm_table)) {
    stop("`msm_table` must be a data frame.", call. = FALSE)
  }
  if (!is.data.frame(x$pooling_table)) {
    stop("`pooling_table` must be a data frame.", call. = FALSE)
  }
  qgcompmulti_validate_labels(x$labels)
  invisible(x)
}

#' @keywords internal
#' @noRd
build_qgcompmulti_mi_summary <- function(object) {
  validate_qgcompmulti_mi(object)

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
      bootstrap_requested = object$analysis$B,
      MCsize = object$analysis$MCsize,
      has_clusters = object$data_info$has_clusters,
      cluster_var = object$data_info$cluster_var,
      n_clusters = object$data_info$n_clusters
    ),
    mi_overview = list(
      m = object$mi_info$m,
      input_type = qgcompmulti_mi_input_label(object$mi_info$input_type),
      keep_fits = object$mi_info$keep_fits,
      seed = object$mi_info$seed,
      fit_seed_count = if (is.null(object$mi_info$fit_seeds)) 0L else length(object$mi_info$fit_seeds)
    ),
    mixtures = list(
      labels = object$labels$mixture_labels,
      mix1 = object$mixtures$mix1,
      mix2 = object$mixtures$mix2
    ),
    msm_table = qgcompmulti_labeled_coef_table(object$results, object$labels),
    pooling_table = build_qgcompmulti_mi_pooling_table(object),
    labels = object$labels
  )

  validate_summary_qgcompmulti_mi(summary_obj)
  summary_obj
}
