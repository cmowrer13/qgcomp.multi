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

    if (!identical(
      qgcompmulti_deparse_one_line(fit$formula),
      qgcompmulti_deparse_one_line(reference$formula)
    )) {
      stop("All per-imputation fits must use the same formula.", call. = FALSE)
    }

    if (!identical(fit$mixtures, reference$mixtures)) {
      stop("All per-imputation fits must use the same mixture definitions and coding settings.", call. = FALSE)
    }

    for (field_name in c("interaction", "family_name", "link", "B", "id", "MCsize")) {
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
