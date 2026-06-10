#' @keywords internal
#' @noRd
qgcompmulti_mi_required_top_components <- function() {
  c(
    "call",
    "formula",
    "data_info",
    "mixtures",
    "mi_info",
    "analysis",
    "fits",
    "results",
    "labels"
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_component_fields <- function() {
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
    mi_info = c(
      "m",
      "input_type",
      "keep_fits",
      "seed",
      "fit_seeds",
      "n_input_per_imputation",
      "n_used_per_imputation"
    ),
    analysis = c(
      "interaction",
      "family",
      "family_name",
      "link",
      "B",
      "id",
      "MCsize"
    ),
    fits = c(
      "imputation_fits"
    ),
    results = c(
      "coefficients",
      "std_error",
      "vcov",
      "coef_table",
      "within_var",
      "between_var",
      "total_var",
      "df",
      "riv",
      "fmi"
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
qgcompmulti_mi_init_component <- function(component_name) {
  fields <- qgcompmulti_mi_component_fields()[[component_name]]
  if (is.null(fields)) {
    stop(sprintf("Unknown qgcompmulti_mi component: `%s`.", component_name), call. = FALSE)
  }
  setNames(vector("list", length(fields)), fields)
}

#' @keywords internal
#' @noRd
new_qgcompmulti_mi <- function(call = NULL,
                               formula = NULL,
                               data_info = qgcompmulti_mi_init_component("data_info"),
                               mixtures = qgcompmulti_mi_init_component("mixtures"),
                               mi_info = qgcompmulti_mi_init_component("mi_info"),
                               analysis = qgcompmulti_mi_init_component("analysis"),
                               fits = qgcompmulti_mi_init_component("fits"),
                               results = qgcompmulti_mi_init_component("results"),
                               labels = qgcompmulti_mi_init_component("labels")) {
  object <- list(
    call = call,
    formula = formula,
    data_info = data_info,
    mixtures = mixtures,
    mi_info = mi_info,
    analysis = analysis,
    fits = fits,
    results = results,
    labels = labels
  )

  validate_qgcompmulti_mi(object)

  structure(object, class = "qgcompmulti_mi")
}

#' @keywords internal
#' @noRd
validate_qgcompmulti_mi <- function(object) {
  if (!is.list(object)) {
    stop("`object` must be a list.", call. = FALSE)
  }

  required_top <- qgcompmulti_mi_required_top_components()
  actual_top <- names(object)

  if (is.null(actual_top) || !identical(actual_top, required_top)) {
    stop(
      sprintf(
        "`qgcompmulti_mi` objects must contain top-level components in this order: %s.",
        paste(required_top, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(object$call) && !is.call(object$call)) {
    stop("`call` must be `NULL` or a matched call.", call. = FALSE)
  }

  if (!is.null(object$formula) && !inherits(object$formula, "formula")) {
    stop("`formula` must be `NULL` or a model formula.", call. = FALSE)
  }

  component_fields <- qgcompmulti_mi_component_fields()

  for (component_name in names(component_fields)) {
    qgcompmulti_validate_component(
      x = object[[component_name]],
      component_name = component_name,
      required_fields = component_fields[[component_name]]
    )
  }

  qgcompmulti_validate_labels(object$labels)
  qgcompmulti_validate_mi_info(object$mi_info)
  qgcompmulti_validate_mi_fits(object$fits, object$mi_info)
  qgcompmulti_validate_mi_results(object$results, object$labels)

  invisible(object)
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_mi_info <- function(mi_info) {
  if (!is_scalar_whole_number(mi_info$m) || mi_info$m < 1L) {
    stop("`mi_info$m` must be a single integer greater than or equal to 1.", call. = FALSE)
  }

  if (!is.character(mi_info$input_type) ||
      length(mi_info$input_type) != 1L ||
      is.na(mi_info$input_type) ||
      !mi_info$input_type %in% c("completed_list", "mids")) {
    stop("`mi_info$input_type` must be either \"completed_list\" or \"mids\".", call. = FALSE)
  }

  if (!is.logical(mi_info$keep_fits) || length(mi_info$keep_fits) != 1L || is.na(mi_info$keep_fits)) {
    stop("`mi_info$keep_fits` must be either `TRUE` or `FALSE`.", call. = FALSE)
  }

  if (!is.null(mi_info$seed) && !is_scalar_whole_number(mi_info$seed)) {
    stop("`mi_info$seed` must be `NULL` or a single integer.", call. = FALSE)
  }

  if (!is.null(mi_info$fit_seeds)) {
    if (!is.numeric(mi_info$fit_seeds) ||
        anyNA(mi_info$fit_seeds) ||
        any(mi_info$fit_seeds != as.integer(mi_info$fit_seeds)) ||
        length(mi_info$fit_seeds) != mi_info$m) {
      stop(
        "`mi_info$fit_seeds` must be `NULL` or an integer vector with one entry per imputation.",
        call. = FALSE
      )
    }
  }

  if (!is.null(mi_info$n_input_per_imputation)) {
    if (!is.numeric(mi_info$n_input_per_imputation) ||
        anyNA(mi_info$n_input_per_imputation) ||
        any(mi_info$n_input_per_imputation != as.integer(mi_info$n_input_per_imputation)) ||
        any(mi_info$n_input_per_imputation < 1L) ||
        length(mi_info$n_input_per_imputation) != mi_info$m) {
      stop(
        "`mi_info$n_input_per_imputation` must be `NULL` or a positive integer vector with one entry per imputation.",
        call. = FALSE
      )
    }
  }

  if (!is.null(mi_info$n_used_per_imputation)) {
    if (!is.numeric(mi_info$n_used_per_imputation) ||
        anyNA(mi_info$n_used_per_imputation) ||
        any(mi_info$n_used_per_imputation != as.integer(mi_info$n_used_per_imputation)) ||
        any(mi_info$n_used_per_imputation < 1L) ||
        length(mi_info$n_used_per_imputation) != mi_info$m) {
      stop(
        "`mi_info$n_used_per_imputation` must be `NULL` or a positive integer vector with one entry per imputation.",
        call. = FALSE
      )
    }
  }

  invisible(mi_info)
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_mi_fits <- function(fits, mi_info) {
  stored_fits <- fits$imputation_fits

  if (isTRUE(mi_info$keep_fits)) {
    if (!is.list(stored_fits) || length(stored_fits) != mi_info$m) {
      stop(
        "When `mi_info$keep_fits` is `TRUE`, `fits$imputation_fits` must be a list with one fitted object per imputation.",
        call. = FALSE
      )
    }
    invalid_idx <- which(!vapply(stored_fits, inherits, logical(1), what = "qgcompmulti"))
    if (length(invalid_idx) > 0L) {
      stop(
        sprintf(
          "All retained fits must inherit from `qgcompmulti`. Problem at position(s): %s",
          paste(invalid_idx, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  } else if (!is.null(stored_fits)) {
    stop("`fits$imputation_fits` must be `NULL` when `mi_info$keep_fits` is `FALSE`.", call. = FALSE)
  }

  invisible(fits)
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_mi_results <- function(results, labels) {
  coef_names <- labels$coefficient_names

  if (!is.null(results$coefficients)) {
    if (!is.numeric(results$coefficients) || is.null(names(results$coefficients))) {
      stop("`results$coefficients` must be `NULL` or a named numeric vector.", call. = FALSE)
    }
    if (!is.null(coef_names) && !identical(names(results$coefficients), coef_names)) {
      stop(
        "Names of `results$coefficients` must match `labels$coefficient_names`.",
        call. = FALSE
      )
    }
  }

  for (field_name in c("std_error", "df", "riv", "fmi")) {
    value <- results[[field_name]]
    if (!is.null(value)) {
      if (!is.numeric(value) || is.null(names(value))) {
        stop(sprintf("`results$%s` must be `NULL` or a named numeric vector.", field_name), call. = FALSE)
      }
      if (!is.null(coef_names) && !identical(names(value), coef_names)) {
        stop(
          sprintf("Names of `results$%s` must match `labels$coefficient_names`.", field_name),
          call. = FALSE
        )
      }
    }
  }

  for (field_name in c("vcov", "within_var", "between_var", "total_var")) {
    value <- results[[field_name]]
    if (!is.null(value)) {
      if (!is.matrix(value)) {
        stop(sprintf("`results$%s` must be `NULL` or a matrix.", field_name), call. = FALSE)
      }
      rn <- rownames(value)
      cn <- colnames(value)
      if (is.null(rn) || is.null(cn) || !identical(rn, cn)) {
        stop(sprintf("`results$%s` must have identical row and column names.", field_name), call. = FALSE)
      }
      if (!is.null(coef_names) && !identical(rn, coef_names)) {
        stop(
          sprintf("Row and column names of `results$%s` must match `labels$coefficient_names`.", field_name),
          call. = FALSE
        )
      }
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
