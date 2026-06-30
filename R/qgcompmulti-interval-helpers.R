# ------------------------------------------------------------------------------
# Interval-method helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_interval_methods <- function(context = c("single_fit", "mi", "prediction")) {
  context <- match.arg(context)

  switch(
    context,
    single_fit = c("wald", "percentile", "basic"),
    mi = "wald",
    prediction = c("bootstrap_percentile", "bootstrap_basic")
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_interval_method <- function(method,
                                                 context = c("single_fit", "mi"),
                                                 allow_null = FALSE) {
  context <- match.arg(context)

  if (is.null(method)) {
    if (isTRUE(allow_null)) {
      return(invisible(NULL))
    }
    stop("`method` must not be `NULL`.", call. = FALSE)
  }

  allowed <- qgcompmulti_interval_methods(context = context)

  if (!is.character(method) ||
      length(method) != 1L ||
      is.na(method) ||
      !method %in% allowed) {
    stop(
      sprintf(
        "`method` must be one of: %s.",
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(method)
}

#' @keywords internal
#' @noRd
qgcompmulti_resolve_interval_method <- function(method,
                                                default_method,
                                                context = c("single_fit", "mi")) {
  context <- match.arg(context)
  if (is.null(method)) {
    method <- default_method
  }
  if (identical(context, "mi") && !identical(method, "wald")) {
    qgcompmulti_require_wald_interval_method(
      method = method,
      object_label = "pooled qgcompmulti multiple-imputation coefficient"
    )
  }
  qgcompmulti_validate_interval_method(
    method = method,
    context = context,
    allow_null = FALSE
  )
  method
}
#' @keywords internal
#' @noRd
qgcompmulti_require_wald_interval_method <- function(method,
                                                     object_label = "qgcompmulti") {
  if (!identical(method, "wald")) {
    stop(
      sprintf(
        "%s reporting currently supports only `method = \"wald\"`; `method = \"%s\"` is not implemented in this phase.",
        object_label,
        method
      ),
      call. = FALSE
    )
  }
  invisible(method)
}

#' @keywords internal
#' @noRd
qgcompmulti_prediction_interval_methods <- function() {
  setdiff(qgcompmulti_interval_methods(context = "single_fit"), "wald")
}
#' @keywords internal
#' @noRd
qgcompmulti_resolve_prediction_interval_method <- function(method,
                                                           default_method) {
  if (is.null(method)) {
    method <- default_method
    if (identical(method, "wald")) {
      method <- "percentile"
    }
  }
  qgcompmulti_validate_interval_method(
    method = method,
    context = "single_fit",
    allow_null = FALSE
  )
  if (identical(method, "wald")) {
    stop(
      "MSM prediction intervals support `method = \"percentile\"` or `method = \"basic\"`; Wald intervals are currently coefficient-only.",
      call. = FALSE
    )
  }
  method
}

#' @keywords internal
#' @noRd
qgcompmulti_prediction_interval_types <- function() {
  qgcompmulti_interval_methods(context = "prediction")
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_prediction_interval_type <- function(interval_type,
                                                          allow_null = TRUE) {
  if (is.null(interval_type)) {
    if (isTRUE(allow_null)) {
      return(invisible(NULL))
    }
    stop("`interval_type` must not be `NULL`.", call. = FALSE)
  }

  allowed <- qgcompmulti_prediction_interval_types()

  if (!is.character(interval_type) ||
      length(interval_type) != 1L ||
      is.na(interval_type) ||
      !interval_type %in% allowed) {
    stop(
      sprintf(
        "`interval_type` must be `NULL` or one of: %s.",
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(interval_type)
}

#' @keywords internal
#' @noRd
qgcompmulti_prediction_interval_type_for_method <- function(method) {
  qgcompmulti_validate_interval_method(
    method = method,
    context = "single_fit",
    allow_null = FALSE
  )

  switch(
    method,
    wald = NULL,
    percentile = "bootstrap_percentile",
    basic = "bootstrap_basic"
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_bootstrap_interval_limits <- function(draws,
                                                  estimates,
                                                  level = 0.95,
                                                  method = c("percentile", "basic"),
                                                  margin = c("row", "column")) {
  method <- match.arg(method)
  margin <- match.arg(margin)
  qgcompmulti_validate_conf_level(level)
  draws <- as.matrix(draws)
  estimates <- as.numeric(estimates)
  expected_length <- if (identical(margin, "row")) nrow(draws) else ncol(draws)
  if (length(estimates) != expected_length) {
    stop("`estimates` must align with the requested bootstrap interval margin.", call. = FALSE)
  }
  alpha <- (1 - level) / 2
  probs <- c(alpha, 1 - alpha)
  quantile_fun <- function(x) {
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }
  quantiles <- if (identical(margin, "row")) {
    t(apply(draws, 1, quantile_fun))
  } else {
    t(apply(draws, 2, quantile_fun))
  }
  if (identical(method, "percentile")) {
    limits <- quantiles
  } else {
    limits <- cbind(
      2 * estimates - quantiles[, 2],
      2 * estimates - quantiles[, 1]
    )
  }
  colnames(limits) <- c("lower", "upper")
  limits
}
#' @keywords internal
#' @noRd
qgcompmulti_build_single_fit_confint <- function(coefficients,
                                                 std_error,
                                                 coef_draws = NULL,
                                                 level = 0.95,
                                                 method = c("wald", "percentile", "basic")) {
  method <- match.arg(method)
  if (identical(method, "wald")) {
    return(build_qgcompmulti_confint(
      coefficients = coefficients,
      std_error = std_error,
      level = level
    ))
  }
  if (is.null(coef_draws)) {
    stop("Stored bootstrap coefficient draws are required for bootstrap interval methods.", call. = FALSE)
  }
  if (!is.numeric(coefficients) || is.null(names(coefficients))) {
    stop("`coefficients` must be a named numeric vector.", call. = FALSE)
  }
  draws <- as.matrix(coef_draws)
  missing_coef <- setdiff(names(coefficients), colnames(draws))
  if (length(missing_coef) > 0L) {
    stop(
      sprintf(
        "Stored bootstrap coefficient draws are missing coefficient(s): %s.",
        paste(missing_coef, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  draws <- draws[, names(coefficients), drop = FALSE]
  limits <- qgcompmulti_bootstrap_interval_limits(
    draws = draws,
    estimates = coefficients,
    level = level,
    method = method,
    margin = "column"
  )
  out <- cbind(limits[, "lower"], limits[, "upper"])
  rownames(out) <- names(coefficients)
  colnames(out) <- qgcompmulti_confint_colnames(level)
  out
}
