# ------------------------------------------------------------------------------
# Interval-method helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_interval_methods <- function(context = c("single_fit", "mi", "prediction")) {
  context <- match.arg(context)

  switch(
    context,
    single_fit = c("wald", "percentile", "percentile_type2"),
    mi = "wald",
    prediction = c("bootstrap_percentile", "bootstrap_percentile_type2")
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
    percentile_type2 = "bootstrap_percentile_type2"
  )
}
