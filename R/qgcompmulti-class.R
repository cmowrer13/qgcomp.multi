#' @keywords internal
#' @noRd

new_qgcompmulti <- function(call = NULL,
                            formula = NULL,
                            data_info = qgcompmulti_init_component("data_info"),
                            mixtures = qgcompmulti_init_component("mixtures"),
                            analysis = qgcompmulti_init_component("analysis"),
                            fits = qgcompmulti_init_component("fits"),
                            prediction = qgcompmulti_init_component("prediction"),
                            bootstrap = qgcompmulti_init_component("bootstrap"),
                            results = qgcompmulti_init_component("results"),
                            labels = qgcompmulti_init_component("labels")) {

  object <- list(
    call = call,
    formula = formula,
    data_info = data_info,
    mixtures = mixtures,
    analysis = analysis,
    fits = fits,
    prediction = prediction,
    bootstrap = bootstrap,
    results = results,
    labels = labels
  )

  validate_qgcompmulti(object)

  structure(object, class = "qgcompmulti")
}

#' @keywords internal
#' @noRd

validate_qgcompmulti <- function(object) {
  if (!is.list(object)) {
    stop("`object` must be a list.", call. = FALSE)
  }

  required_top <- qgcompmulti_required_top_components()
  actual_top <- names(object)

  if (is.null(actual_top) || !identical(actual_top, required_top)) {
    stop(
      sprintf(
        "`qgcompmulti` objects must contain top-level components in this order: %s.",
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

  component_fields <- qgcompmulti_component_fields()

  for (component_name in names(component_fields)) {
    qgcompmulti_validate_component(
      x = object[[component_name]],
      component_name = component_name,
      required_fields = component_fields[[component_name]]
    )
  }

  qgcompmulti_validate_analysis(
    object$analysis,
    allow_seed = TRUE,
    interval_context = "single_fit"
  )

  qgcompmulti_validate_labels(object$labels)
  qgcompmulti_validate_results(object$results, object$labels)
  qgcompmulti_validate_prediction(object$prediction)
  qgcompmulti_validate_bootstrap(object$bootstrap, object$labels)
  qgcompmulti_validate_fits(object$fits)

  invisible(object)
}
