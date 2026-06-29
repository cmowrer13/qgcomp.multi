# ------------------------------------------------------------------------------
# Estimand-scale helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_estimand_scales <- function() {
  c(
    "mean_difference",
    "risk_difference",
    "odds_ratio",
    "rate_ratio"
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_family_name_link <- function(family_name, link) {
  if (!is.character(family_name) ||
      length(family_name) != 1L ||
      is.na(family_name) ||
      !nzchar(family_name)) {
    stop("`family_name` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(link) ||
      length(link) != 1L ||
      is.na(link) ||
      !nzchar(link)) {
    stop("`link` must be a single non-empty character string.", call. = FALSE)
  }

  invisible(list(family_name = family_name, link = link))
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_estimand_scale_name <- function(estimand_scale,
                                                     allow_null = FALSE) {
  if (is.null(estimand_scale)) {
    if (isTRUE(allow_null)) {
      return(invisible(NULL))
    }
    stop("`estimand_scale` must not be `NULL`.", call. = FALSE)
  }

  allowed <- qgcompmulti_estimand_scales()

  if (!is.character(estimand_scale) ||
      length(estimand_scale) != 1L ||
      is.na(estimand_scale) ||
      !estimand_scale %in% allowed) {
    stop(
      sprintf(
        "`estimand_scale` must be one of: %s.",
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(estimand_scale)
}

#' @keywords internal
#' @noRd
qgcompmulti_allowed_estimand_scales <- function(family_name, link) {
  qgcompmulti_validate_family_name_link(family_name = family_name, link = link)

  if (identical(family_name, "gaussian")) {
    return("mean_difference")
  }

  if (identical(family_name, "binomial")) {
    if (identical(link, "logit")) {
      return(c("risk_difference", "odds_ratio"))
    }
    return("risk_difference")
  }

  if (identical(family_name, "poisson")) {
    if (identical(link, "log")) {
      return(c("mean_difference", "rate_ratio"))
    }
    return("mean_difference")
  }

  "mean_difference"
}

#' @keywords internal
#' @noRd
qgcompmulti_default_estimand_scale <- function(family_name,
                                               link,
                                               mode = c("planned", "current")) {
  qgcompmulti_validate_family_name_link(family_name = family_name, link = link)
  mode <- match.arg(mode)

  if (identical(mode, "current")) {
    if (identical(family_name, "binomial")) {
      return("risk_difference")
    }
    return("mean_difference")
  }

  if (identical(family_name, "binomial")) {
    if (identical(link, "logit")) {
      return("odds_ratio")
    }
    return("risk_difference")
  }

  if (identical(family_name, "poisson")) {
    if (identical(link, "log")) {
      return("rate_ratio")
    }
    return("mean_difference")
  }

  "mean_difference"
}

#' @keywords internal
#' @noRd
qgcompmulti_validate_estimand_scale <- function(family_name,
                                                link,
                                                estimand_scale,
                                                allow_null = FALSE) {
  qgcompmulti_validate_family_name_link(family_name = family_name, link = link)
  qgcompmulti_validate_estimand_scale_name(
    estimand_scale = estimand_scale,
    allow_null = allow_null
  )

  if (is.null(estimand_scale)) {
    return(invisible(NULL))
  }

  allowed <- qgcompmulti_allowed_estimand_scales(
    family_name = family_name,
    link = link
  )

  if (!estimand_scale %in% allowed) {
    stop(
      sprintf(
        "`estimand_scale = \"%s\"` is not supported for `%s(link = \"%s\")`. Allowed values are: %s.",
        estimand_scale,
        family_name,
        link,
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(estimand_scale)
}

#' @keywords internal
#' @noRd
qgcompmulti_is_ratio_estimand <- function(estimand_scale) {
  qgcompmulti_validate_estimand_scale_name(estimand_scale)
  estimand_scale %in% c("odds_ratio", "rate_ratio")
}

#' @keywords internal
#' @noRd
qgcompmulti_msm_fitting_scale <- function(estimand_scale) {
  qgcompmulti_validate_estimand_scale_name(estimand_scale)

  switch(
    estimand_scale,
    mean_difference = "identity",
    risk_difference = "identity",
    odds_ratio = "logit",
    rate_ratio = "log"
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_estimand_label <- function(estimand_scale) {
  qgcompmulti_validate_estimand_scale_name(estimand_scale)

  switch(
    estimand_scale,
    mean_difference = "Mean difference",
    risk_difference = "Risk difference",
    odds_ratio = "Odds ratio",
    rate_ratio = "Rate ratio"
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_msm_fitting_scale_label <- function(msm_fitting_scale) {
  if (!is.character(msm_fitting_scale) ||
      length(msm_fitting_scale) != 1L ||
      is.na(msm_fitting_scale) ||
      !msm_fitting_scale %in% c("identity", "logit", "log")) {
    stop(
      "`msm_fitting_scale` must be one of: identity, logit, log.",
      call. = FALSE
    )
  }

  switch(
    msm_fitting_scale,
    identity = "identity",
    logit = "logit",
    log = "log"
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_transform_msm_surface <- function(values,
                                              msm_fitting_scale,
                                              direction = c("to_fitting", "to_response")) {
  direction <- match.arg(direction)
  qgcompmulti_msm_fitting_scale_label(msm_fitting_scale)

  if (!is.numeric(values)) {
    stop("`values` must be numeric.", call. = FALSE)
  }

  transform_fun <- switch(
    paste(msm_fitting_scale, direction, sep = "::"),
    "identity::to_fitting" = identity,
    "identity::to_response" = identity,
    "logit::to_fitting" = qlogis,
    "logit::to_response" = plogis,
    "log::to_fitting" = log,
    "log::to_response" = exp
  )

  transform_fun(values)
}

#' @keywords internal
#' @noRd
qgcompmulti_transform_msm_coefficients <- function(values,
                                                   estimand_scale,
                                                   direction = c("to_fitting", "to_display")) {
  direction <- match.arg(direction)
  qgcompmulti_validate_estimand_scale_name(estimand_scale)

  if (!is.numeric(values)) {
    stop("`values` must be numeric.", call. = FALSE)
  }

  if (!qgcompmulti_is_ratio_estimand(estimand_scale)) {
    return(values)
  }

  if (identical(direction, "to_fitting")) {
    return(log(values))
  }

  exp(values)
}

#' @keywords internal
#' @noRd
qgcompmulti_resolve_estimand_spec <- function(family,
                                              estimand_scale = NULL,
                                              mode = c("planned", "current")) {
  if (!inherits(family, "family")) {
    stop("`family` must inherit from class `family`.", call. = FALSE)
  }

  mode <- match.arg(mode)
  family_name <- family$family
  link <- family$link

  if (is.null(estimand_scale)) {
    chosen <- qgcompmulti_default_estimand_scale(
      family_name = family_name,
      link = link,
      mode = mode
    )
    defaulted <- TRUE
  } else {
    qgcompmulti_validate_estimand_scale(
      family_name = family_name,
      link = link,
      estimand_scale = estimand_scale
    )
    chosen <- estimand_scale
    defaulted <- FALSE
  }

  list(
    estimand_scale = chosen,
    estimand_scale_defaulted = defaulted,
    msm_fitting_scale = qgcompmulti_msm_fitting_scale(chosen)
  )
}
