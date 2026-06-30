# ------------------------------------------------------------------------------
# Plot and presentation helpers
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
qgcompmulti_labeled_coef_table <- function(results, labels, analysis = NULL) {
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

  out <- if (is.null(analysis) ||
             !qgcompmulti_is_ratio_estimand(analysis$estimand_scale)) {
    coef_table
  } else {
    display_estimate <- qgcompmulti_transform_msm_coefficients(
      coef_table$Estimate,
      estimand_scale = analysis$estimand_scale,
      direction = "to_display"
    )
    display_name <- qgcompmulti_estimand_label(analysis$estimand_scale)
    out <- data.frame(
      display_estimate,
      Coefficient = coef_table$Estimate,
      `Std. Error` = coef_table$`Std. Error`,
      check.names = FALSE
    )
    names(out)[1L] <- display_name
    statistic_col <- grep(" value$", names(coef_table), value = TRUE)
    if (length(statistic_col) == 1L) {
      out[[statistic_col]] <- coef_table[[statistic_col]]
    }
    if ("df" %in% names(coef_table)) {
      out$df <- coef_table$df
    }
    p_col <- grep("^Pr\\(>", names(coef_table), value = TRUE)
    if (length(p_col) == 1L) {
      out[[p_col]] <- coef_table[[p_col]]
    }
    row.names(out) <- rownames(coef_table)
    out
  }

  rownames(out) <- unname(coef_labels)
  out
}

#' @keywords internal
#' @noRd
qgcompmulti_scale_description <- function(analysis) {
  estimand <- qgcompmulti_estimand_label(analysis$estimand_scale)
  fitting_scale <- qgcompmulti_msm_fitting_scale_label(analysis$msm_fitting_scale)
  if (isTRUE(analysis$estimand_scale_defaulted)) {
    estimand <- paste0(estimand, " (default)")
  }
  list(
    estimand = estimand,
    fitting_scale = fitting_scale,
    default_interval_method = analysis$default_interval_method
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_print_scale_info <- function(analysis, mi = FALSE) {
  scale_info <- qgcompmulti_scale_description(analysis)
  cat("  Estimand: ", scale_info$estimand, "
", sep = "")
  cat("  MSM fitting scale: ", scale_info$fitting_scale, "
", sep = "")
  if (isTRUE(mi)) {
    cat("  Interval method: wald (pooled multiple imputation)
", sep = "")
  } else {
    cat("  Default interval method: ", scale_info$default_interval_method, "
", sep = "")
  }
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

#' @keywords internal
#' @noRd
qgcompmulti_default_surface_labels <- function(object) {
  validate_qgcompmulti(object)
  if (!is.null(object$mixtures$q)) {
    return(list(
      xlab = "Mixture 1 intervention (quantile index)",
      ylab = "Mixture 2 intervention (quantile index)"
    ))
  }
  list(
    xlab = "Mixture 1 intervention",
    ylab = "Mixture 2 intervention"
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_format_axis_values <- function(x, digits = 3) {
  formatC(x, digits = digits, format = "fg", flag = "#")
}
#' @keywords internal
#' @noRd
qgcompmulti_surface_plot_data <- function(prediction_result) {
  if (!is.list(prediction_result) || is.null(prediction_result$estimates)) {
    stop("`prediction_result` must be a structured prediction result.", call. = FALSE)
  }
  estimates <- prediction_result$estimates
  required_cols <- c("psi1", "psi2", "estimate")
  if (!all(required_cols %in% names(estimates))) {
    stop(
      sprintf(
        "Prediction estimates must contain columns: %s.",
        paste(required_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (anyDuplicated(estimates[c("psi1", "psi2")]) > 0L) {
    stop(
      "Surface plotting requires at most one prediction per (`psi1`, `psi2`) pair.",
      call. = FALSE
    )
  }
  x_vals <- sort(unique(estimates$psi1))
  y_vals <- sort(unique(estimates$psi2))
  expected_n <- length(x_vals) * length(y_vals)
  if (nrow(estimates) != expected_n) {
    stop(
      "Surface plotting requires a rectangular prediction grid covering all (`psi1`, `psi2`) combinations.",
      call. = FALSE
    )
  }
  z_mat <- matrix(
    NA_real_,
    nrow = length(x_vals),
    ncol = length(y_vals),
    dimnames = list(as.character(x_vals), as.character(y_vals))
  )
  row_idx <- match(estimates$psi1, x_vals)
  col_idx <- match(estimates$psi2, y_vals)
  z_mat[cbind(row_idx, col_idx)] <- estimates$estimate
  if (anyNA(z_mat)) {
    stop(
      "Surface plotting requires a complete rectangular prediction grid.",
      call. = FALSE
    )
  }
  list(
    x = x_vals,
    y = y_vals,
    z = z_mat,
    estimates = estimates
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_surface_axis_spec <- function(object,
                                          prediction_result,
                                          plot_data) {
  x_at <- plot_data$x
  y_at <- plot_data$y
  x_labels <- qgcompmulti_format_axis_values(plot_data$x)
  y_labels <- qgcompmulti_format_axis_values(plot_data$y)
  if (!is.null(object$mixtures$q) ||
      !identical(prediction_result$grid_type, "stored_fit_grid")) {
    return(list(
      x_at = x_at,
      x_labels = x_labels,
      y_at = y_at,
      y_labels = y_labels
    ))
  }
  intervention_grid <- object$prediction$intervention_grid
  msm_grid <- object$prediction$msm_grid
  if (is.null(intervention_grid) || is.null(msm_grid)) {
    return(list(
      x_at = x_at,
      x_labels = x_labels,
      y_at = y_at,
      y_labels = y_labels
    ))
  }
  x_map <- unique(data.frame(
    msm = msm_grid$psi1,
    intervention = intervention_grid$psi1,
    row.names = NULL
  ))
  y_map <- unique(data.frame(
    msm = msm_grid$psi2,
    intervention = intervention_grid$psi2,
    row.names = NULL
  ))
  x_map <- x_map[order(x_map$msm), , drop = FALSE]
  y_map <- y_map[order(y_map$msm), , drop = FALSE]
  if (length(plot_data$x) == nrow(x_map) &&
      length(plot_data$y) == nrow(y_map) &&
      isTRUE(all.equal(plot_data$x, x_map$msm, tolerance = 1e-8)) &&
      isTRUE(all.equal(plot_data$y, y_map$msm, tolerance = 1e-8))) {
    x_labels <- qgcompmulti_format_axis_values(x_map$intervention)
    y_labels <- qgcompmulti_format_axis_values(y_map$intervention)
  }
  list(
    x_at = x_at,
    x_labels = x_labels,
    y_at = y_at,
    y_labels = y_labels
  )
}
#' @keywords internal
#' @noRd
qgcompmulti_plot_heatmap <- function(plot_data,
                                     object,
                                     prediction_result,
                                     xlab,
                                     ylab,
                                     main,
                                     ...) {
  axis_spec <- qgcompmulti_surface_axis_spec(object, prediction_result, plot_data)
  dots <- list(...)
  cols <- dots$col
  if (is.null(cols)) {
    cols <- grDevices::hcl.colors(64L, palette = "YlOrRd", rev = TRUE)
  }
  dots$col <- NULL
  zlim <- range(plot_data$z, finite = TRUE)
  old_mar <- graphics::par("mar")
  graphics::layout(matrix(c(1L, 2L), nrow = 1L), widths = c(4, 1))
  on.exit({
    graphics::layout(1)
    graphics::par(mar = old_mar)
  }, add = TRUE)
  graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  do.call(
    graphics::image,
    c(
      list(
        x = plot_data$x,
        y = plot_data$y,
        z = plot_data$z,
        xlab = xlab,
        ylab = ylab,
        main = main,
        col = cols,
        axes = FALSE
      ),
      dots
    )
  )
  graphics::axis(1, at = axis_spec$x_at, labels = axis_spec$x_labels)
  graphics::axis(2, at = axis_spec$y_at, labels = axis_spec$y_labels)
  graphics::box()
  legend_breaks <- seq(zlim[1], zlim[2], length.out = length(cols) + 1L)
  graphics::par(mar = c(5.1, 1.5, 4.1, 4.6))
  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0, 1),
    ylim = zlim,
    xaxs = "i",
    yaxs = "i"
  )
  graphics::rect(
    xleft = 0,
    ybottom = legend_breaks[-length(legend_breaks)],
    xright = 1,
    ytop = legend_breaks[-1L],
    col = cols,
    border = NA
  )
  graphics::axis(4, at = pretty(zlim), labels = qgcompmulti_format_axis_values(pretty(zlim)))
  graphics::mtext("Predicted outcome", side = 4, line = 2.8)
  graphics::box()
  invisible(plot_data)
}
#' @keywords internal
#' @noRd
qgcompmulti_plot_contour <- function(plot_data,
                                     object,
                                     prediction_result,
                                     xlab,
                                     ylab,
                                     main,
                                     ...) {
  axis_spec <- qgcompmulti_surface_axis_spec(object, prediction_result, plot_data)
  graphics::contour(
    x = plot_data$x,
    y = plot_data$y,
    z = plot_data$z,
    xlab = xlab,
    ylab = ylab,
    main = main,
    axes = FALSE,
    ...
  )
  graphics::axis(1, at = axis_spec$x_at, labels = axis_spec$x_labels)
  graphics::axis(2, at = axis_spec$y_at, labels = axis_spec$y_labels)
  graphics::box()
  invisible(plot_data)
}
#' @keywords internal
#' @noRd
qgcompmulti_validate_plot_slice <- function(slice) {
  if (!is.list(slice) || is.data.frame(slice)) {
    stop(
      "`slice` must be a list with elements `var` and `value`.",
      call. = FALSE
    )
  }
  if (!all(c("var", "value") %in% names(slice))) {
    stop("`slice` must contain `var` and `value`.", call. = FALSE)
  }
  if (!is.character(slice$var) ||
      length(slice$var) != 1L ||
      is.na(slice$var) ||
      !slice$var %in% c("psi1", "psi2")) {
    stop("`slice$var` must be either \"psi1\" or \"psi2\".", call. = FALSE)
  }
  if (!is.numeric(slice$value) ||
      length(slice$value) != 1L ||
      is.na(slice$value) ||
      !is.finite(slice$value)) {
    stop("`slice$value` must be a single finite numeric value.", call. = FALSE)
  }
  invisible(slice)
}
#' @keywords internal
#' @noRd
qgcompmulti_build_slice_grid <- function(object, slice) {
  validate_qgcompmulti(object)
  qgcompmulti_validate_plot_slice(slice)
  support_grid <- object$prediction$msm_grid
  if (is.null(support_grid)) {
    stop("No stored MSM grid is available for slice plotting.", call. = FALSE)
  }
  if (slice$var == "psi1") {
    psi2_vals <- sort(unique(support_grid$psi2))
    grid <- data.frame(
      psi1 = rep(slice$value, length(psi2_vals)),
      psi2 = psi2_vals,
      row.names = NULL
    )
  } else {
    psi1_vals <- sort(unique(support_grid$psi1))
    grid <- data.frame(
      psi1 = psi1_vals,
      psi2 = rep(slice$value, length(psi1_vals)),
      row.names = NULL
    )
  }
  qgcompmulti_validate_grid_support(object, grid, scale = "msm")
  qgcompmulti_standardize_grid(grid)
}
#' @keywords internal
#' @noRd
qgcompmulti_plot_interval_slice <- function(prediction_result,
                                            slice,
                                            xlab,
                                            ylab,
                                            main,
                                            ...) {
  qgcompmulti_validate_plot_slice(slice)
  estimates <- prediction_result$estimates
  intervals <- prediction_result$intervals
  if (is.null(intervals)) {
    stop(
      "Interval slice plotting requires prediction intervals.",
      call. = FALSE
    )
  }
  varying_var <- if (slice$var == "psi1") "psi2" else "psi1"
  x_vals <- estimates[[varying_var]]
  ord <- order(x_vals)
  x_vals <- x_vals[ord]
  est_vals <- estimates$estimate[ord]
  lower_vals <- intervals$lower[ord]
  upper_vals <- intervals$upper[ord]
  graphics::plot(
    x = x_vals,
    y = est_vals,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    main = main,
    ylim = range(c(lower_vals, upper_vals), na.rm = TRUE),
    ...
  )
  graphics::lines(x_vals, lower_vals, lty = 2)
  graphics::lines(x_vals, upper_vals, lty = 2)
  invisible(prediction_result)
}
