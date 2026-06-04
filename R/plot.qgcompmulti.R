#' Plot a qgcompmulti fit
#'
#' Produces a base-graphics display of the fitted marginal structural model
#' (MSM) surface from a [qgcomp.glm.multi()] fit. The default display is a
#' heatmap over the stored fit-time MSM grid. Contour plotting is also
#' supported, along with a slice-based interval display for MSM predictions.
#'
#' Plotting in `qgcomp.multi` is intentionally prediction-driven: the plotting
#' method delegates the scientific computation to [predict.qgcompmulti()] and
#' only handles rendering. This helps keep the plotted quantities aligned with
#' the documented public prediction interface.
#'
#' @param x A fitted `"qgcompmulti"` object.
#' @param style Character string specifying the surface display style.
#' Supported values are `"heatmap"` and `"contour"`.
#' @param grid Optional user-specified MSM grid with columns `psi1` and `psi2`.
#' If omitted, the stored fit-time MSM grid is used.
#' @param interval Logical; if `TRUE`, produces a slice-based line display with
#' bootstrap percentile intervals instead of a 2D surface plot.
#' @param slice Optional list with elements `var` and `value` describing the
#' fixed intervention coordinate for interval plotting. For example,
#' `list(var = "psi2", value = 1)` fixes `psi2` and varies `psi1` over the
#' stored MSM grid support.
#' @param level Confidence level used when `interval = TRUE`.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param main Optional plot title.
#' @param ... Additional graphical parameters passed to the underlying base
#' graphics call.
#'
#' @return Invisibly returns the structured prediction result used to draw the
#' plot.
#'
#' @details
#' The default `plot(fit)` call visualizes the fitted MSM surface over the
#' stored fit-time MSM grid. For quantized fits, the heatmap axes correspond to
#' the quantile-index intervention levels. For stored `q = NULL` heatmaps, the
#' plotted axis labels are mapped back to the pooled percentile intervention
#' values used at fit time, even when the MSM itself was fit on centered
#' coordinates.
#'
#' `interval = TRUE` switches to a slice-based line display rather than trying
#' to overlay uncertainty on the full two-dimensional surface. This keeps the
#' uncertainty display readable and aligned with the currently supported MSM
#' interval calculations.
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
#' fit <- qgcomp.glm.multi(
#'   f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
#'   data = dat,
#'   mix1 = c("X1", "X2", "X3"),
#'   mix2 = c("W1", "W2", "W3"),
#'   q = 4,
#'   B = 100,
#'   seed = 13
#' )
#'
#' # Default heatmap of the fitted MSM surface
#' plot(fit)
#'
#' # Contour rendering of the same stored MSM surface
#' plot(fit, style = "contour")
#'
#' # Slice-based interval display
#' plot(
#'   fit,
#'   interval = TRUE,
#'   slice = list(var = "psi2", value = 1)
#' )
#' }
#'
#' @export
plot.qgcompmulti <- function(x,
                             style = c("heatmap", "contour"),
                             grid = NULL,
                             interval = FALSE,
                             slice = NULL,
                             level = 0.95,
                             xlab = NULL,
                             ylab = NULL,
                             main = NULL,
                             ...) {
  validate_qgcompmulti(x)
  style <- match.arg(style)
  if (!is.logical(interval) || length(interval) != 1L || is.na(interval)) {
    stop("`interval` must be either `TRUE` or `FALSE`.", call. = FALSE)
  }
  default_labels <- qgcompmulti_default_surface_labels(x)
  if (isTRUE(interval)) {
    if (!is.null(grid)) {
      stop(
        "Interval slice plotting does not support `grid`; use `slice` to define the fixed coordinate.",
        call. = FALSE
      )
    }
    if (is.null(slice)) {
      stop(
        "Interval plotting requires `slice` to specify the fixed intervention coordinate.",
        call. = FALSE
      )
    }
    qgcompmulti_validate_conf_level(level)
    qgcompmulti_validate_plot_slice(slice)
    slice_grid <- qgcompmulti_build_slice_grid(x, slice)
    prediction_result <- predict(
      x,
      type = "msm",
      grid = slice_grid,
      interval = TRUE,
      level = level
    )
    varying_var <- if (slice$var == "psi1") "psi2" else "psi1"
    if (is.null(xlab)) {
      xlab <- paste0("MSM ", varying_var)
    }
    if (is.null(ylab)) {
      ylab <- "Predicted outcome"
    }
    if (is.null(main)) {
      main <- "MSM prediction with bootstrap interval"
    }
    qgcompmulti_plot_interval_slice(
      prediction_result = prediction_result,
      slice = slice,
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
    return(invisible(prediction_result))
  }
  prediction_result <- predict(
    x,
    type = "msm",
    grid = grid,
    interval = FALSE
  )
  plot_data <- qgcompmulti_surface_plot_data(prediction_result)
  if (is.null(xlab)) {
    xlab <- default_labels$xlab
  }
  if (is.null(ylab)) {
    ylab <- default_labels$ylab
  }
  if (is.null(main)) {
    main <- "Fitted MSM surface"
  }
  if (identical(style, "heatmap")) {
    qgcompmulti_plot_heatmap(
      plot_data = plot_data,
      object = x,
      prediction_result = prediction_result,
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
  } else {
    qgcompmulti_plot_contour(
      plot_data = plot_data,
      object = x,
      prediction_result = prediction_result,
      xlab = xlab,
      ylab = ylab,
      main = main,
      ...
    )
  }
  invisible(prediction_result)
}
