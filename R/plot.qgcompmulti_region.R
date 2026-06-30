#' Plot a qgcomp.multi confidence region
#'
#' Draws the two-parameter confidence-region ellipsoid produced by
#' [confregion()].
#'
#' @param x A `qgcompmulti_region` object returned by [confregion()].
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param main Optional plot title.
#' @param col Line color for the ellipsoid boundary.
#' @param lwd Line width for the ellipsoid boundary.
#' @param pch Plotting character for the point estimate.
#' @param point_col Color for the point estimate.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#' @param asp Numeric aspect ratio passed to [graphics::plot()]. The default
#' `1` preserves the ellipsoid geometry on the fitting coefficient scale.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#'
#' @details
#' Direct plotting is supported only for two-parameter confidence regions. For
#' higher-dimensional regions, inspect the returned object fields directly.
#' Ratio-estimand regions are plotted on the log coefficient scale used for
#' fitting and covariance estimation.
#'
#' @export
plot.qgcompmulti_region <- function(x,
                                    xlab = NULL,
                                    ylab = NULL,
                                    main = NULL,
                                    col = "#2f6f9f",
                                    lwd = 2,
                                    pch = 19,
                                    point_col = "#1f2933",
                                    xlim = NULL,
                                    ylim = NULL,
                                    asp = 1,
                                    ...) {
  validate_qgcompmulti_region(x)
  if (length(x$parm) != 2L) {
    stop(
      "Direct plotting is supported only for two-parameter confidence regions.",
      call. = FALSE
    )
  }
  plot_data <- x$plot_data
  if (is.null(plot_data)) {
    plot_data <- qgcompmulti_region_plot_data(
      center = x$center,
      covariance = x$covariance,
      threshold = x$threshold
    )
  }
  if (is.null(xlab)) {
    xlab <- qgcompmulti_region_axis_label(
      unname(x$coefficient_labels[[1L]]),
      estimand_scale = x$estimand_scale,
      msm_fitting_scale = x$msm_fitting_scale
    )
  }
  if (is.null(ylab)) {
    ylab <- qgcompmulti_region_axis_label(
      unname(x$coefficient_labels[[2L]]),
      estimand_scale = x$estimand_scale,
      msm_fitting_scale = x$msm_fitting_scale
    )
  }
  if (is.null(main)) {
    main <- sprintf(
      "%s%% confidence region (%s scale)",
      formatC(100 * x$level, format = "f", digits = 1),
      qgcompmulti_msm_fitting_scale_label(x$msm_fitting_scale)
    )
  }
  if (is.null(xlim)) {
    xlim <- range(c(plot_data$x, x$center[[1L]]), finite = TRUE)
  }
  if (is.null(ylim)) {
    ylim <- range(c(plot_data$y, x$center[[2L]]), finite = TRUE)
  }
  graphics::plot(
    plot_data$x,
    plot_data$y,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    main = main,
    col = col,
    lwd = lwd,
    xlim = xlim,
    ylim = ylim,
    asp = asp,
    ...
  )
  graphics::abline(v = x$center[[1L]], h = x$center[[2L]], col = "gray80", lty = 3)
  graphics::lines(plot_data$x, plot_data$y, col = col, lwd = lwd)
  graphics::points(x$center[[1L]], x$center[[2L]], pch = pch, col = point_col)
  invisible(x)
}
