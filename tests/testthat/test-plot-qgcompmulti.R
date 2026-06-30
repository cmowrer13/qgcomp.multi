with_plot_device <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  force(code)
}
test_that("plot.qgcompmulti() draws the default MSM heatmap", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  expect_no_error(
    with_plot_device(plot(fit))
  )
})
test_that("plot.qgcompmulti() supports repeated heatmaps on the same device", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  expect_no_error(
    with_plot_device({
      plot(fit)
      plot(fit)
      plot(fit)
    })
  )
})
test_that("plot.qgcompmulti() draws q = NULL heatmaps without error", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 10
  )
  expect_no_error(
    with_plot_device(plot(fit))
  )
})
test_that("plot.qgcompmulti() supports contour output", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  expect_no_error(
    with_plot_device(plot(fit, style = "contour"))
  )
})
test_that("plot.qgcompmulti() keeps transformed-scale fits on the response scale", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    B = 10,
    family = binomial(link = "logit")
  )
  result <- with_plot_device(plot(fit, scale = "response"))
  expect_identical(result$estimate_scale, "response")
  expect_identical(result$estimand_scale, "response")
  expect_identical(result$fit_estimand_scale, "odds_ratio")
})
test_that("plot.qgcompmulti() supports an in-support user MSM grid", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  grid <- expand.grid(
    psi1 = c(0, 1.5, 3),
    psi2 = c(0, 1.5, 3),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  expect_no_error(
    with_plot_device(plot(fit, grid = grid))
  )
})
test_that("plot.qgcompmulti() supports slice-based interval plotting", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  expect_no_error(
    with_plot_device(
      plot(
        fit,
        interval = TRUE,
        slice = list(var = "psi2", value = 1),
        level = 0.95
      )
    )
  )
})
test_that("plot.qgcompmulti() rejects invalid plotting inputs clearly", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  non_rectangular_grid <- data.frame(
    psi1 = c(0, 1, 2),
    psi2 = c(0, 1, 0),
    row.names = NULL
  )
  expect_error(
    with_plot_device(plot(fit, interval = TRUE)),
    "requires `slice`"
  )
  expect_error(
    with_plot_device(plot(fit, interval = TRUE, slice = list(var = "psi3", value = 1))),
    "`slice\\$var` must be either"
  )
  expect_error(
    with_plot_device(plot(fit, interval = TRUE, slice = list(var = "psi2", value = 1), grid = data.frame(psi1 = 0, psi2 = 0))),
    "does not support `grid`"
  )
  expect_error(
    with_plot_device(plot(fit, grid = non_rectangular_grid)),
    "rectangular prediction grid"
  )
  expect_error(
    with_plot_device(plot(fit, scale = "estimand")),
    "supports only `scale = \"response\"`"
  )
})
test_that("stored q = NULL surface plots use intervention-scale axis labels", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 10
  )
  prediction_result <- predict(fit, type = "msm")
  plot_data <- qgcompmulti_surface_plot_data(prediction_result)
  axis_spec <- qgcompmulti_surface_axis_spec(fit, prediction_result, plot_data)
  expected_x <- qgcompmulti_format_axis_values(
    sort(unique(fit$prediction$intervention_grid$psi1))
  )
  expected_y <- qgcompmulti_format_axis_values(
    sort(unique(fit$prediction$intervention_grid$psi2))
  )
  expect_identical(axis_spec$x_labels, expected_x)
  expect_identical(axis_spec$y_labels, expected_y)
  expect_false(
    identical(
      axis_spec$x_labels,
      qgcompmulti_format_axis_values(sort(unique(fit$prediction$msm_grid$psi1)))
    )
  )
})
test_that("default plot labels identify response-scale outcome quantities", {
  gaussian_fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  binomial_fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    B = 10,
    family = binomial(link = "logit")
  )
  poisson_fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    B = 10,
    family = poisson(link = "log")
  )
  expect_identical(
    qgcompmulti_default_surface_labels(gaussian_fit)$outcome_label,
    "Predicted mean outcome (response scale)"
  )
  expect_identical(
    qgcompmulti_default_surface_labels(binomial_fit)$outcome_label,
    "Predicted risk (response scale)"
  )
  expect_identical(
    qgcompmulti_default_surface_labels(poisson_fit)$outcome_label,
    "Predicted expected count (response scale)"
  )
})
