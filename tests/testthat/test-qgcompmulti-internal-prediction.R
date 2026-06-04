test_that("internal MSM grid builders support stored and user grids", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  stored <- qgcompmulti_build_msm_grid(fit)
  user <- qgcompmulti_build_msm_grid(
    fit,
    grid = data.frame(psi1 = c(0, 1.5), psi2 = c(0.5, 2), row.names = NULL)
  )
  expect_identical(stored$grid_type, "stored_fit_grid")
  expect_identical(user$grid_type, "user_grid")
  expect_identical(names(stored$grid), EXPECTED_GRID_COLUMNS)
  expect_identical(names(user$grid), EXPECTED_GRID_COLUMNS)
})
test_that("internal MSM surface predictions return structured results", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  result <- qgcompmulti_predict_msm_surface(fit)
  expect_identical(names(result), EXPECTED_INTERNAL_PREDICTION_FIELDS)
  expect_identical(result$prediction_type, "msm_surface")
  expect_identical(result$grid_type, "stored_fit_grid")
  expect_identical(result$grid_scale, "msm")
  expect_identical(result$estimand_scale, "response")
  expect_null(result$intervals)
  expect_identical(names(result$estimates), c("grid_id", "psi1", "psi2", "estimate"))
  expect_equal(nrow(result$estimates), nrow(fit$prediction$msm_grid))
})
test_that("internal MSM interval helpers use stored bootstrap draws", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  result <- qgcompmulti_predict_msm_surface(fit, level = 0.95)
  expect_identical(result$interval_type, "bootstrap_percentile")
  expect_identical(result$uncertainty_source, "stored_bootstrap_draws")
  expect_true(is.data.frame(result$intervals))
  expect_identical(
    names(result$intervals),
    c("grid_id", "psi1", "psi2", "lower", "upper")
  )
  expect_equal(nrow(result$intervals), nrow(result$estimates))
})
test_that("internal MSM point and contrast helpers work", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  point <- qgcompmulti_predict_msm_point(fit, psi1 = 1, psi2 = 2, level = 0.95)
  contrast <- qgcompmulti_predict_msm_contrast(
    fit,
    from = c(psi1 = 0, psi2 = 0),
    to = c(psi1 = 2, psi2 = 1),
    level = 0.95
  )
  expect_identical(point$prediction_type, "msm_point")
  expect_identical(point$grid_type, "point_regime")
  expect_equal(nrow(point$estimates), 1)
  expect_true(is.data.frame(point$intervals))
  expect_identical(contrast$prediction_type, "msm_contrast")
  expect_true(contrast$contrast)
  expect_equal(nrow(contrast$estimates), 1)
  expect_true(is.data.frame(contrast$intervals))
})
test_that("exact fit-grid surface extraction reuses stored fit-time surface", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  result <- qgcompmulti_exact_fit_surface(fit)
  expect_identical(result$prediction_type, "exact_fit_surface")
  expect_identical(result$grid_type, "stored_fit_grid")
  expect_identical(result$grid_scale, "intervention")
  expect_false(result$data_supplied)
  expect_equal(result$estimates, fit$prediction$counterfactual_surface)
})
test_that("exact arbitrary prediction requires explicit data and preserves scale distinction", {
  fit <- fit_test_model(interaction = TRUE, q = NULL, centering = "median")
  dat <- make_test_data()
  point_idx <- which(
    fit$prediction$msm_grid$psi1 == 0 &
      fit$prediction$msm_grid$psi2 == 0
  )[1]
  stored_at <- c(
    psi1 = fit$prediction$intervention_grid$psi1[point_idx],
    psi2 = fit$prediction$intervention_grid$psi2[point_idx]
  )
  expect_error(qgcompmulti_exact_predict_data(fit))
  result <- qgcompmulti_exact_predict_data(
    fit,
    data = dat,
    at = stored_at
  )
  expect_identical(result$prediction_type, "exact_arbitrary")
  expect_true(result$data_supplied)
  expect_equal(nrow(result$estimates), 1)
  expect_false(identical(
    result$estimates$intervention_psi1,
    result$estimates$msm_psi1
  ))
})
test_that("exact arbitrary contrast helper works on user-supplied data", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  dat <- make_test_data()
  contrast <- qgcompmulti_exact_contrast_data(
    fit,
    data = dat,
    from = c(psi1 = 0, psi2 = 0),
    to = c(psi1 = 1, psi2 = 2)
  )
  expect_identical(contrast$prediction_type, "exact_contrast")
  expect_true(contrast$data_supplied)
  expect_true(contrast$contrast)
  expect_equal(nrow(contrast$estimates), 1)
})
