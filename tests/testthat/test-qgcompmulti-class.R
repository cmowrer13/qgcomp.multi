test_that("qgcomp.glm.multi returns a qgcompmulti object with expected structure", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_s3_class(fit, "qgcompmulti")
  expect_true(is.list(fit))
  expect_identical(names(fit), EXPECTED_TOP_COMPONENTS)
  for (component_name in names(EXPECTED_COMPONENT_FIELDS)) {
    expect_true(is.list(fit[[component_name]]))
    expect_identical(
      names(fit[[component_name]]),
      EXPECTED_COMPONENT_FIELDS[[component_name]]
    )
  }
})

test_that("fitted object stores expected result, fit, and label components", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_true(is.numeric(fit$results$coefficients))
  expect_true(is.numeric(fit$results$std_error))
  expect_true(is.matrix(fit$results$vcov))
  expect_true(is.data.frame(fit$results$coef_table))
  expect_s3_class(fit$fits$outcome_fit, "glm")
  expect_s3_class(fit$fits$msm_fit, "glm")
  expect_true(is.data.frame(fit$prediction$intervention_grid))
  expect_true(is.data.frame(fit$prediction$msm_grid))
  expect_true(is.data.frame(fit$prediction$counterfactual_surface))
  expect_true(is.data.frame(fit$prediction$msm_surface))
  expect_true(is.data.frame(fit$prediction$surface_comparison))
  expect_true(is.data.frame(fit$prediction$counterfactual_surface_target))
  expect_true(is.data.frame(fit$prediction$msm_surface_target))
  expect_true(is.data.frame(fit$prediction$surface_comparison_target))
  expect_true(is.matrix(fit$bootstrap$coef_draws))
  expect_equal(fit$bootstrap$B_requested, 10)
  expect_equal(fit$bootstrap$B_success, 10)
  expect_equal(fit$bootstrap$B_failed, 0)
  expect_null(fit$bootstrap$failure_log)
  expect_true(is.character(fit$labels$coefficient_names))
  expect_true(is.character(fit$labels$coefficient_labels))
  expect_true(is.character(fit$labels$mixture_labels))
  expect_identical(names(fit$labels$mixture_labels), c("mix1", "mix2"))
  expect_identical(fit$analysis$estimand_scale, "mean_difference")
  expect_true(fit$analysis$estimand_scale_defaulted)
  expect_identical(fit$analysis$msm_fitting_scale, "identity")
  expect_identical(fit$analysis$default_interval_method, "wald")
})

test_that("stored fit-time grids and surfaces align on a common grid", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_identical(names(fit$prediction$intervention_grid), EXPECTED_GRID_COLUMNS)
  expect_identical(names(fit$prediction$msm_grid), EXPECTED_GRID_COLUMNS)
  expect_identical(
    names(fit$prediction$counterfactual_surface),
    EXPECTED_COUNTERFACTUAL_SURFACE_COLUMNS
  )
  expect_identical(
    names(fit$prediction$msm_surface),
    EXPECTED_MSM_SURFACE_COLUMNS
  )
  expect_identical(
    names(fit$prediction$surface_comparison),
    EXPECTED_SURFACE_COMPARISON_COLUMNS
  )
  expect_identical(
    names(fit$prediction$counterfactual_surface_target),
    EXPECTED_COUNTERFACTUAL_TARGET_COLUMNS
  )
  expect_identical(
    names(fit$prediction$msm_surface_target),
    EXPECTED_MSM_TARGET_COLUMNS
  )
  expect_identical(
    names(fit$prediction$surface_comparison_target),
    EXPECTED_SURFACE_COMPARISON_TARGET_COLUMNS
  )
  expect_identical(
    fit$prediction$intervention_grid$grid_id,
    fit$prediction$msm_grid$grid_id
  )
  expect_identical(
    fit$prediction$counterfactual_surface$grid_id,
    fit$prediction$msm_surface$grid_id
  )
  expect_identical(
    fit$prediction$counterfactual_surface$grid_id,
    fit$prediction$surface_comparison$grid_id
  )
  expect_identical(
    fit$prediction$counterfactual_surface_target$grid_id,
    fit$prediction$msm_surface_target$grid_id
  )
  expect_identical(
    fit$prediction$counterfactual_surface_target$grid_id,
    fit$prediction$surface_comparison_target$grid_id
  )
})
