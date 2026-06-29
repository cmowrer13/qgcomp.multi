test_that("qgcompmulti_msm_fit returns the expected component structure", {
  engine <- fit_engine_result(interaction = TRUE, q = 4)
  expect_true(is.list(engine))
  expect_identical(
    names(engine),
    c(
      "outcome_fit",
      "msm_fit",
      "coefficients",
      "n_used",
      "intervention_grid",
      "msm_grid",
      "counterfactual_surface",
      "msm_surface",
      "surface_comparison",
      "counterfactual_surface_target",
      "msm_surface_target",
      "surface_comparison_target"
    )
  )
  expect_s3_class(engine$outcome_fit, "glm")
  expect_s3_class(engine$msm_fit, "glm")
})
test_that("qgcompmulti_msm_fit returns correctly named coefficients by branch", {
  engine_interaction <- fit_engine_result(interaction = TRUE, q = 4)
  engine_no_interaction <- fit_engine_result(interaction = FALSE, q = 4)
  expect_identical(
    names(engine_interaction$coefficients),
    EXPECTED_COEF_NAMES_WITH_INTERACTION
  )
  expect_identical(
    names(engine_no_interaction$coefficients),
    EXPECTED_COEF_NAMES_NO_INTERACTION
  )
})

test_that("qgcompmulti_msm_fit works for original-scale q = NULL fits", {
  engine <- fit_engine_result(
    interaction = TRUE,
    q = NULL,
    centering = "median"
  )
  expect_true(is.list(engine))
  expect_identical(names(engine$coefficients), EXPECTED_COEF_NAMES_WITH_INTERACTION)
  expect_equal(engine$n_used, nrow(make_test_data()))
})

test_that("qgcompmulti_msm_fit retains aligned fit-time grids and surfaces", {
  engine <- fit_engine_result(interaction = TRUE, q = 4)
  expect_identical(names(engine$intervention_grid), EXPECTED_GRID_COLUMNS)
  expect_identical(names(engine$msm_grid), EXPECTED_GRID_COLUMNS)
  expect_identical(
    names(engine$counterfactual_surface),
    EXPECTED_COUNTERFACTUAL_SURFACE_COLUMNS
  )
  expect_identical(
    names(engine$msm_surface),
    EXPECTED_MSM_SURFACE_COLUMNS
  )
  expect_identical(
    names(engine$surface_comparison),
    EXPECTED_SURFACE_COMPARISON_COLUMNS
  )
  expect_identical(
    names(engine$counterfactual_surface_target),
    EXPECTED_COUNTERFACTUAL_TARGET_COLUMNS
  )
  expect_identical(
    names(engine$msm_surface_target),
    EXPECTED_MSM_TARGET_COLUMNS
  )
  expect_identical(
    names(engine$surface_comparison_target),
    EXPECTED_SURFACE_COMPARISON_TARGET_COLUMNS
  )
  expect_equal(nrow(engine$intervention_grid), nrow(engine$msm_grid))
  expect_equal(nrow(engine$counterfactual_surface), nrow(engine$msm_surface))
  expect_equal(nrow(engine$counterfactual_surface), nrow(engine$surface_comparison))
  expect_equal(
    nrow(engine$counterfactual_surface_target),
    nrow(engine$msm_surface_target)
  )
  expect_equal(
    nrow(engine$counterfactual_surface_target),
    nrow(engine$surface_comparison_target)
  )
  expect_identical(
    engine$counterfactual_surface$grid_id,
    engine$surface_comparison$grid_id
  )
})

test_that("qgcompmulti_msm_fit stores target-scale surfaces for ratio estimands", {
  binomial_engine <- fit_engine_result(
    interaction = TRUE,
    q = 4,
    family = binomial(link = "logit")
  )
  poisson_engine <- fit_engine_result(
    interaction = TRUE,
    q = 4,
    family = poisson(link = "log")
  )
  expect_false(identical(
    binomial_engine$counterfactual_surface$exact_mean,
    binomial_engine$counterfactual_surface_target$exact_target
  ))
  expect_false(identical(
    poisson_engine$counterfactual_surface$exact_mean,
    poisson_engine$counterfactual_surface_target$exact_target
  ))
})
