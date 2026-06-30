test_that("quantized fits store intervention and MSM grids on the same scale", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_equal(
    fit$prediction$intervention_grid$psi1,
    fit$prediction$msm_grid$psi1
  )
  expect_equal(
    fit$prediction$intervention_grid$psi2,
    fit$prediction$msm_grid$psi2
  )
  expect_equal(
    fit$prediction$counterfactual_surface_target$exact_target,
    fit$prediction$counterfactual_surface$exact_mean
  )
  expect_equal(
    fit$prediction$msm_surface_target$msm_target,
    fit$prediction$msm_surface$msm_mean
  )
})
test_that("original-scale centered fits retain both intervention and MSM grids", {
  fit <- fit_test_model(interaction = TRUE, q = NULL, centering = "median")
  expect_false(identical(
    fit$prediction$intervention_grid$psi1,
    fit$prediction$msm_grid$psi1
  ))
  expect_false(identical(
    fit$prediction$intervention_grid$psi2,
    fit$prediction$msm_grid$psi2
  ))
  expect_equal(
    fit$prediction$surface_comparison$exact_mean,
    fit$prediction$counterfactual_surface$exact_mean
  )
  expect_equal(
    fit$prediction$surface_comparison$msm_mean,
    fit$prediction$msm_surface$msm_mean
  )
  expect_equal(
    fit$prediction$surface_comparison_target$exact_target,
    fit$prediction$counterfactual_surface_target$exact_target
  )
  expect_equal(
    fit$prediction$surface_comparison_target$msm_target,
    fit$prediction$msm_surface_target$msm_target
  )
})
test_that("bootstrap metadata remains coherent after fit-time surface retention", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_equal(
    fit$bootstrap$B_requested,
    fit$bootstrap$B_success + fit$bootstrap$B_failed
  )
  expect_equal(nrow(fit$bootstrap$coef_draws), fit$bootstrap$B_success)
})

test_that("stored prediction contract validates retained target-scale surfaces", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_true(is.data.frame(fit$prediction$counterfactual_surface_target))
  expect_true(is.data.frame(fit$prediction$msm_surface_target))
  expect_true(is.data.frame(fit$prediction$surface_comparison_target))
  bad <- fit
  bad$prediction$counterfactual_surface_target <- data.frame(grid_id = 1L)
  expect_error(
    validate_qgcompmulti(bad),
    "missing required columns"
  )
})
