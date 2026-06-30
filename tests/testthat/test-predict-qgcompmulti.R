test_that("predict defaults to MSM surface predictions", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  result <- predict(fit)
  expect_identical(names(result), EXPECTED_INTERNAL_PREDICTION_FIELDS)
  expect_identical(result$prediction_type, "msm_surface")
  expect_identical(result$grid_type, "stored_fit_grid")
  expect_equal(nrow(result$estimates), nrow(fit$prediction$msm_grid))
})
test_that("predict supports public MSM and exact workflows for q = NULL fits", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 10
  )
  dat <- make_test_data()
  msm_result <- predict(fit)
  exact_result <- predict(fit, type = "exact")
  exact_arbitrary <- predict(
    fit,
    type = "exact",
    data = dat,
    at = c(psi1 = fit$prediction$intervention_grid$psi1[1],
           psi2 = fit$prediction$intervention_grid$psi2[1])
  )
  expect_identical(msm_result$prediction_type, "msm_surface")
  expect_identical(exact_result$prediction_type, "exact_fit_surface")
  expect_identical(exact_arbitrary$prediction_type, "exact_arbitrary")
})
test_that("predict supports user-grid MSM surface predictions", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  result <- predict(
    fit,
    type = "msm",
    grid = data.frame(psi1 = c(0, 1), psi2 = c(2, 3), row.names = NULL)
  )
  expect_identical(result$prediction_type, "msm_surface")
  expect_identical(result$grid_type, "user_grid")
  expect_equal(nrow(result$estimates), 2)
})
test_that("predict allows interpolation within support but rejects extrapolation", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  interpolated <- predict(
    fit,
    type = "msm",
    grid = data.frame(psi1 = c(1.5, 2.5), psi2 = c(0.25, 2.75), row.names = NULL)
  )
  expect_identical(interpolated$grid_type, "user_grid")
  expect_equal(nrow(interpolated$estimates), 2)
  expect_error(
    predict(
      fit,
      type = "msm",
      grid = data.frame(psi1 = c(-1, 1), psi2 = c(0, 2), row.names = NULL)
    )
  )
  expect_error(
    predict(
      fit,
      type = "msm_point",
      at = c(psi1 = 5, psi2 = 1)
    )
  )
})
test_that("predict supports MSM point predictions and intervals", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  result <- predict(
    fit,
    type = "msm_point",
    at = c(psi1 = 1, psi2 = 2),
    interval = TRUE
  )
  expect_identical(result$prediction_type, "msm_point")
  expect_identical(result$grid_type, "point_regime")
  expect_true(is.data.frame(result$intervals))
  expect_equal(nrow(result$estimates), 1)
})
test_that("predict supports MSM contrasts and intervals", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  result <- predict(
    fit,
    type = "msm_contrast",
    from = c(psi1 = 0, psi2 = 0),
    to = c(psi1 = 2, psi2 = 1),
    interval = TRUE
  )
  expect_identical(result$prediction_type, "msm_contrast")
  expect_true(result$contrast)
  expect_true(is.data.frame(result$intervals))
  expect_equal(nrow(result$estimates), 1)
})
test_that("predict supports no-interaction fits through the public API", {
  fit <- fit_test_model(interaction = FALSE, q = 4, B = 10)
  result <- predict(fit)
  point <- predict(fit, type = "msm_point", at = c(psi1 = 1, psi2 = 2))
  expect_identical(result$prediction_type, "msm_surface")
  expect_equal(nrow(result$estimates), nrow(fit$prediction$msm_grid))
  expect_identical(point$prediction_type, "msm_point")
})
test_that("predict exposes exact fit-time surface through the public API", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  result <- predict(fit, type = "exact")
  expect_identical(result$prediction_type, "exact_fit_surface")
  expect_equal(result$estimates, fit$prediction$counterfactual_surface)
})
test_that("predict supports exact arbitrary prediction and exact arbitrary contrasts", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  dat <- make_test_data()
  exact_pred <- predict(
    fit,
    type = "exact",
    data = dat,
    at = c(psi1 = 1, psi2 = 2)
  )
  exact_contrast <- predict(
    fit,
    type = "exact_contrast",
    data = dat,
    from = c(psi1 = 0, psi2 = 0),
    to = c(psi1 = 1, psi2 = 2)
  )
  expect_identical(exact_pred$prediction_type, "exact_arbitrary")
  expect_true(exact_pred$data_supplied)
  expect_identical(exact_contrast$prediction_type, "exact_contrast")
  expect_true(exact_contrast$data_supplied)
  expect_true(exact_contrast$contrast)
})
test_that("predict fails clearly on invalid public combinations", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  dat <- make_test_data()
  expect_error(predict(fit, type = "msm_point"))
  expect_error(predict(fit, type = "msm_contrast", from = c(psi1 = 0, psi2 = 0)))
  expect_error(predict(fit, type = "exact", at = c(psi1 = 1, psi2 = 2)))
  expect_error(predict(fit, type = "exact", data = dat, interval = TRUE),
               "Version 0.4.0"
               )
  expect_error(predict(fit, type = "exact_contrast", from = c(psi1 = 0, psi2 = 0), to = c(psi1 = 1, psi2 = 2)))
  expect_error(predict(fit, type = "msm", at = c(psi1 = 1, psi2 = 2)))
  expect_error(
    predict(
      fit,
      type = "exact",
      data = dat,
      at = c(psi1 = 99, psi2 = 99)
    )
  )
})

test_that("predict keeps transformed-fit MSM surfaces and points on the response scale", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    family = binomial(link = "logit")
  )
  surface <- predict(fit, type = "msm", interval = TRUE)
  point <- predict(
    fit,
    type = "msm_point",
    at = c(psi1 = 1, psi2 = 2),
    interval = TRUE
  )
  expect_identical(surface$prediction_type, "msm_surface")
  expect_identical(surface$estimand_scale, "response")
  expect_identical(surface$estimate_scale, "response")
  expect_identical(surface$fit_estimand_scale, "odds_ratio")
  expect_true(all(surface$estimates$estimate > 0))
  expect_true(all(surface$estimates$estimate < 1))
  expect_true(all(surface$intervals$lower > 0))
  expect_true(all(surface$intervals$upper < 1))
  expect_identical(point$prediction_type, "msm_point")
  expect_identical(point$estimate_scale, "response")
  expect_equal(nrow(point$estimates), 1)
  expect_true(point$estimates$estimate > 0)
  expect_true(point$estimates$estimate < 1)
})
test_that("predict MSM contrast defaults to response scale and supports estimand scale", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    family = poisson(link = "log")
  )
  from <- c(psi1 = 0, psi2 = 0)
  to <- c(psi1 = 2, psi2 = 1)
  from_grid <- qgcompmulti_build_msm_grid(fit, at = from)$grid
  to_grid <- qgcompmulti_build_msm_grid(fit, at = to)$grid
  from_target <- qgcompmulti_predict_msm_target(fit, from_grid)
  to_target <- qgcompmulti_predict_msm_target(fit, to_grid)
  default_contrast <- predict(
    fit,
    type = "msm_contrast",
    from = from,
    to = to
  )
  ratio_contrast <- predict(
    fit,
    type = "msm_contrast",
    from = from,
    to = to,
    contrast_scale = "estimand",
    interval = TRUE
  )
  expect_identical(default_contrast$estimand_scale, "response")
  expect_identical(default_contrast$estimate_scale, "response")
  expect_identical(default_contrast$contrast_scale, "response")
  expect_equal(default_contrast$estimates$estimate, exp(to_target) - exp(from_target))
  expect_identical(ratio_contrast$estimand_scale, "rate_ratio")
  expect_identical(ratio_contrast$estimate_scale, "estimand")
  expect_identical(ratio_contrast$contrast_scale, "estimand")
  expect_equal(ratio_contrast$estimates$estimate, exp(to_target - from_target))
  expect_true(all(ratio_contrast$intervals$lower > 0))
  expect_true(all(ratio_contrast$intervals$upper > 0))
})
test_that("predict validates contrast scale usage", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  dat <- make_test_data()
  expect_error(
    predict(fit, type = "msm", contrast_scale = "estimand"),
    "only supported"
  )
  expect_error(
    predict(
      fit,
      type = "exact_contrast",
      data = dat,
      from = c(psi1 = 0, psi2 = 0),
      to = c(psi1 = 1, psi2 = 2),
      contrast_scale = "estimand"
    ),
    "only supported"
  )
  expect_error(
    predict(
      fit,
      type = "msm_contrast",
      from = c(psi1 = 0, psi2 = 0),
      to = c(psi1 = 1, psi2 = 2),
      contrast_scale = "not_a_scale"
    ),
    "should be one of"
  )
})
