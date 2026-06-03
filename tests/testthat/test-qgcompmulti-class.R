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
  expect_true(all(vapply(fit$prediction, is.null, logical(1))))
  expect_true(is.matrix(fit$bootstrap$coef_draws))
  expect_equal(fit$bootstrap$B_requested, 10)
  expect_equal(fit$bootstrap$B_success, 10)
  expect_equal(fit$bootstrap$B_failed, 0)
  expect_null(fit$bootstrap$failure_log)
  expect_true(is.character(fit$labels$coefficient_names))
  expect_true(is.character(fit$labels$coefficient_labels))
  expect_true(is.character(fit$labels$mixture_labels))
  expect_identical(names(fit$labels$mixture_labels), c("mix1", "mix2"))
})
