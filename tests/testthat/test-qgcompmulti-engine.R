test_that("qgcompmulti_msm_fit returns the expected component structure", {
  engine <- fit_engine_result(interaction = TRUE, q = 4)
  expect_true(is.list(engine))
  expect_identical(
    names(engine),
    c("outcome_fit", "msm_fit", "coefficients", "n_used")
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
