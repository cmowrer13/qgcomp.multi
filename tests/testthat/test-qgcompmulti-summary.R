test_that("summary returns the expected structured summary object", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  s <- summary(fit)
  expect_s3_class(s, "summary.qgcompmulti")
  expect_identical(names(s), EXPECTED_SUMMARY_COMPONENTS)
  expect_true(is.list(s$fit_overview))
  expect_true(is.list(s$mixtures))
  expect_true(is.data.frame(s$msm_table))
  expect_true(is.list(s$outcome_model_info))
})

test_that("print.qgcompmulti returns invisibly and shows key sections", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  printed <- withVisible(print(fit))
  expect_false(printed$visible)
  expect_identical(printed$value, fit)
  expect_output(print(fit), "qgcompmulti fit")
  expect_output(print(fit), "MSM coefficients")
  expect_output(print(fit), "Mixtures")
  expect_output(print(fit), "Estimand")
  expect_output(print(fit), "MSM fitting scale")
  expect_output(print(fit), "Default interval method")
})

test_that("print.summary.qgcompmulti returns invisibly and shows key sections", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  s <- summary(fit)
  printed <- withVisible(print(s))
  expect_false(printed$visible)
  expect_identical(printed$value, s)
  expect_output(print(s), "Summary of qgcompmulti fit")
  expect_output(print(s), "Model overview")
  expect_output(print(s), "Estimand")
  expect_output(print(s), "MSM fitting scale")
  expect_output(print(s), "Outcome model context")
})

test_that("ratio-scale print and summary tables include both display and fitting-scale estimates", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    family = binomial(link = "logit")
  )
  s <- summary(fit)
  expect_true("Odds ratio" %in% names(s$msm_table))
  expect_true("Coefficient" %in% names(s$msm_table))
  expect_true("Std. Error" %in% names(s$msm_table))
  expect_equal(
    unname(s$msm_table[["Odds ratio"]]),
    unname(exp(fit$results$coef_table$Estimate))
  )
  expect_equal(
    unname(s$msm_table[["Coefficient"]]),
    unname(fit$results$coef_table$Estimate)
  )
  expect_output(print(fit), "Odds ratio")
  expect_output(print(fit), "Coefficient")
})
