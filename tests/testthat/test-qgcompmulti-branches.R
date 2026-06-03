test_that("interaction and quantized branch produces four MSM coefficients", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_s3_class(fit, "qgcompmulti")
  expect_identical(names(coef(fit)), EXPECTED_COEF_NAMES_WITH_INTERACTION)
  expect_identical(dim(vcov(fit)), c(4L, 4L))
  expect_s3_class(summary(fit), "summary.qgcompmulti")
})

test_that("no-interaction and quantized branch produces three MSM coefficients", {
  fit <- fit_test_model(interaction = FALSE, q = 4)
  expect_s3_class(fit, "qgcompmulti")
  expect_identical(names(coef(fit)), EXPECTED_COEF_NAMES_NO_INTERACTION)
  expect_identical(dim(vcov(fit)), c(3L, 3L))
  expect_s3_class(summary(fit), "summary.qgcompmulti")
})

test_that("interaction and original-scale branch works with q = NULL", {
  fit <- fit_test_model(interaction = TRUE, q = NULL, centering = "median")
  expect_s3_class(fit, "qgcompmulti")
  expect_null(fit$mixtures$q)
  expect_identical(names(coef(fit)), EXPECTED_COEF_NAMES_WITH_INTERACTION)
  expect_identical(dim(vcov(fit)), c(4L, 4L))
  expect_s3_class(summary(fit), "summary.qgcompmulti")
})

test_that("no-interaction and original-scale branch works with q = NULL", {
  fit <- fit_test_model(interaction = FALSE, q = NULL, centering = "median")
  expect_s3_class(fit, "qgcompmulti")
  expect_null(fit$mixtures$q)
  expect_identical(names(coef(fit)), EXPECTED_COEF_NAMES_NO_INTERACTION)
  expect_identical(dim(vcov(fit)), c(3L, 3L))
  expect_s3_class(summary(fit), "summary.qgcompmulti")
})
