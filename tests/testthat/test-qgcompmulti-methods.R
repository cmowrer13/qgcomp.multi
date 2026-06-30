test_that("core extractor methods return the expected types", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_type(coef(fit), "double")
  expect_true(is.matrix(vcov(fit)))
  expect_true(is.matrix(confint(fit)))
  expect_type(residuals(fit), "double")
  expect_s3_class(formula(fit), "formula")
  expect_equal(nobs(fit), fit$data_info$n_used)
  expect_type(df.residual(fit), "double")
})

test_that("extractor methods agree with stored results", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_equal(coef(fit), fit$results$coefficients)
  expect_equal(vcov(fit), fit$results$vcov)
  expect_equal(residuals(fit), stats::residuals(fit$fits$outcome_fit))
  expect_equal(nobs(fit), fit$data_info$n_used)
  expect_equal(df.residual(fit), nobs(fit) - length(coef(fit)))
})

test_that("coefficient, covariance, and confidence interval outputs are aligned", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  coef_names <- names(coef(fit))
  vc <- vcov(fit)
  ci <- confint(fit)
  expect_identical(coef_names, rownames(vc))
  expect_identical(coef_names, colnames(vc))
  expect_identical(coef_names, rownames(ci))
  expect_identical(colnames(ci), c("2.5 %", "97.5 %"))
})

test_that("confidence interval subsetting works by name and position", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  ci_name <- confint(fit, parm = "psi1")
  ci_pos <- confint(fit, parm = 2)
  expect_true(is.matrix(ci_name))
  expect_true(is.matrix(ci_pos))
  expect_identical(rownames(ci_name), "psi1")
  expect_identical(rownames(ci_pos), "psi1")
  expect_equal(ci_name, ci_pos)
})

test_that("confidence interval validation catches bad inputs", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_error(confint(fit, parm = "not_a_parameter"))
  expect_error(confint(fit, parm = 99))
  expect_error(confint(fit, level = 1))
  expect_error(confint(fit, level = 0))
  expect_error(confint(fit, level = NA_real_))
})

test_that("confidence intervals display ratio estimands on the estimand scale", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    family = binomial(link = "logit")
  )
  coef_names <- names(coef(fit))
  fitting_ci <- build_qgcompmulti_confint(
    coefficients = coef(fit),
    std_error = setNames(sqrt(diag(vcov(fit))), coef_names),
    level = 0.95
  )
  display_ci <- confint(fit, method = "wald")
  expect_identical(fit$analysis$estimand_scale, "odds_ratio")
  expect_equal(display_ci, exp(fitting_ci))
  expect_true(all(display_ci > 0))
})
test_that("single-fit confidence interval method override is explicit while non-Wald reporting is deferred", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    default_interval_method = "percentile_type2"
  )
  expect_error(confint(fit), "currently supports only")
  expect_true(is.matrix(confint(fit, method = "wald")))
})

test_that("ratio-scale fits keep coefficients on the MSM fitting scale", {
  fit_binomial <- fit_test_model(
    interaction = TRUE,
    q = 4,
    family = binomial(link = "logit")
  )
  fit_poisson <- fit_test_model(
    interaction = TRUE,
    q = 4,
    family = poisson(link = "log")
  )
  expect_identical(fit_binomial$analysis$estimand_scale, "odds_ratio")
  expect_identical(fit_binomial$analysis$msm_fitting_scale, "logit")
  expect_equal(
    unname(coef(fit_binomial)),
    unname(stats::coef(fit_binomial$fits$msm_fit)[names(coef(fit_binomial))])
  )
  expect_identical(fit_poisson$analysis$estimand_scale, "rate_ratio")
  expect_identical(fit_poisson$analysis$msm_fitting_scale, "log")
  expect_equal(
    unname(coef(fit_poisson)),
    unname(stats::coef(fit_poisson$fits$msm_fit)[names(coef(fit_poisson))])
  )
})

test_that("explicit additive estimands remain available for binomial and Poisson fits", {
  fit_binomial <- fit_test_model(
    interaction = TRUE,
    q = 4,
    family = binomial(link = "logit"),
    estimand_scale = "risk_difference"
  )
  fit_poisson <- fit_test_model(
    interaction = TRUE,
    q = 4,
    family = poisson(link = "log"),
    estimand_scale = "mean_difference"
  )
  expect_identical(fit_binomial$analysis$estimand_scale, "risk_difference")
  expect_identical(fit_binomial$analysis$msm_fitting_scale, "identity")
  expect_equal(
    fit_binomial$prediction$counterfactual_surface$exact_mean,
    fit_binomial$prediction$counterfactual_surface_target$exact_target
  )
  expect_identical(fit_poisson$analysis$estimand_scale, "mean_difference")
  expect_identical(fit_poisson$analysis$msm_fitting_scale, "identity")
  expect_equal(
    fit_poisson$prediction$counterfactual_surface$exact_mean,
    fit_poisson$prediction$counterfactual_surface_target$exact_target
  )
})
test_that("supported default interval methods are accepted and stored", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    default_interval_method = "percentile_type2"
  )
  expect_identical(fit$analysis$default_interval_method, "percentile_type2")
})
