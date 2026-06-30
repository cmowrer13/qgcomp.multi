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

test_that("interval helpers expose the final method contract", {
  expect_identical(
    qgcompmulti_interval_methods("single_fit"),
    c("wald", "percentile", "basic")
  )
  expect_identical(qgcompmulti_interval_methods("mi"), "wald")
  expect_identical(
    qgcompmulti_prediction_interval_type_for_method("percentile"),
    "bootstrap_percentile"
  )
  expect_identical(
    qgcompmulti_prediction_interval_type_for_method("basic"),
    "bootstrap_basic"
  )
  expect_error(
    qgcompmulti_validate_interval_method("percentile", context = "mi"),
    "must be one of"
  )
  expect_error(
    qgcompmulti_validate_prediction_interval_type("not_a_type"),
    "must be `NULL` or one of"
  )
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
test_that("basic bootstrap intervals display ratio estimands on the estimand scale", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    B = 10,
    family = binomial(link = "logit")
  )
  fitting_ci <- qgcompmulti_build_single_fit_confint(
    coefficients = coef(fit),
    std_error = setNames(sqrt(diag(vcov(fit))), names(coef(fit))),
    coef_draws = fit$bootstrap$coef_draws,
    level = 0.90,
    method = "basic"
  )
  display_ci <- confint(fit, level = 0.90, method = "basic")
  expect_equal(display_ci, exp(fitting_ci))
  expect_true(all(display_ci > 0))
})
test_that("single-fit confidence intervals support Wald, percentile, and basic methods", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  coef_names <- names(coef(fit))
  draws <- as.matrix(fit$bootstrap$coef_draws)[, coef_names, drop = FALSE]
  level <- 0.80
  alpha <- (1 - level) / 2
  quantiles <- t(apply(
    draws,
    2,
    stats::quantile,
    probs = c(alpha, 1 - alpha),
    na.rm = TRUE,
    names = FALSE
  ))
  percentile_ci <- confint(fit, level = level, method = "percentile")
  basic_ci <- confint(fit, level = level, method = "basic")
  expected_basic <- cbind(
    2 * coef(fit) - quantiles[, 2],
    2 * coef(fit) - quantiles[, 1]
  )
  rownames(expected_basic) <- coef_names
  colnames(expected_basic) <- qgcompmulti_confint_colnames(level)
  expect_true(is.matrix(confint(fit, method = "wald")))
  expect_equal(unname(percentile_ci), unname(quantiles))
  expect_equal(basic_ci, expected_basic)
})
test_that("single-fit confidence interval default method can use basic bootstrap", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    B = 10,
    default_interval_method = "basic"
  )
  expect_identical(fit$analysis$default_interval_method, "basic")
  expect_equal(confint(fit), confint(fit, method = "basic"))
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
    default_interval_method = "basic"
  )
  expect_identical(fit$analysis$default_interval_method, "basic")
})
