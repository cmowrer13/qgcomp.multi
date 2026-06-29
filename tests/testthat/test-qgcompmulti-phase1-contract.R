test_that("estimand helpers distinguish planned defaults from current stored defaults", {
  planned_binomial <- qgcompmulti_resolve_estimand_spec(
    family = binomial(link = "logit"),
    mode = "planned"
  )
  current_binomial <- qgcompmulti_resolve_estimand_spec(
    family = binomial(link = "logit"),
    mode = "current"
  )
  planned_poisson <- qgcompmulti_resolve_estimand_spec(
    family = poisson(link = "log"),
    mode = "planned"
  )

  expect_identical(planned_binomial$estimand_scale, "odds_ratio")
  expect_identical(planned_binomial$msm_fitting_scale, "logit")
  expect_true(planned_binomial$estimand_scale_defaulted)

  expect_identical(current_binomial$estimand_scale, "risk_difference")
  expect_identical(current_binomial$msm_fitting_scale, "identity")

  expect_identical(planned_poisson$estimand_scale, "rate_ratio")
  expect_identical(planned_poisson$msm_fitting_scale, "log")
})

test_that("estimand helper validation enforces family and link compatibility", {
  expect_error(
    qgcompmulti_validate_estimand_scale(
      family_name = "gaussian",
      link = "identity",
      estimand_scale = "odds_ratio"
    ),
    "not supported"
  )

  expect_error(
    qgcompmulti_validate_estimand_scale(
      family_name = "binomial",
      link = "probit",
      estimand_scale = "odds_ratio"
    ),
    "not supported"
  )

  expect_error(
    qgcompmulti_validate_estimand_scale(
      family_name = "poisson",
      link = "sqrt",
      estimand_scale = "rate_ratio"
    ),
    "not supported"
  )
})

test_that("interval helpers expose the planned method contract", {
  expect_identical(
    qgcompmulti_interval_methods("single_fit"),
    c("wald", "percentile", "percentile_type2")
  )
  expect_identical(qgcompmulti_interval_methods("mi"), "wald")
  expect_identical(
    qgcompmulti_prediction_interval_type_for_method("percentile"),
    "bootstrap_percentile"
  )
  expect_identical(
    qgcompmulti_prediction_interval_type_for_method("percentile_type2"),
    "bootstrap_percentile_type2"
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

test_that("region helpers require at least two named coefficients", {
  fit <- fit_test_model(interaction = TRUE, q = 4)

  expect_identical(
    qgcompmulti_resolve_region_parm(c("psi1", "psi2"), names(coef(fit))),
    c("psi1", "psi2")
  )

  expect_error(
    qgcompmulti_resolve_region_parm("psi1", names(coef(fit))),
    "at least two"
  )

  covariance <- diag(c(0.04, 0.03), nrow = 2L)
  rownames(covariance) <- c("psi1", "psi2")
  colnames(covariance) <- c("psi1", "psi2")

  region <- new_qgcompmulti_region(
    parm = c("psi1", "psi2"),
    center = c(psi1 = 0.2, psi2 = 0.1),
    covariance = covariance
  )

  expect_s3_class(region, "qgcompmulti_region")
})

test_that("transformed MSM surface helpers reject boundary values for ratio scales", {
  expect_error(
    qgcompmulti_transform_msm_surface(
      values = c(0, 0.25, 0.75),
      msm_fitting_scale = "logit",
      direction = "to_fitting"
    ),
    "strictly between 0 and 1"
  )
  expect_error(
    qgcompmulti_transform_msm_surface(
      values = c(0, 0.5, 1.25),
      msm_fitting_scale = "log",
      direction = "to_fitting"
    ),
    "strictly positive"
  )
})
