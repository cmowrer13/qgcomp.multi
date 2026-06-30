test_that("support diagnostics summarize stored intervention support", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  diag <- support(fit)
  expect_s3_class(diag, "qgcompmulti_support_diagnostic")
  expect_identical(diag$diagnostic_type, "support")
  expect_identical(diag$mode, "quantized")
  expect_equal(
    diag$support_summary$psi1_values,
    sort(unique(fit$prediction$intervention_grid$psi1))
  )
  expect_equal(
    diag$support_summary$psi2_values,
    sort(unique(fit$prediction$intervention_grid$psi2))
  )
})
test_that("original-scale support diagnostics flag pooled percentile support", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 10
  )
  diag <- diagnostics(fit, type = "support")
  expect_s3_class(diag, "qgcompmulti_support_diagnostic")
  expect_identical(diag$mode, "original_scale")
  expect_true(diag$flags$pooled_percentile_grid)
  expect_true(diag$flags$centered_msm_grid)
  expect_equal(
    diag$support_summary$intervention_psi1_range,
    range(fit$prediction$intervention_grid$psi1)
  )
})
test_that("bootstrap diagnostics summarize retained metadata", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  diag <- diagnostics(fit, type = "bootstrap")
  expect_s3_class(diag, "qgcompmulti_bootstrap_diagnostic")
  expect_identical(diag$B_requested, fit$bootstrap$B_requested)
  expect_identical(diag$B_success, fit$bootstrap$B_success)
  expect_identical(diag$B_failed, fit$bootstrap$B_failed)
  expect_equal(diag$success_rate, diag$B_success / diag$B_requested)
})
test_that("bootstrap diagnostics summarize failure messages when present", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  fit$bootstrap$B_requested <- fit$bootstrap$B_success + 1L
  fit$bootstrap$B_failed <- 1L
  fit$bootstrap$failure_log <- data.frame(
    replicate = 11L,
    message = "glm failed",
    row.names = NULL
  )
  diag <- diagnostics(fit, type = "bootstrap")
  expect_true(diag$flags$any_failures)
  expect_equal(diag$failure_summary$message, "glm failed")
  expect_equal(diag$failure_summary$n, 1L)
})
test_that("adequacy diagnostics summarize stored surface comparison", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  diag <- adequacy(fit)
  comparison <- fit$prediction$surface_comparison
  residual <- comparison$residual
  expect_s3_class(diag, "qgcompmulti_adequacy_diagnostic")
  expect_identical(diag$diagnostic_type, "adequacy")
  expect_equal(diag$comparison, comparison)
  expect_equal(diag$summary_metrics$mae, mean(abs(residual)))
  expect_equal(diag$summary_metrics$rmse, sqrt(mean(residual^2)))
  expect_equal(diag$summary_metrics$max_abs_error, max(abs(residual)))
  expect_equal(
    diag$summary_metrics$correlation,
    stats::cor(comparison$exact_mean, comparison$msm_mean)
  )
  expect_identical(diag$comparison_scale, "response")
  expect_identical(diag$msm_fitting_scale, "identity")
})
test_that("adequacy diagnostics use fitting-scale comparisons for odds-ratio fits", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    B = 10,
    family = binomial(link = "logit")
  )
  diag <- adequacy(fit)
  comparison <- fit$prediction$surface_comparison_target
  residual <- comparison$residual_target
  expect_s3_class(diag, "qgcompmulti_adequacy_diagnostic")
  expect_identical(diag$comparison_scale, "fitting")
  expect_identical(diag$msm_fitting_scale, "logit")
  expect_equal(diag$comparison, comparison)
  expect_equal(diag$response_comparison, fit$prediction$surface_comparison)
  expect_equal(diag$summary_metrics$mae, mean(abs(residual)))
  expect_equal(diag$summary_metrics$rmse, sqrt(mean(residual^2)))
  expect_equal(diag$summary_metrics$max_abs_error, max(abs(residual)))
  expect_equal(
    diag$summary_metrics$correlation,
    stats::cor(comparison$exact_target, comparison$msm_target)
  )
})
test_that("adequacy diagnostics use fitting-scale comparisons for rate-ratio fits", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    B = 10,
    family = poisson(link = "log")
  )
  diag <- adequacy(fit)
  comparison <- fit$prediction$surface_comparison_target
  residual <- comparison$residual_target
  expect_identical(diag$comparison_scale, "fitting")
  expect_identical(diag$msm_fitting_scale, "log")
  expect_equal(diag$comparison, comparison)
  expect_equal(diag$summary_metrics$mae, mean(abs(residual)))
  expect_equal(
    diag$summary_metrics$correlation,
    stats::cor(comparison$exact_target, comparison$msm_target)
  )
})
test_that("combined diagnostics return all three diagnostic families", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  diag <- diagnostics(fit, type = "all")
  expect_s3_class(diag, "qgcompmulti_diagnostics")
  expect_named(diag, c("support", "bootstrap", "adequacy"))
  expect_s3_class(diag$support, "qgcompmulti_support_diagnostic")
  expect_s3_class(diag$bootstrap, "qgcompmulti_bootstrap_diagnostic")
  expect_s3_class(diag$adequacy, "qgcompmulti_adequacy_diagnostic")
})
test_that("diagnostics also work for no-interaction fits", {
  fit <- fit_test_model(interaction = FALSE, q = 4, B = 10)
  diag <- diagnostics(fit, type = "all")
  expect_s3_class(diag$support, "qgcompmulti_support_diagnostic")
  expect_s3_class(diag$bootstrap, "qgcompmulti_bootstrap_diagnostic")
  expect_s3_class(diag$adequacy, "qgcompmulti_adequacy_diagnostic")
})
test_that("diagnostic objects print cleanly", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  expect_output(print(support(fit)), "intervention support diagnostic")
  expect_output(print(diagnostics(fit, type = "bootstrap")), "bootstrap diagnostic")
  expect_output(print(adequacy(fit)), "MSM adequacy diagnostic")
  expect_output(print(adequacy(fit)), "Comparison scale: response")
  expect_output(print(diagnostics(fit, type = "all")), "qgcompmulti diagnostics")
})
test_that("transformed-scale adequacy diagnostic prints the fitting scale", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    B = 10,
    family = binomial(link = "logit")
  )
  expect_output(print(adequacy(fit)), "Comparison scale: fitting")
  expect_output(print(adequacy(fit)), "MSM fitting scale: logit")
})

