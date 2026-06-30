test_that("tidy method returns the expected coefficient-centric columns", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  tidy_method <- getFromNamespace("tidy.qgcompmulti", "qgcomp.multi")
  td <- tidy_method(fit)
  expect_true(is.data.frame(td))
  expect_identical(
    names(td),
    c(
      "term", "estimate", "display.estimate", "std.error", "statistic",
      "p.value", "estimand_scale", "msm_fitting_scale", "estimate_scale",
      "display_scale"
    )
  )
  expect_identical(td$term, names(coef(fit)))
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(td$display.estimate, unname(coef(fit)))
  expect_equal(td$std.error, unname(sqrt(diag(vcov(fit)))))
  expect_equal(td$statistic, td$estimate / td$std.error)
  expect_equal(
    td$p.value,
    unname(2 * stats::pnorm(abs(td$statistic), lower.tail = FALSE))
  )
  expect_true(all(td$estimate_scale == "fitting"))
})
test_that("tidy method adds confidence interval columns when requested", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  tidy_method <- getFromNamespace("tidy.qgcompmulti", "qgcomp.multi")
  td <- tidy_method(fit, conf.int = TRUE, conf.level = 0.95)
  ci <- confint(fit)
  expect_identical(
    names(td),
    c(
      "term", "estimate", "display.estimate", "std.error", "statistic",
      "p.value", "estimand_scale", "msm_fitting_scale", "estimate_scale",
      "display_scale", "conf.low", "conf.high", "display.conf.low",
      "display.conf.high"
    )
  )
  expect_equal(td$conf.low, unname(ci[, 1]))
  expect_equal(td$conf.high, unname(ci[, 2]))
  expect_equal(td$display.conf.low, unname(ci[, 1]))
  expect_equal(td$display.conf.high, unname(ci[, 2]))
  expect_error(tidy_method(fit, conf.int = NA))
  expect_error(tidy_method(fit, conf.int = TRUE, conf.level = 1))
})
test_that("tidy method handles no-interaction and q = NULL fits", {
  tidy_method <- getFromNamespace("tidy.qgcompmulti", "qgcomp.multi")
  fit_no_int <- fit_test_model(interaction = FALSE, q = 4)
  fit_cont <- fit_test_model(interaction = TRUE, q = NULL, centering = "median")
  expect_identical(tidy_method(fit_no_int)$term, EXPECTED_COEF_NAMES_NO_INTERACTION)
  expect_identical(tidy_method(fit_cont)$term, EXPECTED_COEF_NAMES_WITH_INTERACTION)
})
test_that("glance method returns compact fit metadata for q = NULL fits", {
  fit <- fit_test_model(interaction = TRUE, q = NULL, centering = "median")
  glance_method <- getFromNamespace("glance.qgcompmulti", "qgcomp.multi")
  gl <- glance_method(fit)
  expect_true(is.data.frame(gl))
  expect_equal(nrow(gl), 1L)
  expect_identical(
    names(gl),
    c(
      "n_input", "n_used", "family", "link", "estimand_scale",
      "msm_fitting_scale", "default_interval_method", "quantized", "q",
      "centering", "B_requested", "B_success", "B_failed", "MCsize",
      "interaction", "has_clusters", "cluster_var", "n_clusters"
    )
  )
  expect_equal(gl$n_input, fit$data_info$n_input)
  expect_equal(gl$n_used, fit$data_info$n_used)
  expect_false(gl$quantized)
  expect_true(is.na(gl$q))
  expect_identical(gl$centering, "median")
  expect_equal(gl$B_requested, fit$bootstrap$B_requested)
  expect_equal(gl$B_success, fit$bootstrap$B_success)
  expect_equal(gl$B_failed, fit$bootstrap$B_failed)
  expect_false(gl$has_clusters)
  expect_true(is.na(gl$cluster_var))
  expect_true(is.na(gl$n_clusters))
})
test_that("glance method returns clustering metadata when clustering is used", {
  fit <- fit_test_model(interaction = TRUE, q = 4, clustered = TRUE)
  glance_method <- getFromNamespace("glance.qgcompmulti", "qgcomp.multi")
  gl <- glance_method(fit)
  expect_true(gl$has_clusters)
  expect_identical(gl$cluster_var, "cluster_id")
  expect_equal(gl$n_clusters, fit$data_info$n_clusters)
})
test_that("broom generic dispatch works when broom is installed", {
  skip_if_not_installed("broom")
  fit <- fit_test_model(interaction = TRUE, q = 4)
  tidy_method <- getFromNamespace("tidy.qgcompmulti", "qgcomp.multi")
  glance_method <- getFromNamespace("glance.qgcompmulti", "qgcomp.multi")
  expect_equal(broom::tidy(fit), tidy_method(fit))
  expect_equal(broom::glance(fit), glance_method(fit))
})

test_that("tidy exposes display-scale ratio columns while preserving fitting-scale pooling columns", {
  tidy_method <- getFromNamespace("tidy.qgcompmulti", "qgcomp.multi")
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    family = poisson(link = "log")
  )
  td <- tidy_method(fit, conf.int = TRUE, conf.level = 0.95)
  fitting_ci <- build_qgcompmulti_confint(
    coefficients = coef(fit),
    std_error = setNames(sqrt(diag(vcov(fit))), names(coef(fit))),
    level = 0.95
  )
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(td$display.estimate, unname(exp(coef(fit))))
  expect_equal(td$conf.low, unname(fitting_ci[, 1]))
  expect_equal(td$conf.high, unname(fitting_ci[, 2]))
  expect_equal(td$display.conf.low, unname(exp(fitting_ci[, 1])))
  expect_equal(td$display.conf.high, unname(exp(fitting_ci[, 2])))
  expect_true(all(td$estimand_scale == "rate_ratio"))
  expect_true(all(td$msm_fitting_scale == "log"))
})
test_that("tidy supports basic bootstrap interval columns for single fits", {
  tidy_method <- getFromNamespace("tidy.qgcompmulti", "qgcomp.multi")
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  td <- tidy_method(fit, conf.int = TRUE, conf.level = 0.90, method = "basic")
  fitting_ci <- qgcompmulti_build_single_fit_confint(
    coefficients = coef(fit),
    std_error = setNames(sqrt(diag(vcov(fit))), names(coef(fit))),
    coef_draws = fit$bootstrap$coef_draws,
    level = 0.90,
    method = "basic"
  )
  expect_equal(td$conf.low, unname(fitting_ci[, 1]))
  expect_equal(td$conf.high, unname(fitting_ci[, 2]))
  expect_equal(td$display.conf.low, unname(fitting_ci[, 1]))
  expect_equal(td$display.conf.high, unname(fitting_ci[, 2]))
})
