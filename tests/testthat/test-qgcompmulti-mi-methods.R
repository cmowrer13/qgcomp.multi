test_that("summary and print methods expose the pooled MI public interface", {
  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = make_completed_data_list(m = 3, seed = 321),
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 5,
    seed = 812,
    keep_fits = FALSE
  )

  s <- summary(fit)
  expect_s3_class(s, "summary.qgcompmulti_mi")
  expect_identical(
    names(s),
    c(
      "call",
      "formula",
      "fit_overview",
      "mi_overview",
      "mixtures",
      "msm_table",
      "pooling_table",
      "labels"
    )
  )

  printed_fit <- withVisible(print(fit))
  expect_false(printed_fit$visible)
  expect_identical(printed_fit$value, fit)
  expect_output(print(fit), "qgcompmulti multiple-imputation fit")
  expect_output(print(fit), "Multiple imputation")
  expect_output(print(fit), "Imputed datasets")
  expect_output(print(fit), "MSM coefficients")
  expect_output(print(fit), "Estimand")
  expect_output(print(fit), "Interval method: wald")

  printed_summary <- withVisible(print(s))
  expect_false(printed_summary$visible)
  expect_identical(printed_summary$value, s)
  expect_output(print(s), "Summary of qgcompmulti multiple-imputation fit")
  expect_output(print(s), "Multiple imputation overview")
  expect_output(print(s), "Estimand")
  expect_output(print(s), "Interval method: wald")
  expect_output(print(s), "Pooling diagnostics")
})

test_that("pooled MI extractor methods agree with stored results", {
  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = make_completed_data_list(m = 3, seed = 654),
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = FALSE,
    q = NULL,
    centering = "median",
    B = 5,
    keep_fits = TRUE,
    seed = 415
  )

  expect_equal(coef(fit), fit$results$coefficients)
  expect_equal(vcov(fit), fit$results$vcov)

  ci_all <- confint(fit)
  ci_name <- confint(fit, parm = "psi1")
  ci_pos <- confint(fit, parm = 2)

  expect_true(is.matrix(ci_all))
  expect_identical(rownames(ci_all), names(coef(fit)))
  expect_identical(rownames(ci_name), "psi1")
  expect_identical(rownames(ci_pos), "psi1")
  expect_equal(ci_name, ci_pos)
  expect_error(confint(fit, parm = "not_a_parameter"))
  expect_error(confint(fit, level = 1))
  expect_error(confint(fit, method = "percentile"), "supports only")
})

test_that("pooled MI tidy and glance methods return the expected metadata", {
  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = make_completed_data_list(m = 3, seed = 321, clustered = TRUE),
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 5,
    id = "cluster_id",
    seed = 511
  )

  tidy_method <- getFromNamespace("tidy.qgcompmulti_mi", "qgcomp.multi")
  glance_method <- getFromNamespace("glance.qgcompmulti_mi", "qgcomp.multi")

  td <- tidy_method(fit)
  expect_true(is.data.frame(td))
  expect_identical(
    names(td),
    c(
      "term", "estimate", "display.estimate", "std.error", "statistic",
      "df", "p.value", "estimand_scale", "msm_fitting_scale",
      "estimate_scale", "display_scale"
    )
  )
  expect_identical(td$term, names(coef(fit)))
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(td$display.estimate, unname(coef(fit)))
  expect_equal(td$std.error, unname(fit$results$std_error))
  expect_equal(td$df, unname(fit$results$df))
  expect_true(all(td$estimate_scale == "fitting"))

  td_ci <- tidy_method(fit, conf.int = TRUE, conf.level = 0.95)
  expect_identical(
    names(td_ci),
    c(
      "term", "estimate", "display.estimate", "std.error", "statistic",
      "df", "p.value", "estimand_scale", "msm_fitting_scale",
      "estimate_scale", "display_scale", "conf.low", "conf.high",
      "display.conf.low", "display.conf.high"
    )
  )

  gl <- glance_method(fit)
  expect_true(is.data.frame(gl))
  expect_equal(nrow(gl), 1L)
  expect_identical(
    names(gl),
    c(
      "n_input", "n_used", "m", "input_type", "family", "link",
      "estimand_scale", "msm_fitting_scale", "default_interval_method",
      "quantized", "q", "centering", "MCsize", "interaction",
      "keep_fits", "has_clusters", "cluster_var", "n_clusters"
    )
  )
  expect_equal(gl$m, fit$mi_info$m)
  expect_identical(gl$input_type, fit$mi_info$input_type)
  expect_true(gl$has_clusters)
  expect_identical(gl$cluster_var, "cluster_id")
  expect_equal(gl$n_clusters, fit$data_info$n_clusters)
})

test_that("broom generic dispatch works for pooled MI objects when broom is installed", {
  skip_if_not_installed("broom")

  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = make_completed_data_list(m = 3, seed = 321),
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 5,
    seed = 812
  )

  tidy_method <- getFromNamespace("tidy.qgcompmulti_mi", "qgcomp.multi")
  glance_method <- getFromNamespace("glance.qgcompmulti_mi", "qgcomp.multi")

  expect_equal(broom::tidy(fit), tidy_method(fit))
  expect_equal(broom::glance(fit), glance_method(fit))
})

test_that("df.residual supports scientifically reasonable mice pooling workflows", {
  skip_if_not_installed("mice")

  fit <- fit_test_model(interaction = TRUE, q = 4, B = 5, seed = 111)
  expect_equal(df.residual(fit), nobs(fit) - length(coef(fit)))
  expect_equal(residuals(fit), stats::residuals(fit$fits$outcome_fit))

  fits <- fit_test_model_list(B = 5)
  pooled <- mice::pool(fits)
  pooled_summary <- summary(pooled)

  expect_s3_class(pooled, "mipo")
  expect_identical(as.character(pooled_summary$term), names(coef(fits[[1]])))
  expect_equal(
    as.numeric(unique(pooled$pooled$dfcom)),
    as.numeric(df.residual(fits[[1]]))
  )
  expect_true(all(is.finite(pooled_summary$df)))
  expect_true(all(pooled_summary$df > 0))
})

test_that("with.mids and pool can produce a simple pooled coefficient table", {
  skip_if_not_installed("mice")

  mids <- make_test_mids(m = 3, seed = 888)
  mira_fit <- with(
    mids,
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = data.frame(Y, X1, X2, X3, W1, W2, W3, C),
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = TRUE,
      q = 4,
      B = 5
    )
  )
  pooled <- mice::pool(mira_fit)
  pooled_summary <- summary(pooled)

  expect_s3_class(mira_fit, "mira")
  expect_s3_class(pooled, "mipo")
  expect_identical(
    as.character(pooled_summary$term),
    c("(Intercept)", "psi1", "psi2", "psi1:psi2")
  )
  expect_true(all(c("estimate", "std.error", "statistic", "df", "p.value") %in% names(pooled_summary)))
})

test_that("pooled MI ratio reporting transforms final pooled estimates only", {
  completed <- lapply(c(432, 433, 434), make_binomial_test_data)
  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    family = binomial(link = "logit"),
    B = 5,
    seed = 912
  )
  tidy_method <- getFromNamespace("tidy.qgcompmulti_mi", "qgcomp.multi")
  td <- tidy_method(fit, conf.int = TRUE)
  fitting_ci <- build_qgcompmulti_mi_confint(
    coefficients = coef(fit),
    std_error = fit$results$std_error[names(coef(fit))],
    df = fit$results$df[names(coef(fit))],
    level = 0.95
  )
  expect_identical(fit$analysis$estimand_scale, "odds_ratio")
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(td$display.estimate, unname(exp(coef(fit))))
  expect_equal(td$conf.low, unname(fitting_ci[, 1]))
  expect_equal(td$conf.high, unname(fitting_ci[, 2]))
  expect_equal(td$display.conf.low, unname(exp(fitting_ci[, 1])))
  expect_equal(td$display.conf.high, unname(exp(fitting_ci[, 2])))
  expect_equal(confint(fit), exp(fitting_ci))
})
