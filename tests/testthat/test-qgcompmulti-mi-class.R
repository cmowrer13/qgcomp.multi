test_that("qgcompmulti_pool_mi_fits returns a qgcompmulti_mi object with expected structure", {
  fits <- fit_test_model_list()

  pooled <- qgcompmulti_pool_mi_fits(
    fits = fits,
    input_type = "completed_list",
    keep_fits = FALSE,
    seed = 9001,
    fit_seeds = c(1101, 1202, 1303)
  )

  expect_s3_class(pooled, "qgcompmulti_mi")
  expect_identical(names(pooled), EXPECTED_MI_TOP_COMPONENTS)

  for (component_name in names(EXPECTED_MI_COMPONENT_FIELDS)) {
    expect_true(is.list(pooled[[component_name]]))
    expect_identical(
      names(pooled[[component_name]]),
      EXPECTED_MI_COMPONENT_FIELDS[[component_name]]
    )
  }

  expect_null(pooled$fits$imputation_fits)
  expect_identical(pooled$formula, fits[[1]]$formula)
  expect_equal(pooled$mi_info$m, 3)
  expect_identical(pooled$mi_info$input_type, "completed_list")
  expect_false(pooled$mi_info$keep_fits)
  expect_equal(pooled$mi_info$seed, 9001L)
  expect_equal(pooled$mi_info$fit_seeds, c(1101L, 1202L, 1303L))
  expect_identical(pooled$analysis$estimand_scale, "mean_difference")
  expect_true(pooled$analysis$estimand_scale_defaulted)
  expect_identical(pooled$analysis$msm_fitting_scale, "identity")
  expect_identical(pooled$analysis$default_interval_method, "wald")
})

test_that("Rubin pooling helper returns expected pooled variance components", {
  coefficients <- list(
    c("(Intercept)" = 1.00, psi1 = 0.50),
    c("(Intercept)" = 1.20, psi1 = 0.70),
    c("(Intercept)" = 0.80, psi1 = 0.60)
  )
  diag_matrix <- function(values) {
    out <- diag(values, nrow = length(values))
    dimnames(out) <- list(names(coefficients[[1]]), names(coefficients[[1]]))
    out
  }

  vcovs <- list(
    diag_matrix(c(0.04, 0.09)),
    diag_matrix(c(0.05, 0.10)),
    diag_matrix(c(0.03, 0.08))
  )

  pooled <- qgcompmulti_pool_rubin(coefficients = coefficients, vcovs = vcovs)

  coef_mat <- do.call(rbind, lapply(coefficients, unname))
  expected_coef <- setNames(colMeans(coef_mat), names(coefficients[[1]]))
  expected_within <- Reduce(`+`, vcovs) / 3
  expected_between <- stats::cov(coef_mat)
  dimnames(expected_between) <- list(names(expected_coef), names(expected_coef))
  expected_total <- expected_within + (1 + 1 / 3) * expected_between
  expected_riv <- ((1 + 1 / 3) * diag(expected_between)) / diag(expected_within)
  expected_df <- (3 - 1) * (1 + 1 / expected_riv)^2
  expected_fmi <- ((1 + 1 / 3) * diag(expected_between)) / diag(expected_total)

  expect_equal(pooled$coefficients, expected_coef)
  expect_equal(pooled$within_var, expected_within)
  expect_equal(pooled$between_var, expected_between)
  expect_equal(pooled$total_var, expected_total)
  expect_equal(pooled$vcov, expected_total)
  expect_equal(pooled$std_error, sqrt(diag(expected_total)), tolerance = 1e-10)
  expect_equal(pooled$riv, setNames(as.numeric(expected_riv), names(expected_coef)))
  expect_equal(pooled$df, setNames(as.numeric(expected_df), names(expected_coef)))
  expect_equal(pooled$fmi, setNames(as.numeric(expected_fmi), names(expected_coef)))
})

test_that("pooled MI coefficient table uses finite-df t inference and infinite-df fallback", {
  diag_matrix <- function(values) {
    out <- diag(values, nrow = length(values))
    dimnames(out) <- list(c("(Intercept)", "psi1"), c("(Intercept)", "psi1"))
    out
  }

  pooled_manual <- build_qgcompmulti_mi_results(
    list(
      coefficients = c("(Intercept)" = 1.0, psi1 = 0.5),
      std_error = c("(Intercept)" = 0.2, psi1 = 0.1),
      vcov = diag_matrix(c(0.04, 0.01)),
      within_var = diag_matrix(c(0.04, 0.01)),
      between_var = diag_matrix(c(0.01, 0)),
      total_var = diag_matrix(c(0.04, 0.01)),
      df = c("(Intercept)" = 12, psi1 = Inf),
      riv = c("(Intercept)" = 0.25, psi1 = 0),
      fmi = c("(Intercept)" = 0.2, psi1 = 0)
    )
  )

  stat <- c("(Intercept)" = 5, psi1 = 5)
  expected_p_intercept <- 2 * stats::pt(5, df = 12, lower.tail = FALSE)
  expected_p_psi1 <- 2 * stats::pnorm(5, lower.tail = FALSE)

  expect_equal(pooled_manual$coef_table$`t value`, unname(stat))
  expect_equal(pooled_manual$coef_table$df, c(12, Inf))
  expect_equal(pooled_manual$coef_table$`Pr(>|t|)`, c(expected_p_intercept, expected_p_psi1))
})

test_that("qgcompmulti_pool_mi_fits optionally retains per-imputation fits", {
  fits <- fit_test_model_list()

  pooled <- qgcompmulti_pool_mi_fits(
    fits = fits,
    keep_fits = TRUE,
    fit_seeds = c(1101, 1202, 1303)
  )

  expect_true(is.list(pooled$fits$imputation_fits))
  expect_length(pooled$fits$imputation_fits, 3)
  expect_true(all(vapply(pooled$fits$imputation_fits, inherits, logical(1), what = "qgcompmulti")))
  expect_true(is.finite(pooled$results$df[["(Intercept)"]]))
  expect_true(all(pooled$results$fmi >= 0))
})

test_that("qgcompmulti_pool_mi_fits handles no-interaction, q = NULL, and clustered fits", {
  fits_no_interaction <- fit_test_model_list(interaction = FALSE, fit_seeds = c(2101, 2202, 2303))
  pooled_no_interaction <- qgcompmulti_pool_mi_fits(
    fits = fits_no_interaction,
    fit_seeds = c(2101, 2202, 2303)
  )
  expect_identical(names(pooled_no_interaction$results$coefficients), EXPECTED_COEF_NAMES_NO_INTERACTION)

  fits_original_scale <- fit_test_model_list(
    interaction = TRUE,
    q = NULL,
    centering = "median",
    fit_seeds = c(3101, 3202, 3303)
  )
  pooled_original_scale <- qgcompmulti_pool_mi_fits(
    fits = fits_original_scale,
    fit_seeds = c(3101, 3202, 3303)
  )
  expect_null(pooled_original_scale$mixtures$q)
  expect_identical(pooled_original_scale$mixtures$centering, "median")
  expect_false(pooled_original_scale$data_info$quantized)

  fits_clustered <- fit_test_model_list(
    clustered = TRUE,
    fit_seeds = c(4101, 4202, 4303)
  )
  pooled_clustered <- qgcompmulti_pool_mi_fits(
    fits = fits_clustered,
    fit_seeds = c(4101, 4202, 4303)
  )
  expect_true(pooled_clustered$data_info$has_clusters)
  expect_identical(pooled_clustered$data_info$cluster_var, "cluster_id")
  expect_true(pooled_clustered$data_info$n_clusters >= 2L)
})

test_that("qgcompmulti_pool_mi_fits rejects incompatible per-imputation fits", {
  fit_a <- fit_test_model(interaction = TRUE, q = 4, seed = 111)
  fit_b <- fit_test_model(interaction = FALSE, q = 4, seed = 222)

  expect_error(
    qgcompmulti_pool_mi_fits(list(fit_a, fit_b)),
    "analysis\\$interaction"
  )
})

test_that("qgcompmulti_pool_mi_fits rejects non-Wald per-imputation interval defaults", {
  fits <- list(
    fit_test_model(B = 5, seed = 610, default_interval_method = "percentile"),
    fit_test_model(B = 5, seed = 620, default_interval_method = "percentile"),
    fit_test_model(B = 5, seed = 630, default_interval_method = "percentile")
  )

  expect_error(
    qgcompmulti_pool_mi_fits(fits),
    "support only Wald interval methods"
  )
})

test_that("qgcompmulti_mi validation enforces keep_fits and fit_seeds contracts", {
  fits <- fit_test_model_list()
  pooled <- qgcompmulti_pool_mi_fits(
    fits = fits,
    keep_fits = TRUE,
    fit_seeds = c(1101, 1202, 1303)
  )

  bad_fits <- pooled
  bad_fits$mi_info$keep_fits <- FALSE
  expect_error(
    validate_qgcompmulti_mi(bad_fits),
    "must be `NULL` when `mi_info\\$keep_fits` is `FALSE`"
  )

  bad_seeds <- pooled
  bad_seeds$mi_info$fit_seeds <- c(1L, 2L)
  expect_error(
    validate_qgcompmulti_mi(bad_seeds),
    "one entry per imputation"
  )
})
