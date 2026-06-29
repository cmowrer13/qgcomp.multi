test_that("qgcomp.glm.multi.mi pools completed-list fits and stores master-seed metadata", {
  completed <- make_completed_data_list(m = 3, seed = 321)

  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 10,
    seed = 812,
    keep_fits = FALSE
  )

  manual_fits <- lapply(
    seq_along(completed),
    function(i) {
      qgcomp.glm.multi(
        f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
        data = completed[[i]],
        mix1 = c("X1", "X2", "X3"),
        mix2 = c("W1", "W2", "W3"),
        interaction = TRUE,
        q = 4,
        B = 10,
        MCsize = nrow(completed[[i]]),
        seed = fit$mi_info$fit_seeds[[i]]
      )
    }
  )

  manual <- qgcompmulti_pool_mi_fits(
    fits = manual_fits,
    input_type = "completed_list",
    keep_fits = FALSE,
    seed = 812L,
    fit_seeds = fit$mi_info$fit_seeds,
    formula = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C
  )

  expect_s3_class(fit, "qgcompmulti_mi")
  expect_identical(fit$mi_info$input_type, "completed_list")
  expect_false(fit$mi_info$keep_fits)
  expect_equal(fit$analysis$MCsize, nrow(completed[[1]]))
  expect_null(fit$fits$imputation_fits)
  expect_length(fit$mi_info$fit_seeds, 3)
  expect_identical(fit$mi_info$seed, 812L)
  expect_equal(fit$results$coefficients, manual$results$coefficients)
  expect_equal(fit$results$vcov, manual$results$vcov)
  expect_equal(fit$results$df, manual$results$df)
})

test_that("qgcomp.glm.multi.mi optionally retains per-imputation fits and handles q = NULL", {
  completed <- make_completed_data_list(m = 3, seed = 654)

  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = FALSE,
    q = NULL,
    centering = "median",
    B = 10,
    keep_fits = TRUE,
    seed = 415
  )

  expect_s3_class(fit, "qgcompmulti_mi")
  expect_true(fit$mi_info$keep_fits)
  expect_identical(names(fit$results$coefficients), EXPECTED_COEF_NAMES_NO_INTERACTION)
  expect_identical(fit$mixtures$centering, "median")
  expect_length(fit$fits$imputation_fits, 3)
  expect_true(all(vapply(fit$fits$imputation_fits, inherits, logical(1), what = "qgcompmulti")))
})

test_that("qgcomp.glm.multi.mi supports bootstrap-level parallel execution", {
  completed <- make_completed_data_list(m = 3, seed = 321)

  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 6,
    seed = 812,
    keep_fits = FALSE,
    parallel = TRUE,
    workers = 2
  )

  manual_fits <- lapply(
    seq_along(completed),
    function(i) {
      qgcomp.glm.multi(
        f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
        data = completed[[i]],
        mix1 = c("X1", "X2", "X3"),
        mix2 = c("W1", "W2", "W3"),
        interaction = TRUE,
        q = 4,
        B = 6,
        MCsize = nrow(completed[[i]]),
        seed = fit$mi_info$fit_seeds[[i]],
        parallel = TRUE,
        workers = 2
      )
    }
  )

  manual <- qgcompmulti_pool_mi_fits(
    fits = manual_fits,
    input_type = "completed_list",
    keep_fits = FALSE,
    seed = 812L,
    fit_seeds = fit$mi_info$fit_seeds,
    formula = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C
  )

  expect_s3_class(fit, "qgcompmulti_mi")
  expect_null(fit$fits$imputation_fits)
  expect_equal(fit$results$coefficients, manual$results$coefficients)
  expect_equal(fit$results$vcov, manual$results$vcov)
  expect_equal(fit$results$df, manual$results$df)
})

test_that("qgcomp.glm.multi.mi parallel workflows are reproducible within a fixed mode", {
  completed <- make_completed_data_list(m = 3, seed = 321)

  fit1 <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 6,
    seed = 811,
    parallel = TRUE,
    workers = 2
  )

  fit2 <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 6,
    seed = 811,
    parallel = TRUE,
    workers = 2
  )

  expect_equal(fit1$results$coefficients, fit2$results$coefficients)
  expect_equal(fit1$results$vcov, fit2$results$vcov)
  expect_equal(fit1$results$df, fit2$results$df)
  expect_equal(fit1$mi_info$fit_seeds, fit2$mi_info$fit_seeds)
})

test_that("qgcomp.glm.multi.mi parallel workflows handle q = NULL and retained fits", {
  completed <- make_completed_data_list(m = 3, seed = 654)

  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = FALSE,
    q = NULL,
    centering = "median",
    B = 6,
    keep_fits = TRUE,
    seed = 415,
    parallel = TRUE,
    workers = 2
  )

  expect_s3_class(fit, "qgcompmulti_mi")
  expect_true(fit$mi_info$keep_fits)
  expect_identical(names(fit$results$coefficients), EXPECTED_COEF_NAMES_NO_INTERACTION)
  expect_identical(fit$mixtures$centering, "median")
  expect_length(fit$fits$imputation_fits, 3)
  expect_true(all(vapply(fit$fits$imputation_fits, inherits, logical(1), what = "qgcompmulti")))
})

test_that("qgcomp.glm.multi.mi supports clustered completed-list workflows", {
  completed <- make_completed_data_list(m = 3, seed = 777, clustered = TRUE)

  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 10,
    id = "cluster_id",
    seed = 511
  )

  expect_true(fit$data_info$has_clusters)
  expect_identical(fit$data_info$cluster_var, "cluster_id")
  expect_equal(fit$data_info$n_clusters, length(unique(completed[[1]]$cluster_id)))
})

test_that("qgcomp.glm.multi.mi supports clustered parallel completed-list workflows", {
  completed <- make_completed_data_list(m = 3, seed = 777, clustered = TRUE)

  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 6,
    id = "cluster_id",
    seed = 511,
    parallel = TRUE,
    workers = 2
  )

  expect_true(fit$data_info$has_clusters)
  expect_identical(fit$data_info$cluster_var, "cluster_id")
  expect_equal(fit$data_info$n_clusters, length(unique(completed[[1]]$cluster_id)))
})

test_that("qgcomp.glm.multi.mi forwards explicit estimand_scale into completed-data fits", {
  completed <- replicate(3, make_binomial_test_data(seed = 432), simplify = FALSE)
  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    family = binomial(link = "logit"),
    estimand_scale = "risk_difference",
    q = 4,
    B = 5,
    seed = 901,
    keep_fits = TRUE
  )
  expect_identical(fit$analysis$estimand_scale, "risk_difference")
  expect_identical(fit$analysis$msm_fitting_scale, "identity")
  expect_length(fit$fits$imputation_fits, 3)
  expect_true(all(vapply(
    fit$fits$imputation_fits,
    function(x) identical(x$analysis$estimand_scale, "risk_difference"),
    logical(1)
  )))
})
test_that("qgcomp.glm.multi.mi preserves explicit estimand_scale under parallel ordinary fits", {
  completed <- replicate(3, make_poisson_test_data(seed = 654), simplify = FALSE)
  serial_fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    family = poisson(link = "log"),
    estimand_scale = "rate_ratio",
    q = 4,
    B = 5,
    seed = 902,
    keep_fits = FALSE,
    parallel = FALSE
  )
  parallel_fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    family = poisson(link = "log"),
    estimand_scale = "rate_ratio",
    q = 4,
    B = 5,
    seed = 902,
    keep_fits = FALSE,
    parallel = TRUE,
    workers = 2
  )
  expect_identical(parallel_fit$analysis$estimand_scale, "rate_ratio")
  expect_identical(parallel_fit$analysis$msm_fitting_scale, "log")
  expect_equal(serial_fit$results$coefficients, parallel_fit$results$coefficients)
  expect_equal(serial_fit$results$vcov, parallel_fit$results$vcov)
  expect_equal(serial_fit$results$df, parallel_fit$results$df)
})

test_that("qgcomp.glm.multi.mi supports mids input with keep_fits and bootstrap-level parallelism when mice is installed", {
  skip_if_not_installed("mice")

  mids <- make_test_mids(m = 3, seed = 888)
  completed <- mice::complete(mids, action = "all")

  fit_mids <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = mids,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 6,
    seed = 903,
    keep_fits = TRUE,
    parallel = TRUE,
    workers = 2
  )

  fit_list <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 6,
    seed = 903,
    keep_fits = TRUE,
    parallel = TRUE,
    workers = 2
  )

  expect_identical(fit_mids$mi_info$input_type, "mids")
  expect_identical(fit_list$mi_info$input_type, "completed_list")
  expect_true(fit_mids$mi_info$keep_fits)
  expect_identical(fit_mids$mixtures$centering, "median")
  expect_length(fit_mids$fits$imputation_fits, 3L)
  expect_equal(fit_mids$results$coefficients, fit_list$results$coefficients)
  expect_equal(fit_mids$results$vcov, fit_list$results$vcov)
  expect_equal(fit_mids$mi_info$fit_seeds, fit_list$mi_info$fit_seeds)
})

test_that("qgcomp.glm.multi.mi supports mids input when mice is installed", {
  skip_if_not_installed("mice")

  mids <- make_test_mids(m = 3, seed = 888)
  completed <- mice::complete(mids, action = "all")

  fit_mids <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = mids,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 10,
    seed = 902
  )

  fit_list <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = completed,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 10,
    seed = 902
  )

  expect_identical(fit_mids$mi_info$input_type, "mids")
  expect_identical(fit_list$mi_info$input_type, "completed_list")
  expect_equal(fit_mids$results$coefficients, fit_list$results$coefficients)
  expect_equal(fit_mids$results$vcov, fit_list$results$vcov)
  expect_equal(fit_mids$mi_info$fit_seeds, fit_list$mi_info$fit_seeds)
})

test_that("qgcomp.glm.multi.mi parallel mode respects an active future plan when workers is NULL", {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::multisession, workers = 1)

  fit <- qgcomp.glm.multi.mi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = make_completed_data_list(m = 3, seed = 456),
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 4,
    seed = 123,
    parallel = TRUE
  )

  expect_s3_class(fit, "qgcompmulti_mi")
})

test_that("qgcomp.glm.multi.mi parallel mode warns once for progress and errors on unsupported combinations", {
  completed <- make_completed_data_list(m = 3, seed = 123)

  warnings <- character()
  fit <- withCallingHandlers(
    qgcomp.glm.multi.mi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = completed,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = TRUE,
      q = 4,
      B = 4,
      seed = 111,
      progress = TRUE,
      parallel = TRUE,
      workers = 2
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(fit, "qgcompmulti_mi")
  expect_length(warnings, 1L)
  expect_match(warnings[[1]], "progress = TRUE", fixed = TRUE)

  expect_error(
    qgcomp.glm.multi.mi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = completed,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = TRUE,
      q = 4,
      B = 4,
      seed = 222,
      parallel = FALSE,
      workers = 2
    ),
    "`workers` is only supported when `parallel = TRUE`"
  )

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 1)

  expect_error(
    qgcomp.glm.multi.mi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = completed,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = TRUE,
      q = 4,
      B = 4,
      seed = 223,
      parallel = TRUE,
      workers = 2
    ),
    "non-sequential future plan is already active"
  )
})

test_that("qgcomp.glm.multi.mi rejects incomplete or incompatible completed-list inputs", {
  completed <- make_completed_data_list(m = 3, seed = 999)

  bad_missing <- completed
  bad_missing[[2]]$X1[1] <- NA_real_
  expect_error(
    qgcomp.glm.multi.mi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = bad_missing,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = TRUE,
      q = 4,
      B = 10
    ),
    "still contains missing values"
  )

  bad_rows <- completed
  bad_rows[[3]] <- bad_rows[[3]][-1, , drop = FALSE]
  expect_error(
    qgcomp.glm.multi.mi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = bad_rows,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = TRUE,
      q = 4,
      B = 10
    ),
    "same number of rows"
  )

  bad_cluster <- make_completed_data_list(m = 3, seed = 111, clustered = TRUE)
  bad_cluster[[2]] <- bad_cluster[[2]][c(2:nrow(bad_cluster[[2]]), 1), , drop = FALSE]
  rownames(bad_cluster[[2]]) <- NULL
  expect_error(
    qgcomp.glm.multi.mi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = bad_cluster,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = TRUE,
      q = 4,
      B = 10,
      id = "cluster_id"
    ),
    "identical cluster IDs and ordering"
  )
})
