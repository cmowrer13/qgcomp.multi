test_that("parallel ordinary fits preserve the qgcompmulti object contract", {
  dat <- make_test_data(seed = 123)

  serial_fit <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 6,
    MCsize = 40,
    seed = 2024,
    parallel = FALSE
  )

  parallel_fit <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 6,
    MCsize = 40,
    seed = 2024,
    parallel = TRUE,
    workers = 2
  )

  expect_s3_class(parallel_fit, "qgcompmulti")
  expect_identical(names(parallel_fit), EXPECTED_TOP_COMPONENTS)

  for (component_name in names(EXPECTED_COMPONENT_FIELDS)) {
    expect_identical(
      names(parallel_fit[[component_name]]),
      EXPECTED_COMPONENT_FIELDS[[component_name]]
    )
  }

  expect_equal(serial_fit$results$coefficients, parallel_fit$results$coefficients)
  expect_equal(serial_fit$results$vcov, parallel_fit$results$vcov)
  expect_equal(serial_fit$bootstrap$coef_draws, parallel_fit$bootstrap$coef_draws)
  expect_equal(serial_fit$bootstrap$B_requested, parallel_fit$bootstrap$B_requested)
  expect_equal(serial_fit$bootstrap$B_success, parallel_fit$bootstrap$B_success)
  expect_equal(serial_fit$bootstrap$B_failed, parallel_fit$bootstrap$B_failed)
})

test_that("parallel ordinary fits are reproducible within a fixed mode", {
  dat <- make_test_data(seed = 123)

  fit1 <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 6,
    MCsize = 40,
    seed = 811,
    parallel = TRUE,
    workers = 2
  )

  fit2 <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 6,
    MCsize = 40,
    seed = 811,
    parallel = TRUE,
    workers = 2
  )

  expect_equal(fit1$results$coefficients, fit2$results$coefficients)
  expect_equal(fit1$results$vcov, fit2$results$vcov)
  expect_equal(fit1$bootstrap$coef_draws, fit2$bootstrap$coef_draws)
})

test_that("parallel ordinary fits cover no-interaction and q = NULL branches", {
  dat <- make_test_data(seed = 321)

  no_interaction_fit <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = FALSE,
    q = 4,
    B = 5,
    MCsize = 40,
    seed = 500,
    parallel = TRUE,
    workers = 2
  )

  qnull_fit <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 5,
    MCsize = 40,
    seed = 501,
    parallel = TRUE,
    workers = 2
  )

  expect_identical(names(coef(no_interaction_fit)), EXPECTED_COEF_NAMES_NO_INTERACTION)
  expect_identical(names(coef(qnull_fit)), EXPECTED_COEF_NAMES_WITH_INTERACTION)
  expect_null(qnull_fit$mixtures$q)
})

test_that("parallel ordinary fits support clustered bootstrap", {
  dat <- make_clustered_test_data(seed = 555)

  fit <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 6,
    id = "cluster_id",
    MCsize = nrow(dat),
    seed = 900,
    parallel = TRUE,
    workers = 2
  )

  expect_true(fit$data_info$has_clusters)
  expect_identical(fit$data_info$cluster_var, "cluster_id")
  expect_equal(fit$bootstrap$B_requested, 6L)
  expect_equal(nrow(fit$bootstrap$coef_draws), fit$bootstrap$B_success)
})

test_that("parallel mode disables progress with an explicit warning", {
  dat <- make_test_data(seed = 123)

  expect_warning(
    out <- capture.output(
      fit <- qgcomp.glm.multi(
        f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
        data = dat,
        mix1 = c("X1", "X2", "X3"),
        mix2 = c("W1", "W2", "W3"),
        interaction = TRUE,
        q = 4,
        B = 4,
        MCsize = 40,
        seed = 111,
        progress = TRUE,
        parallel = TRUE,
        workers = 2
      )
    ),
    "progress = TRUE"
  )

  expect_s3_class(fit, "qgcompmulti")
  expect_length(out, 0L)
})

test_that("parallel mode respects an active future plan when workers is NULL", {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::multisession, workers = 1)

  fit <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = make_test_data(seed = 456),
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 4,
    MCsize = 40,
    seed = 123,
    parallel = TRUE
  )

  expect_s3_class(fit, "qgcompmulti")
})

test_that("parallel mode uses the default internal backend when workers is NULL", {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::sequential)

  dat <- make_test_data(seed = 789)

  fit1 <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 4,
    MCsize = 40,
    seed = 333,
    parallel = TRUE,
    workers = NULL
  )

  fit2 <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 4,
    MCsize = 40,
    seed = 333,
    parallel = TRUE,
    workers = NULL
  )

  expect_s3_class(fit1, "qgcompmulti")
  expect_equal(fit1$bootstrap$B_requested, 4L)
  expect_equal(nrow(fit1$bootstrap$coef_draws), fit1$bootstrap$B_success)
  expect_equal(fit1$results$coefficients, fit2$results$coefficients)
  expect_equal(fit1$results$vcov, fit2$results$vcov)
  expect_equal(fit1$bootstrap$coef_draws, fit2$bootstrap$coef_draws)
})

test_that("unsupported parallel combinations error clearly", {
  dat <- make_test_data(seed = 123)

  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = TRUE,
      q = 4,
      B = 4,
      MCsize = 40,
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
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = TRUE,
      q = 4,
      B = 4,
      MCsize = 40,
      seed = 223,
      parallel = TRUE,
      workers = 2
    ),
    "non-sequential future plan is already active"
  )
})
