test_that("Monte Carlo subsampling works without clustering", {
  fit_mc <- fit_test_model(
    interaction = TRUE,
    q = 4,
    mcsize = 40,
    clustered = FALSE
  )
  expect_s3_class(fit_mc, "qgcompmulti")
  expect_equal(fit_mc$data_info$n_input, nrow(make_test_data()))
  expect_equal(fit_mc$data_info$n_used, 40)
  expect_equal(nobs(fit_mc), 40)
  expect_identical(names(coef(fit_mc)), EXPECTED_COEF_NAMES_WITH_INTERACTION)
})
test_that("Monte Carlo subsampling works with clustered resampling", {
  fit_cluster_mc <- fit_test_model(
    interaction = TRUE,
    q = 4,
    mcsize = 4,
    clustered = TRUE
  )
  expect_s3_class(fit_cluster_mc, "qgcompmulti")
  expect_true(isTRUE(fit_cluster_mc$data_info$has_clusters))
  expect_equal(fit_cluster_mc$data_info$n_clusters, 12)
  expect_equal(fit_cluster_mc$data_info$n_input, nrow(make_clustered_test_data()))
  expect_equal(fit_cluster_mc$data_info$n_used, 40)
  expect_equal(nobs(fit_cluster_mc), 40)
})

test_that("clustered engine treats MCsize above the current row count as full evaluation", {
  dat <- make_uneven_clustered_test_data(
    seed = 123,
    cluster_sizes = c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 7, 8)
  )
  original_n <- nrow(dat)

  split_data <- split(dat, dat$cluster_id)
  sampled_ids <- c("1", "1", "2", "2", "3", "4", "5", "6", "7", "8")
  boot_dat <- do.call(rbind, split_data[sampled_ids])
  row.names(boot_dat) <- NULL

  expect_lt(nrow(boot_dat), original_n)

  fit_engine <- qgcompmulti_msm_fit(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = boot_dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    id = "cluster_id",
    MCsize = original_n,
    seed = 123
  )

  expect_equal(fit_engine$n_used, nrow(boot_dat))
})

test_that("clustered bootstrap does not fail when full-sample MCsize is carried into smaller resamples", {
  dat <- make_uneven_clustered_test_data(
    seed = 123,
    cluster_sizes = c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 7, 8)
  )

  fit <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 20,
    id = "cluster_id",
    MCsize = nrow(dat),
    seed = 123
  )

  expect_s3_class(fit, "qgcompmulti")
  expect_equal(fit$bootstrap$B_failed, 0L)
  expect_equal(fit$bootstrap$B_success, 20L)
})
