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
