test_that("qgcomp.glm.multi is reproducible when seed is supplied", {
  dat <- make_test_data(seed = 123)
  fit1 <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 10,
    MCsize = 40,
    seed = 2024
  )
  fit2 <- qgcomp.glm.multi(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    B = 10,
    MCsize = 40,
    seed = 2024
  )
  expect_identical(fit1$analysis$seed, 2024L)
  expect_identical(fit2$analysis$seed, 2024L)
  expect_equal(fit1$data_info$n_used, 40)
  expect_equal(fit2$data_info$n_used, 40)
  expect_equal(fit1$results$coefficients, fit2$results$coefficients)
  expect_equal(fit1$results$vcov, fit2$results$vcov)
  expect_equal(fit1$bootstrap$coef_draws, fit2$bootstrap$coef_draws)
})
test_that("qgcompmulti_msm_fit is reproducible when seed is supplied", {
  dat <- make_test_data(seed = 123)
  engine1 <- qgcompmulti_msm_fit(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    MCsize = 40,
    seed = 2024
  )
  engine2 <- qgcompmulti_msm_fit(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    interaction = TRUE,
    q = 4,
    MCsize = 40,
    seed = 2024
  )
  expect_equal(engine1$n_used, 40)
  expect_equal(engine2$n_used, 40)
  expect_equal(engine1$coefficients, engine2$coefficients)
  expect_equal(stats::coef(engine1$msm_fit), stats::coef(engine2$msm_fit))
})
