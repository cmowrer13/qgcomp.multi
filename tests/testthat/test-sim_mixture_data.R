test_that("sim_mixture_data returns expected columns and dimensions", {
  dat <- sim_mixture_data(
    n = 50,
    pA = 2,
    pB = 2,
    rho_within_A = 0.3,
    rho_within_B = 0.3,
    rho_between = 0.2,
    psi1 = 0.5,
    psi2 = 0.3,
    psi12 = 0.2,
    return_quantized = FALSE,
    seed = 123
  )
  expect_s3_class(dat, "data.frame")
  expect_identical(nrow(dat), 50L)
  expect_true(all(c("Y", "X1", "X2", "W1", "W2", "C") %in% names(dat)))
})

test_that("sim_mixture_data is reproducible with a fixed seed", {
  dat1 <- sim_mixture_data(
    n = 50,
    pA = 2,
    pB = 2,
    rho_within_A = 0.3,
    rho_within_B = 0.3,
    rho_between = 0.2,
    psi1 = 0.5,
    psi2 = 0.3,
    psi12 = 0.2,
    return_quantized = FALSE,
    seed = 999
  )
  dat2 <- sim_mixture_data(
    n = 50,
    pA = 2,
    pB = 2,
    rho_within_A = 0.3,
    rho_within_B = 0.3,
    rho_between = 0.2,
    psi1 = 0.5,
    psi2 = 0.3,
    psi12 = 0.2,
    return_quantized = FALSE,
    seed = 999
  )
  expect_equal(dat1, dat2)
})

test_that("sim_mixture_data can return quantized exposures", {
  dat_q <- sim_mixture_data(
    n = 50,
    pA = 2,
    pB = 2,
    rho_within_A = 0.3,
    rho_within_B = 0.3,
    rho_between = 0.2,
    psi1 = 0.5,
    psi2 = 0.3,
    psi12 = 0.2,
    q = 4,
    return_quantized = TRUE,
    seed = 123
  )
  exposure_values <- unlist(dat_q[c("X1", "X2", "W1", "W2")], use.names = FALSE)
  expect_true(all(exposure_values %in% 0:3))
})

test_that("sim_mixture_data adapts output dimensions to pA and pB", {
  dat <- sim_mixture_data(
    n = 30,
    pA = 3,
    pB = 2,
    rho_within_A = 0.3,
    rho_within_B = 0.3,
    rho_between = 0.2,
    psi1 = 0.5,
    psi2 = 0.3,
    psi12 = 0.2,
    return_quantized = FALSE,
    seed = 321
  )
  expect_true(all(c("X1", "X2", "X3", "W1", "W2") %in% names(dat)))
  expect_false("W3" %in% names(dat))
})
