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

test_that("sim_mixture_data preserves legacy Gaussian defaults", {
  dat_default <- sim_mixture_data(
    n = 80,
    pA = 2,
    pB = 2,
    rho_within_A = 0.3,
    rho_within_B = 0.3,
    rho_between = 0.2,
    psi1 = 0.5,
    psi2 = 0.3,
    psi12 = 0.2,
    seed = 321
  )

  dat_explicit <- sim_mixture_data(
    n = 80,
    pA = 2,
    pB = 2,
    rho_within_A = 0.3,
    rho_within_B = 0.3,
    rho_between = 0.2,
    psi1 = 0.5,
    psi2 = 0.3,
    psi12 = 0.2,
    seed = 321,
    family = gaussian(),
    estimand_scale = "mean_difference",
    baseline_mean = 0,
    beta_C = 1
  )

  expect_equal(dat_default, dat_explicit)
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

test_that("sim_mixture_data generates binary outcomes with natural-scale odds ratios", {
  dat <- sim_mixture_data(
    n = 80000,
    pA = 1,
    pB = 1,
    rho_within_A = 0,
    rho_within_B = 0,
    rho_between = 0,
    psi1 = 2,
    psi2 = 1,
    psi12 = 1,
    q = 4,
    return_quantized = TRUE,
    seed = 42,
    family = binomial(),
    baseline_risk = 0.2,
    beta_C = 0
  )

  risk0 <- mean(dat$Y[dat$X1 == 0])
  risk1 <- mean(dat$Y[dat$X1 == 1])

  odds_ratio <- (risk1 / (1 - risk1)) / (risk0 / (1 - risk0))

  expect_true(all(dat$Y %in% 0:1))
  expect_equal(odds_ratio, 2, tolerance = 0.2)
})

test_that("sim_mixture_data generates count outcomes with natural-scale rate ratios", {
  dat <- sim_mixture_data(
    n = 80000,
    pA = 1,
    pB = 1,
    rho_within_A = 0,
    rho_within_B = 0,
    rho_between = 0,
    psi1 = 1.5,
    psi2 = 1,
    psi12 = 1,
    q = 4,
    return_quantized = TRUE,
    seed = 84,
    family = poisson(),
    baseline_rate = 1,
    beta_C = 0
  )

  rate0 <- mean(dat$Y[dat$X1 == 0])
  rate1 <- mean(dat$Y[dat$X1 == 1])

  expect_true(all(dat$Y >= 0))
  expect_true(all(dat$Y == floor(dat$Y)))
  expect_equal(rate1 / rate0, 1.5, tolerance = 0.12)
})

test_that("sim_mixture_data supports additive binomial generation when probabilities are valid", {
  dat <- sim_mixture_data(
    n = 80000,
    pA = 1,
    pB = 1,
    rho_within_A = 0,
    rho_within_B = 0,
    rho_between = 0,
    psi1 = 0.05,
    psi2 = 0,
    psi12 = 0,
    q = 4,
    return_quantized = TRUE,
    seed = 21,
    family = binomial(),
    estimand_scale = "risk_difference",
    baseline_risk = 0.2,
    beta_C = 0
  )

  risk0 <- mean(dat$Y[dat$X1 == 0])
  risk1 <- mean(dat$Y[dat$X1 == 1])

  expect_lt(abs((risk1 - risk0) - 0.05), 0.005)
})

test_that("sim_mixture_data supports additive Poisson generation when means are valid", {
  dat <- sim_mixture_data(
    n = 80000,
    pA = 1,
    pB = 1,
    rho_within_A = 0,
    rho_within_B = 0,
    rho_between = 0,
    psi1 = 0.75,
    psi2 = 0,
    psi12 = 0,
    q = 4,
    return_quantized = TRUE,
    seed = 22,
    family = poisson(),
    estimand_scale = "mean_difference",
    baseline_rate = 2,
    beta_C = 0
  )

  rate0 <- mean(dat$Y[dat$X1 == 0])
  rate1 <- mean(dat$Y[dat$X1 == 1])

  expect_equal(rate1 - rate0, 0.75, tolerance = 0.08)
})

test_that("sim_mixture_data rejects invalid baseline and ratio inputs", {
  expect_error(
    sim_mixture_data(
      n = 50,
      pA = 2,
      pB = 2,
      rho_within_A = 0.3,
      rho_within_B = 0.3,
      rho_between = 0.2,
      psi1 = 1.5,
      psi2 = 1.2,
      psi12 = 1.1,
      family = binomial(),
      baseline_risk = 1.1
    ),
    "baseline_risk"
  )

  expect_error(
    sim_mixture_data(
      n = 50,
      pA = 2,
      pB = 2,
      rho_within_A = 0.3,
      rho_within_B = 0.3,
      rho_between = 0.2,
      psi1 = 1.5,
      psi2 = 1.2,
      psi12 = 1.1,
      family = poisson(),
      baseline_rate = 0
    ),
    "baseline_rate"
  )

  expect_error(
    sim_mixture_data(
      n = 50,
      pA = 2,
      pB = 2,
      rho_within_A = 0.3,
      rho_within_B = 0.3,
      rho_between = 0.2,
      psi1 = 0,
      psi2 = 1.2,
      psi12 = 1.1,
      family = poisson()
    ),
    "strictly positive"
  )
})

test_that("sim_mixture_data rejects impossible additive binomial and Poisson requests", {
  expect_error(
    sim_mixture_data(
      n = 100,
      pA = 1,
      pB = 1,
      rho_within_A = 0,
      rho_within_B = 0,
      rho_between = 0,
      psi1 = 0.5,
      psi2 = 0,
      psi12 = 0,
      family = binomial(),
      estimand_scale = "risk_difference",
      baseline_risk = 0.8,
      beta_C = 0
    ),
    "outside \\[0, 1\\]"
  )

  expect_error(
    sim_mixture_data(
      n = 100,
      pA = 1,
      pB = 1,
      rho_within_A = 0,
      rho_within_B = 0,
      rho_between = 0,
      psi1 = -2,
      psi2 = 0,
      psi12 = 0,
      family = poisson(),
      estimand_scale = "mean_difference",
      baseline_rate = 1,
      beta_C = 0
    ),
    "nonpositive means"
  )
})

test_that("sim_mixture_data rejects unsupported families and invalid correlation structures", {
  expect_error(
    sim_mixture_data(
      n = 50,
      pA = 2,
      pB = 2,
      rho_within_A = 0.3,
      rho_within_B = 0.3,
      rho_between = 0.2,
      psi1 = 0.5,
      psi2 = 0.3,
      psi12 = 0.2,
      family = quasipoisson()
    ),
    "currently supports only"
  )

  expect_error(
    sim_mixture_data(
      n = 50,
      pA = 3,
      pB = 1,
      rho_within_A = -0.8,
      rho_within_B = 0.3,
      rho_between = 0,
      psi1 = 0.5,
      psi2 = 0.3,
      psi12 = 0.2
    ),
    "positive-definite"
  )
})
