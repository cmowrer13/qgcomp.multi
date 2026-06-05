test_that("sim_biomarker_data returns expected columns and positive biomarkers", {
  dat <- sim_biomarker_data(n = 120, seed = 123)

  expected_names <- c(
    "inflammation_score",
    "mep_ng_ml", "mibp_ng_ml", "mehpp_ng_ml",
    "bpa_ng_ml", "bp3_ng_ml", "ppb_ng_ml",
    "age_years", "bmi_kg_m2", "smoker"
  )

  expect_s3_class(dat, "data.frame")
  expect_identical(nrow(dat), 120L)
  expect_true(all(expected_names %in% names(dat)))
  expect_true(all(dat[c("mep_ng_ml", "mibp_ng_ml", "mehpp_ng_ml", "bpa_ng_ml", "bp3_ng_ml", "ppb_ng_ml")] > 0))
  expect_true(all(dat$age_years >= 20 & dat$age_years <= 75))
  expect_true(all(dat$bmi_kg_m2 >= 18 & dat$bmi_kg_m2 <= 45))
  expect_true(all(dat$smoker %in% 0:1))
})

test_that("sim_biomarker_data can include log-transformed biomarkers", {
  dat <- sim_biomarker_data(n = 80, include_log = TRUE, seed = 456)

  expect_true(all(c("ln_mep", "ln_mibp", "ln_mehpp", "ln_bpa", "ln_bp3", "ln_ppb") %in% names(dat)))
  expect_equal(dat$ln_mep, log(dat$mep_ng_ml))
  expect_equal(dat$ln_bpa, log(dat$bpa_ng_ml))
})

test_that("sim_biomarker_data is reproducible and restores RNG state", {
  set.seed(999)
  reference <- stats::runif(2)

  set.seed(999)
  dat1 <- sim_biomarker_data(n = 60, seed = 2024)
  after1 <- stats::runif(1)

  set.seed(999)
  dat2 <- sim_biomarker_data(n = 60, seed = 2024)
  after2 <- stats::runif(1)

  expect_equal(dat1, dat2)
  expect_equal(after1, after2)
  expect_equal(after1, reference[1])
})

test_that("sim_biomarker_data has right-skewed biomarkers and stronger within-mixture correlation", {
  dat <- sim_biomarker_data(n = 600, include_log = TRUE, seed = 321)

  raw_vars <- c("mep_ng_ml", "mibp_ng_ml", "mehpp_ng_ml", "bpa_ng_ml", "bp3_ng_ml", "ppb_ng_ml")
  skew_flags <- vapply(dat[raw_vars], function(x) mean(x) > stats::median(x), logical(1))
  expect_true(sum(skew_flags) >= 5L)

  corr <- stats::cor(dat[c("ln_mep", "ln_mibp", "ln_mehpp", "ln_bpa", "ln_bp3", "ln_ppb")])
  within_phthalates <- mean(c(corr["ln_mep", "ln_mibp"], corr["ln_mep", "ln_mehpp"], corr["ln_mibp", "ln_mehpp"]))
  within_phenols <- mean(c(corr["ln_bpa", "ln_bp3"], corr["ln_bpa", "ln_ppb"], corr["ln_bp3", "ln_ppb"]))
  between <- mean(corr[1:3, 4:6])

  expect_gt(within_phthalates, between)
  expect_gt(within_phenols, between)
  expect_gt(within_phthalates, 0.30)
  expect_gt(within_phenols, 0.15)
})

test_that("sim_biomarker_data supports quantized and q = NULL workflows", {
  dat <- sim_biomarker_data(n = 180, include_log = TRUE, seed = 111)

  fit_q <- qgcomp.glm.multi(
    f = inflammation_score ~ ln_mep + ln_mibp + ln_mehpp + ln_bpa + ln_bp3 + ln_ppb +
      age_years + bmi_kg_m2 + smoker,
    data = dat,
    mix1 = c("ln_mep", "ln_mibp", "ln_mehpp"),
    mix2 = c("ln_bpa", "ln_bp3", "ln_ppb"),
    interaction = TRUE,
    q = 4,
    B = 4,
    seed = 2024
  )

  fit_cont <- qgcomp.glm.multi(
    f = inflammation_score ~ ln_mep + ln_mibp + ln_mehpp + ln_bpa + ln_bp3 + ln_ppb +
      age_years + bmi_kg_m2 + smoker,
    data = dat,
    mix1 = c("ln_mep", "ln_mibp", "ln_mehpp"),
    mix2 = c("ln_bpa", "ln_bp3", "ln_ppb"),
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 4,
    seed = 2024
  )

  expect_s3_class(fit_q, "qgcompmulti")
  expect_s3_class(fit_cont, "qgcompmulti")

  support_q <- support(fit_q)
  support_cont <- support(fit_cont)
  adequacy_q <- adequacy(fit_q)
  adequacy_cont <- adequacy(fit_cont)

  expect_identical(support_q$flags$is_quantized, TRUE)
  expect_identical(support_cont$flags$is_original_scale, TRUE)
  expect_true(all(is.finite(unlist(adequacy_q$summary_metrics))))
  expect_true(all(is.finite(unlist(adequacy_cont$summary_metrics))))
})

test_that("reference vignette fit shows positive mixture effects and interaction", {
  dat <- sim_biomarker_data(n = 400, include_log = TRUE, seed = 123)

  fit_q <- qgcomp.glm.multi(
    f = inflammation_score ~ ln_mep + I(ln_mep^2) + ln_mibp + ln_mehpp + ln_bpa + ln_bp3 + ln_ppb +
      age_years + bmi_kg_m2 + smoker,
    data = dat,
    mix1 = c("ln_mep", "ln_mibp", "ln_mehpp"),
    mix2 = c("ln_bpa", "ln_bp3", "ln_ppb"),
    interaction = TRUE,
    q = 4,
    B = 4,
    seed = 2024
  )

  fit_cont <- qgcomp.glm.multi(
    f = inflammation_score ~ ln_mep + I(ln_mep^2) + ln_mibp + ln_mehpp + ln_bpa + ln_bp3 + ln_ppb +
      age_years + bmi_kg_m2 + smoker,
    data = dat,
    mix1 = c("ln_mep", "ln_mibp", "ln_mehpp"),
    mix2 = c("ln_bpa", "ln_bp3", "ln_ppb"),
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 4,
    seed = 2024
  )

  expect_gt(unname(coef(fit_q)["psi1"]), 0)
  expect_gt(unname(coef(fit_q)["psi2"]), 0)
  expect_gt(unname(coef(fit_q)["psi1:psi2"]), 0)

  expect_gt(unname(coef(fit_cont)["psi1"]), 0)
  expect_gt(unname(coef(fit_cont)["psi2"]), 0)
  expect_gt(unname(coef(fit_cont)["psi1:psi2"]), 0)
})

test_that("reference vignette fit with a mild nonlinear term has informative adequacy error", {
  dat <- sim_biomarker_data(n = 400, include_log = TRUE, seed = 123)

  fit_q <- qgcomp.glm.multi(
    f = inflammation_score ~ ln_mep + I(ln_mep^2) + ln_mibp + ln_mehpp + ln_bpa + ln_bp3 + ln_ppb +
      age_years + bmi_kg_m2 + smoker,
    data = dat,
    mix1 = c("ln_mep", "ln_mibp", "ln_mehpp"),
    mix2 = c("ln_bpa", "ln_bp3", "ln_ppb"),
    interaction = TRUE,
    q = 4,
    B = 4,
    seed = 2024
  )

  fit_cont <- qgcomp.glm.multi(
    f = inflammation_score ~ ln_mep + I(ln_mep^2) + ln_mibp + ln_mehpp + ln_bpa + ln_bp3 + ln_ppb +
      age_years + bmi_kg_m2 + smoker,
    data = dat,
    mix1 = c("ln_mep", "ln_mibp", "ln_mehpp"),
    mix2 = c("ln_bpa", "ln_bp3", "ln_ppb"),
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 4,
    seed = 2024
  )

  metrics_q <- adequacy(fit_q)$summary_metrics
  metrics_cont <- adequacy(fit_cont)$summary_metrics

  expect_gt(metrics_q$mae, 0.005)
  expect_lt(metrics_q$mae, 0.20)
  expect_gt(metrics_q$max_abs_error, 0.005)
  expect_gt(metrics_cont$mae, 0.005)
  expect_lt(metrics_cont$mae, 0.20)
})
