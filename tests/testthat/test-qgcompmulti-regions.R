with_region_plot_device <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  force(code)
}
test_that("confregion() builds the default two-parameter region", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  region <- confregion(fit)
  expect_s3_class(region, "qgcompmulti_region")
  expect_identical(region$parm, c("psi1", "psi2"))
  expect_equal(region$center, coef(fit)[region$parm])
  expect_equal(region$covariance, vcov(fit)[region$parm, region$parm])
  expect_equal(region$threshold, stats::qchisq(0.95, df = 2))
  expect_identical(region$df, 2L)
  expect_identical(region$method, "bootstrap_chisq")
  expect_identical(region$estimand_scale, fit$analysis$estimand_scale)
  expect_identical(region$msm_fitting_scale, fit$analysis$msm_fitting_scale)
  expect_s3_class(region$plot_data, "data.frame")
  expect_named(region$plot_data, c("x", "y"))
  expect_true(all(is.finite(region$plot_data$x)))
  expect_true(all(is.finite(region$plot_data$y)))
  centered_boundary <- sweep(as.matrix(region$plot_data), 2L, region$center, FUN = "-")
  quadratic_form <- rowSums((centered_boundary %*% solve(region$covariance)) * centered_boundary)
  expect_equal(quadratic_form, rep(region$threshold, length(quadratic_form)), tolerance = 1e-8)
})
test_that("confregion() supports named, positional, and all-parameter selections", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  by_name <- confregion(fit, parm = c("psi1", "psi1:psi2"), level = 0.90)
  by_position <- confregion(fit, parm = c(2, 4), level = 0.90)
  all_parameters <- confregion(fit, parm = NULL, level = 0.90)
  expect_identical(by_name$parm, c("psi1", "psi1:psi2"))
  expect_equal(by_name$center, by_position$center)
  expect_equal(by_name$covariance, by_position$covariance)
  expect_equal(by_name$threshold, stats::qchisq(0.90, df = 2))
  expect_identical(all_parameters$parm, names(coef(fit)))
  expect_identical(all_parameters$df, length(coef(fit)))
  expect_null(all_parameters$plot_data)
})
test_that("confregion() validates region inputs clearly", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  expect_error(confregion(fit, parm = "psi1"), "at least two")
  expect_error(confregion(fit, parm = "not_a_parameter"), "Unknown coefficient")
  expect_error(confregion(fit, parm = 99), "valid coefficient positions")
  expect_error(confregion(fit, level = 1), "strictly between 0 and 1")
  expect_error(confregion(fit, method = "wald"), "method")
  expect_error(confregion(fit, npoints = 3), "npoints")
})
test_that("region helpers require at least two named coefficients", {
  fit <- fit_test_model(interaction = TRUE, q = 4)
  expect_identical(
    qgcompmulti_resolve_region_parm(c("psi1", "psi2"), names(coef(fit))),
    c("psi1", "psi2")
  )
  expect_error(
    qgcompmulti_resolve_region_parm("psi1", names(coef(fit))),
    "at least two"
  )
  covariance <- diag(c(0.04, 0.03), nrow = 2L)
  rownames(covariance) <- c("psi1", "psi2")
  colnames(covariance) <- c("psi1", "psi2")
  region <- new_qgcompmulti_region(
    parm = c("psi1", "psi2"),
    center = c(psi1 = 0.2, psi2 = 0.1),
    covariance = covariance
  )
  expect_s3_class(region, "qgcompmulti_region")
})
test_that("confidence-region objects print scientific scale metadata", {
  fit <- fit_test_model(
    interaction = TRUE,
    q = 4,
    B = 10,
    family = binomial(link = "logit")
  )
  region <- confregion(fit)
  printed <- capture.output(print(region))
  expect_identical(region$estimand_scale, "odds_ratio")
  expect_identical(region$msm_fitting_scale, "logit")
  expect_equal(region$center, coef(fit)[region$parm])
  expect_match(paste(printed, collapse = "\n"), "log coefficient scale")
  expect_match(paste(printed, collapse = "\n"), "not geometrically exponentiated")
})
test_that("plot.qgcompmulti_region() draws two-parameter regions only", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  region <- confregion(fit)
  higher_dimensional <- confregion(fit, parm = NULL)
  expect_no_error(with_region_plot_device(plot(region)))
  expect_error(
    with_region_plot_device(plot(higher_dimensional)),
    "two-parameter confidence regions"
  )
})
test_that("pooled MI confidence regions remain unsupported", {
  fit <- fit_test_model(interaction = TRUE, q = 4, B = 10)
  mi_like <- fit
  class(mi_like) <- "qgcompmulti_mi"
  expect_error(
    confregion(mi_like),
    "not supported for pooled qgcompmulti multiple-imputation objects"
  )
})
