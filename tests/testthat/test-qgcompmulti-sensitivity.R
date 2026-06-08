test_that("mcsize_sensitivity() returns a structured sensitivity object", {
  dat <- make_test_data()
  sens <- mcsize_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    MCsize_values = c(60, 120),
    interaction = TRUE,
    q = 4,
    B = 10,
    seed = 123
  )
  expect_s3_class(sens, "qgcompmulti_mcsize_sensitivity")
  expect_equal(sens$MCsize_values, c(60L, 120L))
  expect_equal(sens$results_table$MCsize, c(60L, 120L))
  expect_identical(sens$settings_fixed$q, 4)
  expect_length(sens$fits, 2L)
})
test_that("mcsize_sensitivity() can be used with q = NULL fits", {
  dat <- make_test_data()
  sens <- mcsize_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    MCsize_values = c(60, 120),
    interaction = TRUE,
    q = NULL,
    centering = "median",
    B = 10,
    seed = 123
  )
  expect_s3_class(sens, "qgcompmulti_mcsize_sensitivity")
  expect_identical(sens$settings_fixed$q, NULL)
  expect_identical(sens$settings_fixed$centering, "median")
})
test_that("q_sensitivity() returns a structured sensitivity object", {
  dat <- make_test_data()
  sens <- q_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    q_values = c(3, 4),
    interaction = TRUE,
    B = 10,
    MCsize = nrow(dat),
    seed = 123
  )
  expect_s3_class(sens, "qgcompmulti_q_sensitivity")
  expect_equal(sens$q_values, c(3L, 4L))
  expect_equal(sens$results_table$q, c(3L, 4L))
  expect_true(is.character(sens$comparability_note))
  expect_match(sens$comparability_note, "not directly comparable across different choices of q")
  expect_length(sens$fits, 2L)
  expect_identical(sens$settings_fixed$MCsize, nrow(dat))
})
test_that("sensitivity helpers can drop stored fits", {
  dat <- make_test_data()
  mc_sens <- mcsize_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    MCsize_values = c(120),
    B = 10,
    keep_fits = FALSE
  )
  q_sens <- q_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    q_values = c(4),
    B = 10,
    keep_fits = FALSE
  )
  expect_null(mc_sens$fits)
  expect_null(q_sens$fits)
})
test_that("sensitivity helpers validate requested values", {
  dat <- make_test_data()
  expect_error(
    mcsize_sensitivity(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      MCsize_values = c(0, 10),
      B = 10
    ),
    "positive integers"
  )
  sens_mc_cap <- mcsize_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    MCsize_values = c(nrow(dat), nrow(dat) + 1, nrow(dat) + 10),
    B = 10,
    seed = 123
  )
  expect_equal(sens_mc_cap$MCsize_values, nrow(dat))
  expect_equal(sens_mc_cap$results_table$MCsize, nrow(dat))
  expect_error(
    q_sensitivity(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      q_values = c(1, 4),
      B = 10
    ),
    "greater than or equal to 2"
  )
})
test_that("sensitivity print methods run and q sensitivity warns about comparability", {
  dat <- make_test_data()
  mc_sens <- mcsize_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    MCsize_values = c(120),
    B = 10
  )
  q_sens <- q_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    q_values = c(4),
    B = 10
  )
  expect_output(print(mc_sens), "MCsize sensitivity")
  expect_output(print(q_sens), "Comparability note")
  expect_output(print(q_sens), "not directly comparable")
})
test_that("q_sensitivity() works for no-interaction fits", {
  dat <- make_test_data()
  sens <- q_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    q_values = c(3, 4),
    interaction = FALSE,
    B = 10,
    seed = 123
  )
  expect_s3_class(sens, "qgcompmulti_q_sensitivity")
  expect_false("psi12" %in% names(sens$results_table))
})
test_that("b_sensitivity() returns a structured sensitivity object", {
  dat <- make_test_data()
  sens <- b_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    B_values = c(6, 10),
    interaction = TRUE,
    q = 4,
    MCsize = nrow(dat),
    seed = 123
  )
  expect_s3_class(sens, "qgcompmulti_b_sensitivity")
  expect_equal(sens$B_values, c(6L, 10L))
  expect_equal(sens$results_table$B, c(6L, 10L))
  expect_identical(sens$settings_fixed$MCsize, nrow(dat))
  expect_length(sens$fits, 2L)

  small_draws <- sens$fits[[1]]$bootstrap$coef_draws
  large_draws_prefix <- sens$fits[[2]]$bootstrap$coef_draws[
    seq_len(nrow(small_draws)),
    ,
    drop = FALSE
  ]
  expect_false(isTRUE(all.equal(
    small_draws,
    large_draws_prefix,
    check.attributes = FALSE
  )))
})

test_that("b_sensitivity() can be used with q = NULL fits", {
  dat <- make_test_data()
  sens <- b_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    B_values = c(6, 10),
    interaction = TRUE,
    q = NULL,
    centering = "median",
    MCsize = nrow(dat),
    seed = 123
  )
  expect_s3_class(sens, "qgcompmulti_b_sensitivity")
  expect_identical(sens$settings_fixed$q, NULL)
  expect_identical(sens$settings_fixed$centering, "median")
})

test_that("b_sensitivity() is reproducible when seed is supplied", {
  dat <- make_test_data()
  sens1 <- b_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    B_values = c(6, 10),
    interaction = TRUE,
    q = 4,
    MCsize = 40,
    seed = 2024
  )
  sens2 <- b_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    B_values = c(6, 10),
    interaction = TRUE,
    q = 4,
    MCsize = 40,
    seed = 2024
  )
  expect_equal(sens1$results_table, sens2$results_table)
  expect_equal(sens1$fits[[1]]$bootstrap$coef_draws, sens2$fits[[1]]$bootstrap$coef_draws)
  expect_equal(sens1$fits[[2]]$bootstrap$coef_draws, sens2$fits[[2]]$bootstrap$coef_draws)
})

test_that("b_sensitivity() can drop stored fits", {
  dat <- make_test_data()
  b_sens <- b_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    B_values = c(10),
    MCsize = nrow(dat),
    keep_fits = FALSE
  )
  expect_null(b_sens$fits)
})

test_that("b_sensitivity() validates requested values", {
  dat <- make_test_data()
  expect_error(
    b_sensitivity(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      B_values = c(1, 10),
      MCsize = nrow(dat)
    ),
    "greater than or equal to 2"
  )

  sens_b_unique <- b_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    B_values = c(10, 10, 20),
    MCsize = nrow(dat),
    seed = 123
  )
  expect_equal(sens_b_unique$B_values, c(10L, 20L))
  expect_equal(sens_b_unique$results_table$B, c(10L, 20L))
})

test_that("b_sensitivity print method runs", {
  dat <- make_test_data()
  b_sens <- b_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    B_values = c(10),
    MCsize = nrow(dat)
  )
  expect_output(print(b_sens), "B sensitivity")
})

test_that("b_sensitivity() works for no-interaction fits", {
  dat <- make_test_data()
  sens <- b_sensitivity(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    B_values = c(6, 10),
    interaction = FALSE,
    q = 4,
    MCsize = nrow(dat),
    seed = 123
  )
  expect_s3_class(sens, "qgcompmulti_b_sensitivity")
  expect_false("psi12" %in% names(sens$results_table))
})
