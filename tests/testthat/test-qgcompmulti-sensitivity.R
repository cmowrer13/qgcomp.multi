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
  expect_error(
    mcsize_sensitivity(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      MCsize_values = c(nrow(dat) + 1),
      B = 10
    ),
    "less than or equal to `nrow\\(data\\)`"
  )
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
