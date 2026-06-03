test_that("quantize_mixtures preserves dimensions and quantizes mixture columns", {
  dat <- make_test_data()
  out <- quantize_mixtures(
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    q = 4
  )
  expect_s3_class(out, "data.frame")
  expect_identical(dim(out), dim(dat))
  expect_identical(names(out), names(dat))
  mixture_values <- unlist(out[c("X1", "X2", "X3", "W1", "W2", "W3")], use.names = FALSE)
  expect_true(all(mixture_values %in% 0:3))
  expect_equal(out$Y, dat$Y)
  expect_equal(out$C, dat$C)
})

test_that("quantize_mixtures returns unchanged data when q = NULL", {
  dat <- make_test_data()
  out <- quantize_mixtures(
    data = dat,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    q = NULL
  )
  expect_equal(out, dat)
})

test_that("quantize_mixtures validates mixture definitions and data columns", {
  dat <- make_test_data()
  expect_error(
    quantize_mixtures(
      data = dat,
      mix1 = c("X1", "X2"),
      mix2 = c("X2", "W1"),
      q = 4
    )
  )
  expect_error(
    quantize_mixtures(
      data = dat,
      mix1 = c("X1", "X1"),
      mix2 = c("W1", "W2", "W3"),
      q = 4
    )
  )
  expect_error(
    quantize_mixtures(
      data = dat,
      mix1 = c("X1", "X_missing", "X3"),
      mix2 = c("W1", "W2", "W3"),
      q = 4
    )
  )
})

test_that("quantize_mixtures errors when a variable lacks enough unique values", {
  dat_bad <- make_test_data()
  dat_bad$X1 <- 1
  expect_error(
    quantize_mixtures(
      data = dat_bad,
      mix1 = c("X1"),
      mix2 = c("W1"),
      q = 4
    )
  )
})
