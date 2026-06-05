test_that("validation rejects overlapping and duplicated mixture definitions", {
  dat <- make_test_data()
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2"),
      mix2 = c("X2", "W1"),
      B = 10
    )
  )
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X1"),
      mix2 = c("W1", "W2", "W3"),
      B = 10
    )
  )
})

test_that("validation rejects missing variables and missing formula terms", {
  dat <- make_test_data()
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X_missing", "X3"),
      mix2 = c("W1", "W2", "W3"),
      B = 10
    )
  )
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      B = 10
    )
  )
})

test_that("validation rejects invalid q, centering, B, MCsize, and seed values", {
  dat <- make_test_data()
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      q = 1,
      B = 10
    )
  )
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      q = NULL,
      centering = "mean",
      B = 10
    )
  )
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      B = 1
    )
  )
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      B = 10,
      MCsize = 0
    )
  )
  expect_silent(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      B = 10,
      MCsize = nrow(dat) + 1
    )
  )
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      B = 10,
      seed = 3.5
    )
  )
})

test_that("validation rejects invalid cluster ID inputs", {
  dat <- make_test_data()
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      B = 10,
      id = 1
    )
  )
  dat_missing_id <- dat
  dat_missing_id$cluster_id <- rep(seq_len(12), each = 10)
  dat_missing_id$cluster_id[1] <- NA
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat_missing_id,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      B = 10,
      id = "cluster_id"
    )
  )
  dat_one_cluster <- dat
  dat_one_cluster$cluster_id <- 1
  expect_error(
    qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat_one_cluster,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      B = 10,
      id = "cluster_id"
    )
  )
})
