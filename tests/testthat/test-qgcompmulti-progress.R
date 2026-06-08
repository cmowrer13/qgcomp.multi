test_that("qgcomp.glm.multi progress output is silent by default", {
  dat <- make_test_data()
  out <- capture.output({
    fit <- qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = TRUE,
      q = 4,
      B = 6,
      MCsize = 40,
      seed = 123,
      progress = FALSE
    )
  })
  expect_length(out, 0L)
})

test_that("qgcomp.glm.multi progress output appears when enabled", {
  dat <- make_test_data()
  out <- paste(
    capture.output({
      fit <- qgcomp.glm.multi(
        f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
        data = dat,
        mix1 = c("X1", "X2", "X3"),
        mix2 = c("W1", "W2", "W3"),
        interaction = TRUE,
        q = 4,
        B = 6,
        MCsize = 40,
        seed = 123,
        progress = TRUE
      )
    }),
    collapse = ""
  )
  expect_match(out, "Bootstrap")
  expect_match(out, "eta")
  expect_false(grepl("failed", out, fixed = TRUE))
})

test_that("progress helper only shows failed counter after a failure", {
  out_no_fail <- paste(
    capture.output({
      state <- qgcompmulti_progress_init(B = 5, enabled = TRUE)
      qgcompmulti_progress_finish(state)
    }),
    collapse = ""
  )
  expect_false(grepl("failed", out_no_fail, fixed = TRUE))

  out_fail <- paste(
    capture.output({
      state <- qgcompmulti_progress_init(B = 5, enabled = TRUE)
      state$attempted <- 3L
      state$failures <- 1L
      state$cumulative_iteration_time <- 0.3
      state <- qgcompmulti_progress_render(state)
      qgcompmulti_progress_finish(state)
    }),
    collapse = ""
  )
  expect_true(grepl("failed 1", out_fail, fixed = TRUE))
})
