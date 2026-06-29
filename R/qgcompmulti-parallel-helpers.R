# ------------------------------------------------------------------------------
# Parallel bootstrap helpers
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
qgcompmulti_bootstrap_seed_plan <- function(B, seed = NULL) {
  if (!is_scalar_whole_number(B) || B < 1L) {
    stop("`B` must be a positive integer.", call. = FALSE)
  }

  if (is.null(seed)) {
    return(list(
      full_fit_seed = NULL,
      worker_seeds = rep(list(NULL), B)
    ))
  }

  derived <- qgcompmulti_with_seed(
    seed,
    sample.int(.Machine$integer.max, size = B + 1L, replace = FALSE)
  )

  list(
    full_fit_seed = as.integer(derived[[1L]]),
    worker_seeds = as.list(as.integer(derived[-1L]))
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_parallel_progress_flag <- function(parallel = FALSE, progress = FALSE) {
  if (!parallel) {
    return(progress)
  }

  if (isTRUE(progress)) {
    warning(
      "`progress = TRUE` is not supported when `parallel = TRUE`; progress display has been disabled.",
      call. = FALSE
    )
  }

  FALSE
}

#' @keywords internal
#' @noRd
qgcompmulti_bootstrap_replicate_fit <- function(replicate,
                                                data,
                                                f,
                                                mix1,
                                                mix2,
                                                interaction,
                                                family,
                                                estimand_scale,
                                                q,
                                                centering,
                                                id,
                                                MCsize,
                                                seed = NULL) {
  boot_fit <- tryCatch(
    qgcompmulti_with_seed(
      seed,
      {
        data_b <- qgcompmulti_resample_data(
          data = data,
          id = id
        )

        qgcompmulti_msm_fit(
          f = f,
          data = data_b,
          mix1 = mix1,
          mix2 = mix2,
          interaction = interaction,
          family = family,
          estimand_scale = estimand_scale,
          q = q,
          centering = centering,
          id = id,
          MCsize = MCsize,
          seed = NULL
        )
      }
    ),
    error = function(e) e
  )

  if (inherits(boot_fit, "error")) {
    return(list(
      replicate = as.integer(replicate),
      coefficients = NULL,
      message = conditionMessage(boot_fit)
    ))
  }

  list(
    replicate = as.integer(replicate),
    coefficients = boot_fit$coefficients,
    message = NULL
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_parallel_worker_count <- function(n_tasks, workers = NULL) {
  if (!is_scalar_whole_number(n_tasks) || n_tasks < 1L) {
    stop("`n_tasks` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(workers)) {
    return(as.integer(workers))
  }

  available <- suppressWarnings(future::availableCores())
  available <- as.integer(available[[1L]])

  if (is.na(available) || available < 1L) {
    available <- 1L
  }

  min(as.integer(n_tasks), available)
}

#' @keywords internal
#' @noRd
qgcompmulti_run_bootstrap_replicates <- function(B,
                                                 worker_fun,
                                                 parallel = FALSE,
                                                 workers = NULL,
                                                 progress = FALSE) {
  if (!is.function(worker_fun)) {
    stop("`worker_fun` must be a function.", call. = FALSE)
  }

  indices <- seq_len(B)

  if (!parallel) {
    progress_state <- qgcompmulti_progress_init(B = B, enabled = progress)
    on.exit(qgcompmulti_progress_finish(progress_state), add = TRUE)

    results <- vector("list", B)

    for (i in indices) {
      results[[i]] <- worker_fun(i)
      progress_state <- qgcompmulti_progress_tick(
        progress_state,
        iteration = i,
        failed = !is.null(results[[i]]$message)
      )
    }

    return(results)
  }

  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop(
      "Parallel execution requires the `future` and `future.apply` packages.",
      call. = FALSE
    )
  }

  current_plan <- future::plan()
  using_active_plan <- !inherits(current_plan, "sequential")

  if (!is.null(workers) && using_active_plan) {
    stop(
      paste(
        "A non-sequential future plan is already active.",
        "Leave `workers = NULL` to use that plan or reset the plan before supplying `workers`."
      ),
      call. = FALSE
    )
  }

  if (using_active_plan) {
    return(
      future.apply::future_lapply(
        indices,
        worker_fun,
        future.packages = "qgcomp.multi",
        future.seed = FALSE,
        future.label = "qgcompmulti-bootstrap-%d"
      )
    )
  }

  local_workers <- qgcompmulti_parallel_worker_count(
    n_tasks = B,
    workers = workers
  )

  on.exit(future::plan(current_plan), add = TRUE)
  future::plan(future::multisession, workers = local_workers)

  future.apply::future_lapply(
    indices,
    worker_fun,
    future.packages = "qgcomp.multi",
    future.seed = FALSE,
    future.label = "qgcompmulti-bootstrap-%d"
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_collect_bootstrap_results <- function(results, coef_names) {
  if (!is.list(results) || length(results) == 0L) {
    stop("`results` must be a non-empty list.", call. = FALSE)
  }

  if (!is.character(coef_names) || length(coef_names) == 0L || anyNA(coef_names)) {
    stop("`coef_names` must be a non-empty character vector.", call. = FALSE)
  }

  successful <- Filter(
    function(x) is.list(x) && !is.null(x$coefficients),
    results
  )

  if (length(successful) == 0L) {
    stop("All bootstrap replications failed.", call. = FALSE)
  }

  coef_draws <- do.call(
    rbind,
    lapply(
      successful,
      function(x) {
        coefs <- x$coefficients[coef_names]
        matrix(
          as.numeric(coefs),
          nrow = 1L,
          dimnames = list(NULL, coef_names)
        )
      }
    )
  )

  failures <- Filter(
    function(x) is.list(x) && !is.null(x$message),
    results
  )

  failure_log <- if (length(failures) == 0L) {
    NULL
  } else {
    do.call(
      rbind,
      lapply(
        failures,
        function(x) {
          data.frame(
            replicate = as.integer(x$replicate),
            message = x$message,
            row.names = NULL,
            check.names = FALSE
          )
        }
      )
    )
  }

  list(
    coef_draws = coef_draws,
    failure_log = failure_log
  )
}
