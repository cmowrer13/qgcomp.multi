# ------------------------------------------------------------------------------
# Multiple-imputation fit helpers
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
qgcompmulti_mi_input_type <- function(data) {
  if (inherits(data, "mids")) {
    return("mids")
  }
  if (is.list(data)) {
    return("completed_list")
  }
  stop(
    "`data` must be either a `mids` object or a non-empty list of completed data frames.",
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_completed_data <- function(data, input_type = qgcompmulti_mi_input_type(data)) {
  completed <- if (identical(input_type, "mids")) {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("The optional `mice` package must be installed to use `mids` inputs.", call. = FALSE)
    }
    mice::complete(data, action = "all")
  } else {
    data
  }

  if (!is.list(completed) || length(completed) == 0L) {
    stop("Completed imputations must be supplied as a non-empty list of data frames.", call. = FALSE)
  }

  invalid_idx <- which(!vapply(completed, is.data.frame, logical(1)))
  if (length(invalid_idx) > 0L) {
    stop(
      sprintf(
        "All completed imputations must be data frames. Problem at position(s): %s",
        paste(invalid_idx, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  completed
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_validate_completed_data <- function(completed,
                                                   f,
                                                   mix1,
                                                   mix2,
                                                   interaction,
                                                   family,
                                                   estimand_scale = NULL,
                                                   q,
                                                   centering,
                                                   id,
                                                   MCsize,
                                                   B,
                                                   progress = FALSE,
                                                   parallel = FALSE,
                                                   workers = NULL) {
  if (!is.list(completed) || length(completed) == 0L) {
    stop("`completed` must be a non-empty list of data frames.", call. = FALSE)
  }

  reference_names <- names(completed[[1]])
  reference_n <- nrow(completed[[1]])
  required_vars <- unique(c(all.vars(f), mix1, mix2, id))
  reference_ids <- if (is.null(id)) NULL else completed[[1]][[id]]

  for (i in seq_along(completed)) {
    dat <- completed[[i]]

    if (!identical(names(dat), reference_names)) {
      stop("All completed imputations must have identical column names and ordering.", call. = FALSE)
    }

    if (nrow(dat) != reference_n) {
      stop("All completed imputations must have the same number of rows.", call. = FALSE)
    }

    missing_vars <- setdiff(required_vars, names(dat))
    if (length(missing_vars) > 0L) {
      stop(
        sprintf(
          "Completed imputation %s is missing required variable(s): %s",
          i,
          paste(missing_vars, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    if (length(required_vars) > 0L && anyNA(dat[, required_vars, drop = FALSE])) {
      stop(
        sprintf(
          "Completed imputation %s still contains missing values in required analysis variables.",
          i
        ),
        call. = FALSE
      )
    }

    if (!is.null(id)) {
      if (anyNA(dat[[id]])) {
        stop(
          sprintf("Completed imputation %s contains missing cluster IDs.", i),
          call. = FALSE
        )
      }

      if (length(unique(dat[[id]])) < 2L) {
        stop("Clustered resampling requires at least two unique cluster IDs.", call. = FALSE)
      }

      if (!identical(dat[[id]], reference_ids)) {
        stop(
          "All completed imputations must share identical cluster IDs and ordering when `id` is supplied.",
          call. = FALSE
        )
      }
    }
  }

  validate_qgcomp_multi_inputs(
    f = f,
    data = completed[[1]],
    mix1 = mix1,
    mix2 = mix2,
    interaction = interaction,
    family = family,
    estimand_scale = estimand_scale,
    q = q,
    centering = centering,
    id = id,
    MCsize = MCsize,
    B = B,
    seed = NULL,
    progress = progress,
    parallel = parallel,
    workers = workers
  )

  invisible(completed)
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_fit_seeds <- function(m, seed = NULL) {
  if (!is_scalar_whole_number(m) || m < 1L) {
    stop("`m` must be a positive integer.", call. = FALSE)
  }

  if (is.null(seed)) {
    return(rep(list(NULL), m))
  }

  derived <- qgcompmulti_with_seed(
    seed,
    sample.int(.Machine$integer.max, size = m, replace = FALSE)
  )

  as.list(as.integer(derived))
}

#' @keywords internal
#' @noRd
qgcompmulti_mi_fit_models <- function(completed,
                                      f,
                                      mix1,
                                      mix2,
                                      interaction,
                                      family,
                                      estimand_scale = NULL,
                                      q,
                                      centering,
                                      B,
                                      id,
                                      MCsize,
                                      fit_seeds,
                                      progress = FALSE,
                                      parallel = FALSE,
                                      workers = NULL) {
  if (!is.list(completed) || length(completed) == 0L) {
    stop("`completed` must be a non-empty list of data frames.", call. = FALSE)
  }

  if (!is.list(fit_seeds) || length(fit_seeds) != length(completed)) {
    stop("`fit_seeds` must be a list with one entry per completed imputation.", call. = FALSE)
  }

  lapply(
    seq_along(completed),
    function(i) {
      qgcomp.glm.multi(
        f = f,
        data = completed[[i]],
        mix1 = mix1,
        mix2 = mix2,
        interaction = interaction,
        family = family,
        estimand_scale = estimand_scale,
        q = q,
        centering = centering,
        B = B,
        id = id,
        MCsize = MCsize,
        seed = fit_seeds[[i]],
        progress = progress,
        parallel = parallel,
        workers = workers
      )
    }
  )
}
