.qgcompmulti_test_cache <- new.env(parent = emptyenv())
EXPECTED_TOP_COMPONENTS <- c(
  "call",
  "formula",
  "data_info",
  "mixtures",
  "analysis",
  "fits",
  "prediction",
  "bootstrap",
  "results",
  "labels"
)
EXPECTED_COMPONENT_FIELDS <- list(
  data_info = c(
    "n_input",
    "n_used",
    "outcome",
    "has_clusters",
    "cluster_var",
    "n_clusters",
    "quantized"
  ),
  mixtures = c(
    "mix1",
    "mix2",
    "q",
    "centering"
  ),
  analysis = c(
    "interaction",
    "family",
    "family_name",
    "link",
    "estimand_scale",
    "estimand_scale_defaulted",
    "msm_fitting_scale",
    "default_interval_method",
    "B",
    "id",
    "MCsize",
    "seed"
  ),
  fits = c(
    "outcome_fit",
    "msm_fit"
  ),
  prediction = c(
    "intervention_grid",
    "msm_grid",
    "counterfactual_surface",
    "msm_surface",
    "surface_comparison",
    "counterfactual_surface_target",
    "msm_surface_target",
    "surface_comparison_target"
  ),
  bootstrap = c(
    "coef_draws",
    "B_requested",
    "B_success",
    "B_failed",
    "failure_log"
  ),
  results = c(
    "coefficients",
    "std_error",
    "vcov",
    "coef_table"
  ),
  labels = c(
    "mixture_labels",
    "coefficient_names",
    "coefficient_labels"
  )
)
EXPECTED_MI_TOP_COMPONENTS <- c(
  "call",
  "formula",
  "data_info",
  "mixtures",
  "mi_info",
  "analysis",
  "fits",
  "results",
  "labels"
)
EXPECTED_MI_COMPONENT_FIELDS <- list(
  data_info = c(
    "n_input",
    "n_used",
    "outcome",
    "has_clusters",
    "cluster_var",
    "n_clusters",
    "quantized"
  ),
  mixtures = c(
    "mix1",
    "mix2",
    "q",
    "centering"
  ),
  mi_info = c(
    "m",
    "input_type",
    "keep_fits",
    "seed",
    "fit_seeds",
    "n_input_per_imputation",
    "n_used_per_imputation"
  ),
  analysis = c(
    "interaction",
    "family",
    "family_name",
    "link",
    "estimand_scale",
    "estimand_scale_defaulted",
    "msm_fitting_scale",
    "default_interval_method",
    "B",
    "id",
    "MCsize"
  ),
  fits = c("imputation_fits"),
  results = c(
    "coefficients",
    "std_error",
    "vcov",
    "coef_table",
    "within_var",
    "between_var",
    "total_var",
    "df",
    "riv",
    "fmi"
  ),
  labels = c(
    "mixture_labels",
    "coefficient_names",
    "coefficient_labels"
  )
)
EXPECTED_SUMMARY_COMPONENTS <- c(
  "call",
  "formula",
  "fit_overview",
  "mixtures",
  "msm_table",
  "outcome_model_info",
  "labels"
)
EXPECTED_COEF_NAMES_WITH_INTERACTION <- c(
  "(Intercept)",
  "psi1",
  "psi2",
  "psi1:psi2"
)
EXPECTED_COEF_NAMES_NO_INTERACTION <- c(
  "(Intercept)",
  "psi1",
  "psi2"
)
EXPECTED_GRID_COLUMNS <- c("grid_id", "psi1", "psi2")
EXPECTED_COUNTERFACTUAL_SURFACE_COLUMNS <- c(
  "grid_id",
  "intervention_psi1",
  "intervention_psi2",
  "msm_psi1",
  "msm_psi2",
  "exact_mean"
)
EXPECTED_MSM_SURFACE_COLUMNS <- c(
  "grid_id",
  "intervention_psi1",
  "intervention_psi2",
  "msm_psi1",
  "msm_psi2",
  "msm_mean"
)
EXPECTED_SURFACE_COMPARISON_COLUMNS <- c(
  "grid_id",
  "intervention_psi1",
  "intervention_psi2",
  "msm_psi1",
  "msm_psi2",
  "exact_mean",
  "msm_mean",
  "residual"
)
EXPECTED_COUNTERFACTUAL_TARGET_COLUMNS <- c(
  "grid_id",
  "intervention_psi1",
  "intervention_psi2",
  "msm_psi1",
  "msm_psi2",
  "exact_target"
)
EXPECTED_MSM_TARGET_COLUMNS <- c(
  "grid_id",
  "intervention_psi1",
  "intervention_psi2",
  "msm_psi1",
  "msm_psi2",
  "msm_target"
)
EXPECTED_SURFACE_COMPARISON_TARGET_COLUMNS <- c(
  "grid_id",
  "intervention_psi1",
  "intervention_psi2",
  "msm_psi1",
  "msm_psi2",
  "exact_target",
  "msm_target",
  "residual_target"
)
EXPECTED_INTERNAL_PREDICTION_FIELDS <- c(
  "prediction_type",
  "grid_type",
  "grid_scale",
  "estimand_scale",
  "estimates",
  "intervals",
  "interval_type",
  "uncertainty_source",
  "data_supplied",
  "contrast"
)

make_test_data <- function(seed = 123, n = 120) {
  key <- sprintf("data_%s_%s", seed, n)
  if (!exists(key, envir = .qgcompmulti_test_cache, inherits = FALSE)) {
    dat <- sim_mixture_data(
      n = n,
      pA = 3,
      pB = 3,
      rho_within_A = 0.3,
      rho_within_B = 0.3,
      rho_between = 0.2,
      psi1 = 0.5,
      psi2 = 0.3,
      psi12 = 0.2,
      return_quantized = FALSE,
      seed = seed
    )
    assign(key, dat, envir = .qgcompmulti_test_cache)
  }
  get(key, envir = .qgcompmulti_test_cache, inherits = FALSE)
}

make_clustered_test_data <- function(seed = 123, n = 120, cluster_size = 10) {
  stopifnot(n %% cluster_size == 0L)
  key <- sprintf("clustered_data_%s_%s_%s", seed, n, cluster_size)
  if (!exists(key, envir = .qgcompmulti_test_cache, inherits = FALSE)) {
    dat <- make_test_data(seed = seed, n = n)
    dat$cluster_id <- rep(seq_len(n / cluster_size), each = cluster_size)
    assign(key, dat, envir = .qgcompmulti_test_cache)
  }
  get(key, envir = .qgcompmulti_test_cache, inherits = FALSE)
}

make_uneven_clustered_test_data <- function(
    seed = 123,
    cluster_sizes = c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 7, 8)
) {
  n <- sum(cluster_sizes)
  key <- sprintf(
    "uneven_clustered_data_%s_%s",
    seed,
    paste(cluster_sizes, collapse = "_")
  )
  if (!exists(key, envir = .qgcompmulti_test_cache, inherits = FALSE)) {
    dat <- make_test_data(seed = seed, n = n)
    dat$cluster_id <- rep(seq_along(cluster_sizes), times = cluster_sizes)
    assign(key, dat, envir = .qgcompmulti_test_cache)
  }
  get(key, envir = .qgcompmulti_test_cache, inherits = FALSE)
}

make_test_fit_args <- function(data) {
  list(
    f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
    data = data,
    mix1 = c("X1", "X2", "X3"),
    mix2 = c("W1", "W2", "W3"),
    family = gaussian(),
    B = 10,
    MCsize = nrow(data)
  )
}

fit_test_model <- function(interaction = TRUE,
                           q = 4,
                           centering = "none",
                           mcsize = NULL,
                           clustered = FALSE,
                           B = 10,
                           seed = 123) {
  dat <- if (clustered) {
    make_clustered_test_data(seed = seed)
  } else {
    make_test_data(seed = seed)
  }
  if (is.null(mcsize)) {
    mcsize <- nrow(dat)
  }
  q_key <- if (is.null(q)) "NULL" else as.character(q)
  key <- paste(
    "fit",
    interaction,
    q_key,
    centering,
    mcsize,
    clustered,
    B,
    seed,
    sep = "_"
  )
  if (!exists(key, envir = .qgcompmulti_test_cache, inherits = FALSE)) {
    fit <- qgcomp.glm.multi(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = interaction,
      q = q,
      centering = centering,
      B = B,
      id = if (clustered) "cluster_id" else NULL,
      MCsize = mcsize
    )
    assign(key, fit, envir = .qgcompmulti_test_cache)
  }
  get(key, envir = .qgcompmulti_test_cache, inherits = FALSE)
}
fit_test_model_list <- function(interaction = TRUE,
                                q = 4,
                                centering = "none",
                                clustered = FALSE,
                                B = 10,
                                data_seeds = c(101, 202, 303),
                                fit_seeds = data_seeds + 1000L) {
  stopifnot(length(data_seeds) == length(fit_seeds))

  key <- paste(
    "fit_list",
    interaction,
    if (is.null(q)) "NULL" else as.character(q),
    centering,
    clustered,
    B,
    paste(data_seeds, collapse = "_"),
    paste(fit_seeds, collapse = "_"),
    sep = "_"
  )

  if (!exists(key, envir = .qgcompmulti_test_cache, inherits = FALSE)) {
    fits <- lapply(
      seq_along(data_seeds),
      function(i) {
        dat <- if (clustered) {
          make_clustered_test_data(seed = data_seeds[[i]])
        } else {
          make_test_data(seed = data_seeds[[i]])
        }
        qgcomp.glm.multi(
          f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
          data = dat,
          mix1 = c("X1", "X2", "X3"),
          mix2 = c("W1", "W2", "W3"),
          interaction = interaction,
          q = q,
          centering = centering,
          B = B,
          id = if (clustered) "cluster_id" else NULL,
          MCsize = nrow(dat),
          seed = fit_seeds[[i]]
        )
      }
    )
    assign(key, fits, envir = .qgcompmulti_test_cache)
  }

  get(key, envir = .qgcompmulti_test_cache, inherits = FALSE)
}

make_mi_incomplete_data <- function(seed = 123, clustered = FALSE) {
  key <- paste("mi_incomplete", seed, clustered, sep = "_")

  if (!exists(key, envir = .qgcompmulti_test_cache, inherits = FALSE)) {
    dat <- if (clustered) {
      make_clustered_test_data(seed = seed)
    } else {
      make_test_data(seed = seed)
    }

    dat$X1[seq(5, nrow(dat), by = 11)] <- NA_real_
    dat$W2[seq(3, nrow(dat), by = 13)] <- NA_real_
    dat$C[seq(7, nrow(dat), by = 17)] <- NA_real_

    assign(key, dat, envir = .qgcompmulti_test_cache)
  }

  get(key, envir = .qgcompmulti_test_cache, inherits = FALSE)
}

make_completed_data_list <- function(m = 3, seed = 123, clustered = FALSE) {
  stopifnot(m >= 1L)

  key <- paste("mi_completed_list", m, seed, clustered, sep = "_")

  if (!exists(key, envir = .qgcompmulti_test_cache, inherits = FALSE)) {
    dat <- make_mi_incomplete_data(seed = seed, clustered = clustered)

    x1_missing <- is.na(dat$X1)
    w2_missing <- is.na(dat$W2)
    c_missing <- is.na(dat$C)

    x1_fill <- mean(dat$X1, na.rm = TRUE)
    w2_fill <- mean(dat$W2, na.rm = TRUE)
    c_fill <- mean(dat$C, na.rm = TRUE)

    completed <- lapply(
      seq_len(m),
      function(i) {
        dat_i <- dat
        dat_i$X1[x1_missing] <- x1_fill + (i - ((m + 1) / 2)) * 0.05
        dat_i$W2[w2_missing] <- w2_fill - (i - ((m + 1) / 2)) * 0.05
        dat_i$C[c_missing] <- c_fill + (i - ((m + 1) / 2)) * 0.02
        dat_i
      }
    )

    assign(key, completed, envir = .qgcompmulti_test_cache)
  }

  get(key, envir = .qgcompmulti_test_cache, inherits = FALSE)
}

make_test_mids <- function(m = 3, seed = 123, clustered = FALSE) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("`mice` must be installed to build test mids objects.", call. = FALSE)
  }

  key <- paste("mi_mids", m, seed, clustered, sep = "_")

  if (!exists(key, envir = .qgcompmulti_test_cache, inherits = FALSE)) {
    dat <- make_mi_incomplete_data(seed = seed, clustered = clustered)
    mids <- mice::mice(
      dat,
      m = m,
      maxit = 1,
      printFlag = FALSE,
      seed = seed
    )
    assign(key, mids, envir = .qgcompmulti_test_cache)
  }

  get(key, envir = .qgcompmulti_test_cache, inherits = FALSE)
}

fit_engine_result <- function(interaction = TRUE,
                              q = 4,
                              centering = "none",
                              mcsize = NULL,
                              clustered = FALSE,
                              seed = 123) {
  dat <- if (clustered) {
    make_clustered_test_data(seed = seed)
  } else {
    make_test_data(seed = seed)
  }
  if (is.null(mcsize)) {
    mcsize <- nrow(dat)
  }
  q_key <- if (is.null(q)) "NULL" else as.character(q)
  key <- paste(
    "engine",
    interaction,
    q_key,
    centering,
    mcsize,
    clustered,
    seed,
    sep = "_"
  )
  if (!exists(key, envir = .qgcompmulti_test_cache, inherits = FALSE)) {
    result <- qgcompmulti_msm_fit(
      f = Y ~ X1 + X2 + X3 + W1 + W2 + W3 + C,
      data = dat,
      mix1 = c("X1", "X2", "X3"),
      mix2 = c("W1", "W2", "W3"),
      interaction = interaction,
      q = q,
      centering = centering,
      id = if (clustered) "cluster_id" else NULL,
      MCsize = mcsize
    )
    assign(key, result, envir = .qgcompmulti_test_cache)
  }
  get(key, envir = .qgcompmulti_test_cache, inherits = FALSE)
}
