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
    "surface_comparison"
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
