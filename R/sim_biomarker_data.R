#' Simulate a biomarker dataset for two-mixture inflammation analyses
#'
#' Generates a narrative environmental-health teaching dataset with two
#' correlated exposure mixtures, a continuous inflammation outcome, and a small
#' amount of intentional model imperfection. The simulated biomarkers are
#' positive, right-skewed concentrations that naturally motivate log
#' transformation before model fitting.
#'
#' Mixture 1 represents urinary phthalate metabolites:
#' \itemize{
#'   \item `mep_ng_ml`
#'   \item `mibp_ng_ml`
#'   \item `mehpp_ng_ml`
#' }
#'
#' Mixture 2 represents urinary phenol biomarkers:
#' \itemize{
#'   \item `bpa_ng_ml`
#'   \item `bp3_ng_ml`
#'   \item `ppb_ng_ml`
#' }
#'
#' The outcome `inflammation_score` is generated from weighted mixture summaries
#' on the log-concentration scale, a moderate positive interaction between the
#' two mixtures, common epidemiologic covariates, and a small amount of
#' deliberate nonlinearity. In the current implementation, that nonlinearity
#' comes from a centered cubic term in the phthalate mixture score and a mild
#' threshold term in the phenol mixture score. This design makes the fitted MSM
#' informative without making the exact fit-time surface unrealistically perfect.
#'
#' @param n Number of observations to generate.
#' @param include_log Logical; if `TRUE`, also include log-transformed biomarker
#'   variables named `ln_mep`, `ln_mibp`, `ln_mehpp`, `ln_bpa`, `ln_bp3`, and
#'   `ln_ppb`.
#' @param seed Optional integer random seed. When supplied, the generated
#'   dataset is reproducible and the caller's RNG state is restored on exit.
#'
#' @return A data frame containing `inflammation_score`, six raw biomarker
#'   concentrations, and the covariates `age_years`, `bmi_kg_m2`, and
#'   `smoker`. If `include_log = TRUE`, log-transformed biomarker variables are
#'   included as well.
#'
#' @details
#' This helper is intended for examples, documentation, and teaching. The
#' exposure concentrations are simulated as correlated log-normal variables with
#' moderate within-mixture correlation, mild between-mixture correlation, and
#' deliberately different raw scales across chemicals.
#'
#' Those scale differences are useful for demonstrating:
#' \itemize{
#'   \item standard quantized fitting,
#'   \item original-scale fitting with `q = NULL`,
#'   \item why log transformation can be scientifically natural for biomarker
#'     concentrations, and
#'   \item why median centering can make the `q = NULL` intercept more
#'     interpretable.
#' }
#'
#' @examples
#' dat <- sim_biomarker_data(n = 200, seed = 123)
#' head(dat)
#'
#' dat_log <- sim_biomarker_data(n = 200, include_log = TRUE, seed = 123)
#' head(dat_log[c("ln_mep", "ln_bpa", "inflammation_score")])
#'
#' @export
sim_biomarker_data <- function(n = 800, include_log = FALSE, seed = NULL) {
  qgcompmulti_validate_sim_biomarker_input(
    n = n,
    include_log = include_log,
    seed = seed
  )

  qgcompmulti_with_seed(
    seed,
    {
      log_names <- c("ln_mep", "ln_mibp", "ln_mehpp", "ln_bpa", "ln_bp3", "ln_ppb")
      raw_names <- c("mep_ng_ml", "mibp_ng_ml", "mehpp_ng_ml", "bpa_ng_ml", "bp3_ng_ml", "ppb_ng_ml")

      log_means <- c(4.1, 2.5, 3.2, 0.8, 2.8, 1.2)
      log_sds <- c(0.75, 0.65, 0.70, 0.85, 0.80, 0.90)

      corr <- qgcompmulti_biomarker_correlation()
      sigma <- diag(log_sds) %*% corr %*% diag(log_sds)

      log_mat <- MASS::mvrnorm(
        n = n,
        mu = log_means,
        Sigma = sigma
      )
      colnames(log_mat) <- log_names

      raw_mat <- exp(log_mat)
      colnames(raw_mat) <- raw_names

      age_years <- round(pmin(pmax(stats::rnorm(n, mean = 46, sd = 13), 20), 75))
      bmi_kg_m2 <- pmin(pmax(stats::rnorm(n, mean = 28, sd = 5.5), 18), 45)
      smoker <- stats::rbinom(n, size = 1, prob = 0.24)

      z_mat <- scale(log_mat)
      phthalate_score <- 0.55 * z_mat[, "ln_mep"] +
        0.20 * z_mat[, "ln_mibp"] +
        0.45 * z_mat[, "ln_mehpp"]
      phenol_score <- 0.50 * z_mat[, "ln_bpa"] +
        0.15 * z_mat[, "ln_bp3"] +
        0.40 * z_mat[, "ln_ppb"]

      phthalate_curvature <- (phthalate_score^3 - mean(phthalate_score^3))
      phenol_threshold <- pmax(phenol_score - 1.00, 0) -
        mean(pmax(phenol_score - 1.00, 0))

      inflammation_score <- as.numeric(
        2.6 +
          0.38 * phthalate_score +
          0.42 * phenol_score +
          0.16 * phthalate_score * phenol_score +
          0.07 * phthalate_curvature +
          0.06 * phenol_threshold +
          0.015 * (age_years - 46) +
          0.085 * (bmi_kg_m2 - 28) +
          0.30 * smoker +
          stats::rnorm(n, mean = 0, sd = 0.80)
      )

      dat <- data.frame(
        inflammation_score = inflammation_score,
        raw_mat,
        age_years = age_years,
        bmi_kg_m2 = bmi_kg_m2,
        smoker = smoker,
        check.names = FALSE
      )

      if (isTRUE(include_log)) {
        dat[log_names] <- log_mat
      }

      dat
    }
  )
}

qgcompmulti_biomarker_correlation <- function() {
  corr <- matrix(0.15, nrow = 6, ncol = 6)
  diag(corr) <- 1

  corr[1:3, 1:3] <- matrix(
    c(
      1.00, 0.55, 0.48,
      0.55, 1.00, 0.42,
      0.48, 0.42, 1.00
    ),
    nrow = 3,
    byrow = TRUE
  )

  corr[4:6, 4:6] <- matrix(
    c(
      1.00, 0.32, 0.38,
      0.32, 1.00, 0.28,
      0.38, 0.28, 1.00
    ),
    nrow = 3,
    byrow = TRUE
  )

  corr[1:3, 4:6] <- matrix(
    c(
      0.18, 0.12, 0.16,
      0.14, 0.10, 0.12,
      0.20, 0.16, 0.18
    ),
    nrow = 3,
    byrow = TRUE
  )
  corr[4:6, 1:3] <- t(corr[1:3, 4:6])

  corr
}

qgcompmulti_validate_sim_biomarker_input <- function(n, include_log, seed) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 10 || n %% 1 != 0) {
    stop("`n` must be a single integer of at least 10.", call. = FALSE)
  }

  if (!is.logical(include_log) || length(include_log) != 1L || is.na(include_log)) {
    stop("`include_log` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed %% 1 != 0) {
      stop("`seed` must be `NULL` or a single integer.", call. = FALSE)
    }
  }

  invisible(NULL)
}
