is_scalar_whole_number <- function(x) {
  is.numeric(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    is.finite(x) &&
    x == as.integer(x)
}

validate_mix_names <- function(mix, arg) {
  if (!is.character(mix) || length(mix) == 0L || anyNA(mix)) {
    stop(sprintf("`%s` must be a non-empty character vector.", arg), call. = FALSE)
  }

  if (anyDuplicated(mix)) {
    stop(sprintf("`%s` must not contain duplicated variable names.", arg), call. = FALSE)
  }
}

validate_q_argument <- function(q) {
  if (is.null(q)) {
    return(invisible(NULL))
  }

  if (!is_scalar_whole_number(q) || q < 2L) {
    stop("`q` must be `NULL` or a single integer greater than or equal to 2.", call. = FALSE)
  }
}

validate_qgcomp_multi_inputs <- function(
    f,
    data,
    mix1,
    mix2,
    interaction,
    family,
    q,
    centering,
    id,
    MCsize,
    B = NULL,
    seed = NULL) {

  if (!inherits(f, "formula")) {
    stop("`f` must be a formula.", call. = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  if (nrow(data) == 0L) {
    stop("`data` must contain at least one row.", call. = FALSE)
  }

  validate_mix_names(mix1, "mix1")
  validate_mix_names(mix2, "mix2")

  if (length(intersect(mix1, mix2)) > 0L) {
    stop("`mix1` and `mix2` must not contain overlapping variables.", call. = FALSE)
  }

  required_vars <- unique(c(all.vars(f), mix1, mix2, id))
  missing_vars <- setdiff(required_vars, names(data))

  if (length(missing_vars) > 0L) {
    stop(
      sprintf(
        "The following variables are missing from `data`: %s",
        paste(missing_vars, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  formula_vars <- all.vars(f)
  missing_mix_terms <- setdiff(c(mix1, mix2), formula_vars)

  if (length(missing_mix_terms) > 0L) {
    stop(
      sprintf(
        "All mixture variables must appear in `f`. Missing from formula: %s",
        paste(missing_mix_terms, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.logical(interaction) || length(interaction) != 1L || is.na(interaction)) {
    stop("`interaction` must be either `TRUE` or `FALSE`.", call. = FALSE)
  }

  if (!inherits(family, "family")) {
    stop("`family` must be a valid GLM family object, such as gaussian(), binomial(), or poisson().", call. = FALSE)
  }

  validate_q_argument(q)

  if (!is.character(centering) || length(centering) != 1L || is.na(centering)) {
    stop("`centering` must be a single character string.", call. = FALSE)
  }

  if (!centering %in% c("none", "median")) {
    stop("`centering` must be one of \"none\" or \"median\".", call. = FALSE)
  }

  if (!is.null(q) && centering != "none") {
    warning("`centering` is ignored unless `q = NULL`.", call. = FALSE)
  }

  if (!is.null(B)) {
    if (!is_scalar_whole_number(B) || B < 2L) {
      stop("`B` must be a single integer greater than or equal to 2.", call. = FALSE)
    }
  }

  if (!is.null(seed)) {
    if (!is_scalar_whole_number(seed)) {
      stop("`seed` must be `NULL` or a single integer.", call. = FALSE)
    }
  }

  if (!is.null(id)) {
    if (!is.character(id) || length(id) != 1L || is.na(id)) {
      stop("`id` must be `NULL` or a single character string.", call. = FALSE)
    }

    if (anyNA(data[[id]])) {
      stop("`id` must not contain missing values.", call. = FALSE)
    }

    if (length(unique(data[[id]])) < 2L) {
      stop("Clustered resampling requires at least two unique cluster IDs.", call. = FALSE)
    }
  }

  if (!is_scalar_whole_number(MCsize) || MCsize < 1L) {
    stop("`MCsize` must be a single integer greater than or equal to 1.", call. = FALSE)
  }

  if (MCsize > nrow(data)) {
    stop("`MCsize` must be less than or equal to `nrow(data)`.", call. = FALSE)
  }

  invisible(NULL)
}

validate_quantize_mixtures_inputs <- function(data, mix1, mix2, q) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  if (nrow(data) == 0L) {
    stop("`data` must contain at least one row.", call. = FALSE)
  }

  validate_mix_names(mix1, "mix1")
  validate_mix_names(mix2, "mix2")

  if (length(intersect(mix1, mix2)) > 0L) {
    stop("`mix1` and `mix2` must not contain overlapping variables.", call. = FALSE)
  }

  required_vars <- unique(c(mix1, mix2))
  missing_vars <- setdiff(required_vars, names(data))

  if (length(missing_vars) > 0L) {
    stop(
      sprintf(
        "The following mixture variables are missing from `data`: %s",
        paste(missing_vars, collapse = ", ")
      ), call. = FALSE
    )
  }

  validate_q_argument(q)

  invisible(NULL)
}
