# ------------------------------------------------------------------------------
# Progress helpers
# ------------------------------------------------------------------------------

#' @keywords internal
#' @noRd
qgcompmulti_progress_update_every <- function(B) {
  B <- as.integer(B)
  if (B <= 50L) {
    return(1L)
  }
  if (B <= 100L) {
    return(2L)
  }
  if (B <= 500L) {
    return(5L)
  }
  max(10L, floor(B / 100L))
}

#' @keywords internal
#' @noRd
qgcompmulti_format_duration <- function(seconds) {
  seconds <- as.integer(max(0, round(seconds)))
  hours <- seconds %/% 3600L
  minutes <- (seconds %% 3600L) %/% 60L
  seconds <- seconds %% 60L
  if (hours > 0L) {
    return(sprintf("%02d:%02d:%02d", hours, minutes, seconds))
  }
  sprintf("%02d:%02d", minutes, seconds)
}

#' @keywords internal
#' @noRd
qgcompmulti_progress_bar <- function(attempted, total, width = 32L) {
  attempted <- as.integer(attempted)
  total <- as.integer(total)
  width <- as.integer(width)
  chars <- rep(".", width)
  if (attempted >= total) {
    chars[] <- "="
  } else {
    filled <- floor(width * attempted / total)
    if (filled > 0L) {
      chars[seq_len(filled)] <- "="
    }
    chars[min(width, filled + 1L)] <- ">"
  }
  paste(chars, collapse = "")
}

#' @keywords internal
#' @noRd
qgcompmulti_progress_eta <- function(state) {
  if (state$attempted < 3L) {
    return(NA_real_)
  }
  average_iter_time <- state$cumulative_iteration_time / state$attempted
  average_iter_time * (state$total - state$attempted)
}

#' @keywords internal
#' @noRd
qgcompmulti_progress_status <- function(state) {
  percent <- floor(100 * state$attempted / state$total)
  eta_seconds <- qgcompmulti_progress_eta(state)
  status <- paste0(
    "Bootstrap ",
    sprintf("%3d%%", percent),
    " [",
    qgcompmulti_progress_bar(
      attempted = state$attempted,
      total = state$total,
      width = state$bar_width
    ),
    "] ",
    state$attempted,
    "/",
    state$total,
    "  elapsed ",
    qgcompmulti_format_duration(unname(proc.time()[["elapsed"]]) - state$start_time),
    "  eta ",
    if (is.na(eta_seconds)) "--:--" else qgcompmulti_format_duration(eta_seconds)
  )
  if (state$failures > 0L) {
    status <- paste0(status, "  failed ", state$failures)
  }
  status
}

#' @keywords internal
#' @noRd
qgcompmulti_progress_render <- function(state) {
  if (!isTRUE(state$enabled)) {
    return(state)
  }
  status <- qgcompmulti_progress_status(state)
  padded_width <- max(state$last_width, nchar(status, type = "width"))
  cat("\r", sprintf(paste0("%-", padded_width, "s"), status), sep = "")
  if (interactive()) {
    utils::flush.console()
  }
  state$last_width <- padded_width
  state
}

#' @keywords internal
#' @noRd
qgcompmulti_progress_init <- function(B, enabled = FALSE) {
  now <- unname(proc.time()[["elapsed"]])
  state <- list(
    enabled = isTRUE(enabled),
    total = as.integer(B),
    attempted = 0L,
    failures = 0L,
    start_time = now,
    last_time = now,
    cumulative_iteration_time = 0,
    update_every = qgcompmulti_progress_update_every(B),
    bar_width = 32L,
    last_width = 0L
  )
  if (state$enabled) {
    state <- qgcompmulti_progress_render(state)
  }
  state
}

#' @keywords internal
#' @noRd
qgcompmulti_progress_tick <- function(state, iteration, failed = FALSE) {
  if (!isTRUE(state$enabled)) {
    return(state)
  }
  now <- unname(proc.time()[["elapsed"]])
  iteration_time <- max(0, now - state$last_time)
  state$last_time <- now
  state$attempted <- as.integer(iteration)
  state$cumulative_iteration_time <- state$cumulative_iteration_time + iteration_time
  if (isTRUE(failed)) {
    state$failures <- state$failures + 1L
  }
  if (iteration %% state$update_every == 0L || iteration == state$total) {
    state <- qgcompmulti_progress_render(state)
  }
  state
}

#' @keywords internal
#' @noRd
qgcompmulti_progress_finish <- function(state) {
  if (!is.list(state) || !isTRUE(state$enabled)) {
    return(invisible(state))
  }
  cat("\n", sep = "")
  invisible(state)
}
