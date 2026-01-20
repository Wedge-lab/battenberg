#' Initialize and Configure Logging
#'
#' Sets up a file-based logger using the `logger` package. It creates the
#' destination directory if it does not already exist and adjusts the
#' logging threshold based on the desired verbosity.
#'
#' @param log_path Character string. The full path to the log file.
#' @param verbose Logical. If `TRUE`, the log level is set to `DEBUG`.
#'   If `FALSE`, it defaults to `INFO`.
#'
#' @export
log_setup <- function(log_path, verbose = FALSE) {
  if (file.info(log_path)$isdir %||% dir.exists(log_path)) {
    log_path <- file.path(log_path, "session.log")
  }

  # Create directory if it doesn't exist
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

  # Set where the log goes
  logger::log_appender(logger::appender_file(log_path))

  # Set sensitivity: if verbose=TRUE, we record DEBUG level
  if (verbose) {
    logger::log_threshold(logger::DEBUG)
  } else {
    logger::log_threshold(logger::INFO)
  }
}

#' Log Informational Messages
#'
#' Displays a formatted message to the console using `cli` and
#' simultaneously records a clean, non-ANSI version of the message to
#' the log file at the `INFO` level.
#'
#' @param msg Character string. The message to be logged and displayed.
#' @param ... Additional arguments passed to `cli` formatting functions.
#'
#' @export
log_info <- function(msg, ...) {
  caller_env <- parent.frame()
  cli::cli_inform(msg, .envir = caller_env, ...)

  formatted_msg <- cli::format_inline(
    msg,
    .envir = parent.frame()
  )
  clean <- cli::ansi_strip(formatted_msg)

  logger::log_info(clean) # Record to file
}

#' Log Debugging Messages
#'
#' Displays a message to the console and records it to the log file
#' specifically at the `DEBUG` level. Note that the message will only
#' appear in the log file if the logger threshold is set to `DEBUG`.
#'
#' @param msg Character string. The message to be logged and displayed.
#' @param ... Additional arguments passed to `cli` formatting functions.
#'
#' @export
log_debug <- function(msg, ...) {
  caller_env <- parent.frame()
  cli::cli_inform(msg, .envir = caller_env, ...)
  formatted_msg <- cli::format_inline(
    msg,
    .envir = parent.frame()
  )
  clean <- cli::ansi_strip(formatted_msg)
  logger::log_debug(clean) # Record to file ONLY if threshold is DEBUG
}

#' Log Failure Messages and Abort
#'
#' Signals a critical failure by calling `cli::cli_abort()`, which stops
#' execution. The error message is stripped of ANSI formatting and
#' recorded to the log file at the `FAILURE` level.
#'
#' @param msg Character string. The error message.
#' @param ... Additional arguments passed to `cli::cli_abort()`.
#'
#' @export
log_failure <- function(msg, ...) {
  caller_env <- parent.frame()
  cli::cli_abort(msg, .envir = caller_env, ...)
  formatted_msg <- cli::format_inline(
    msg,
    .envir = parent.frame()
  )
  clean <- cli::ansi_strip(formatted_msg)
  logger::log_failure(clean)
}

#' Log Warning Messages
#'
#' Displays a warning to the console and records it to the log file
#' at the `WARN` level.
#'
#' @param msg Character string. The warning message.
#' @param ... Additional arguments passed to `cli::cli_warn()`.
#'
#' @export
log_warning <- function(msg, ...) {
  caller_env <- parent.frame()
  cli::cli_warn(msg, .envir = caller_env, ...)
  formatted_msg <- cli::format_inline(
    msg,
    .envir = parent.frame()
  )
  clean <- cli::ansi_strip(formatted_msg)
  logger::log_warn(clean)
}
