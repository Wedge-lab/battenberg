#' @export
log_setup <- function(log_path, verbose = FALSE) {
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

#' @export
log_info <- function(msg, ...) {
  cli::cli_inform(msg, ...) # High-level UI for the human
  formatted_msg <- cli::cli_format_method(
    cli::cli_text(msg),
    .envir = parent.frame()
  )
  clean <- cli::ansi_strip(formatted_msg)

  logger::log_info(clean) # Record to file
}

#' @export
log_debug <- function(msg, ...) {
  cli::cli_inform(msg, ...)
  formatted_msg <- cli::cli_format_method(
    cli::cli_text(msg),
    .envir = parent.frame()
  )
  clean <- cli::ansi_strip(formatted_msg)
  logger::log_debug(clean) # Record to file ONLY if threshold is DEBUG
}

#' @export
log_failure <- function(msg, ...) {
  cli::cli_abort(msg, ...)
  formatted_msg <- cli::cli_format_method(
    cli::cli_text(msg),
    .envir = parent.frame()
  )
  clean <- cli::ansi_strip(formatted_msg)
  logger::log_failure(clean) # Record to file ONLY if threshold is DEBUG
}
