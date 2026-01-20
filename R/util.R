########################################################################################
# Other
########################################################################################
#' Check if a file exists, if it doesn't, exit non-clean
#' @noRd
assert_file_exists <- function(filename) {
  if (!file.exists(filename)) {
    log_failure("Supplied file does not exist: {filename}")
    quit(save = "no", status = 1)
  }
}
