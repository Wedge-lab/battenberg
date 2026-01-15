#' Run code in parallel or serial based on debug status
#'
#' A helper function to abstract the pattern of switching between parallel
#' execution via foreach and serial execution via lapply.
#'
#' @param iterator A vector or list to iterate over (e.g., seq_along(x)).
#' @param func A function to apply to each element of the iterator.
#' @param debug Logical; if TRUE, uses lapply for easier debugging and
#'   tracebacks. If FALSE, uses foreach with the %dopar% operator.
#'
#' @return A list of results from the applied function.
#' @keywords internal
run_parallel_or_serial <- function(iterator, func, debug, libs) {
  if (length(iterator) == 0) {
    log_info("Warning: {iterator} is empty")
    return(list())
  }
  if (debug) {
    # Sequential execution for easier debugging/tracebacks
    lapply(iterator, func)
  } else {
    # Parallel execution
    `%dopar%` <- foreach::`%dopar%`
    foreach::foreach(i = iterator) %dopar% {
      .libPaths(libs)
      func(i)
    }
  }
}
