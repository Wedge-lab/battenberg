#' Run code in parallel or serial based on debug status
#'
#' A helper function to abstract the pattern of switching between parallel
#' execution via foreach and serial execution via lapply.
#'
#' @param iterator A vector or list to iterate over (e.g., seq_along(x)).
#' @param func A function to apply to each element of the iterator.
#' @param libs Path to library paths for workers.
#'
#' @return A list of results from the applied function.
#' @keywords internal
run_parallel_or_serial <- function(iterator, func, libs) {
  if (length(iterator) == 0) {
    return(list())
  }

  # Set up foreach to use the registered backend
  `%dopar%` <- foreach::`%dopar%`

  foreach::foreach(i = iterator) %dopar% {
    .libPaths(libs)

    # Wrap in calling handler to capture more context on failure
    # This remains in parallel but gives us more info if it crashes
    withCallingHandlers(
      {
        func(i)
      },
      error = function(e) {
        # In parallel workers, stdout/stderr are often captured or redirected.
        # By using cat() here, it will go to the cluster's outfile,
        # which we set to the empty string (master's stdout) in battenberg.R.
        msg <- sprintf("!!! BATTENBERG ERROR IN PARALLEL WORKER NODE %s !!!\nMessage: %s\nStack Trace:", i, conditionMessage(e))
        calls <- sys.calls()
        for (j in rev(seq_along(calls))) {
          msg <- paste(msg, sprintf("%d: %s", j, deparse(calls[[j]])), sep = "\n")
        }
        msg <- paste(msg, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", sep = "\n")
        log_failure("{msg}")
      }
    )
  }
}
