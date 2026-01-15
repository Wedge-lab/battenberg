#' This function decides if a segment is "clonal" (= TRUE) or not (= FALSE).
#' (The alternative hypothesis is that the tumour genome segment in question exhibits "sub-clonal" variation.)
#' We test the integer solutions for all 4 corners. Also, along side the hypothesis test for the BAF.
#' We use a decision rule based on LogR (we could use a hypothesis test which takes account of the variance in LogR, or a fixed “tolerance”).
#' If the null hypothesis is accepted for at least one corner, then we accept that
#' the tumour genome segment in question is "clonal".
#' @noRd
is_segment_clonal <- function(
  LogR,
  BAF_req,
  BAF_length,
  BAF_size,
  BAF_mean,
  BAF_sd,
  rho,
  psi,
  gamma_param,
  siglevel_BAF,
  maxdist_BAF
) {
  # if we don't have a value for LogR, fill in 0
  if (is.na(LogR)) {
    LogR <- 0
  }

  nA <- (rho - 1 - (BAF_req - 1) * 2^(LogR / gamma_param) * ((1 - rho) * 2 + rho * psi)) / rho
  nB <- (rho - 1 + BAF_req * 2^(LogR / gamma_param) * ((1 - rho) * 2 + rho * psi)) / rho

  nMajor <- max(nA, nB, na.rm = TRUE)
  nMinor <- min(nA, nB, na.rm = TRUE)

  # check for big shifts in nMajor - if there's a big shift, we shouldn't trust a clonal call
  nMajor.saved <- nMajor

  # DCW - increase nMajor and nMinor together, to avoid impossible combinations (with negative subclonal fractions)
  if (nMinor < 0) {
    if (BAF_req == 1) {
      # avoid calling infinite copy number
      nMajor <- 1000
    } else {
      nMajor <- nMajor + BAF_req * (0.01 - nMinor) / (1 - BAF_req)
      if (nMajor < 0) nMajor <- 1000
    }
    nMinor <- 0.01
  }

  # note that these are sorted in the order of ascending BAF:
  nMaj <- c(floor(nMajor), ceiling(nMajor), floor(nMajor), ceiling(nMajor))
  nMin <- c(ceiling(nMinor), ceiling(nMinor), floor(nMinor), floor(nMinor))

  BAF_levels <- (1 - rho + rho * nMaj) / (2 - 2 * rho + rho * (nMaj + nMin))
  # problem if rho=1 and nMaj=0 and nMin=0
  BAF_levels[nMaj == 0 & nMin == 0] <- 0.5

  # DCW - just test corners on the nearest edge to determine clonality
  # If the segment is called as subclonal, this is the edge that will be used to determine the subclonal proportions that are reported first
  all.edges <- prioritizeCopyNumbers(
    rho = rho,
    psi = psi,
    BAF_req = BAF_req, # The observed BAF value for this segment
    nMajor = nMajor,
    nMinor = nMinor,
    full = TRUE
  )

  nMaj.test <- all.edges[1, c(1, 3)]
  nMin.test <- all.edges[1, c(2, 4)]
  test.BAF_levels <- (1 - rho + rho * nMaj.test) / (2 - 2 * rho + rho * (nMaj.test + nMin.test))
  # problem if rho=1 and nMaj=0 and nMin=0
  test.BAF_levels[nMaj.test == 0 & nMin.test == 0] <- 0.5
  whichclosestlevel.test <- which.min(abs(test.BAF_levels - BAF_req))

  # problem caused by segments with constant BAF (usually 1 or 2)
  if (BAF_sd == 0) {
    pval <- 0
  } else {
    pval <- calc_Pvalue_t_twotailed(BAF_size, BAF_req, BAF_sd, test.BAF_levels[whichclosestlevel.test], maxdist_BAF)
  }

  balanced <- nMaj.test[whichclosestlevel.test] == nMin.test[whichclosestlevel.test]

  is_clonal <- (pval > siglevel_BAF)
  # check for big shifts in nMajor - if there's a big shift, we shouldn't trust a clonal call
  # This is particularly problematic for very high cellularity samples, like some of the ovarian samples
  is_clonal <- (pval > siglevel_BAF & nMajor - nMajor.saved < 1)

  return(list(is_clonal = is_clonal, balanced = balanced, nMaj.test = nMaj.test[whichclosestlevel.test], nMin.test = nMin.test[whichclosestlevel.test]))
}


#' Helper function to find new rho and psi boundaries given a current optimum pair.
#' @noRd
get_new_bounds <- function(input_optimum_pair, initial_bounds) {
  # Define the window sizes (half-ranges)
  psi_half <- 0.05 * (initial_bounds$psi_max - initial_bounds$psi_min)
  rho_half <- 0.05 * input_optimum_pair$rho

  # Calculate raw windows
  psi_bounds <- c(input_optimum_pair$psi - psi_half, input_optimum_pair$psi + psi_half)
  rho_bounds <- c(input_optimum_pair$rho - rho_half, input_optimum_pair$rho + rho_half)

  # Clamp the windows to ensure they stay within initial boundaries
  # If the window hits the bottom, shift it up; if it hits the top, shift it down.
  adjust_bounds <- function(bounds, start, end) {
    range_val <- bounds[2] - bounds[1]
    low <- max(start, min(end - range_val, bounds[1]))
    high <- min(end, max(start + range_val, bounds[2]))
    return(c(low, high))
  }

  psi_final <- adjust_bounds(
    psi_bounds, initial_bounds$psi_min,
    initial_bounds$psi_max
  )
  rho_final <- adjust_bounds(
    rho_bounds, initial_bounds$rho_min,
    initial_bounds$rho_max
  )

  return(list(
    psi_min = psi_final[1], psi_max = psi_final[2],
    rho_min = rho_final[1], rho_max = rho_final[2]
  ))
}
