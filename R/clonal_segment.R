#' This function decides if a segment is "clonal" (= TRUE) or not (= FALSE).
#' (The alternative hypothesis is that the tumour genome segment in question exhibits "sub-clonal" variation.)
#' We test the integer solutions for all 4 corners. Also, along side the hypothesis test for the BAF.
#' We use a decision rule based on LogR (we could use a hypothesis test which takes account of the variance in LogR, or a fixed “tolerance”).
#' If the null hypothesis is accepted for at least one corner, then we accept that
#' the tumour genome segment in question is "clonal".
#' @noRd
is_segment_clonal <- function(
  LogR, BAF_req, BAF_length, BAF_size, BAF_mean, BAF_sd,
  rho, psi, gamma_param, siglevel_BAF, maxdist_BAF
) {
  # Handle NAs in LogR efficiently
  # If LogR is a vector, we modify it in place
  LogR[is.na(LogR)] <- 0

  # Pre-calculate shared terms
  factor <- 2^(LogR / gamma_param)
  term_base <- (rho - 1)
  term_psi <- ((1 - rho) * 2 + rho * psi)

  nA <- (term_base - (BAF_req - 1) * factor * term_psi) / rho
  nB <- (term_base + BAF_req * factor * term_psi) / rho

  nMajor <- pmax(nA, nB, na.rm = TRUE)
  nMinor <- pmin(nA, nB, na.rm = TRUE)

  # Check validation logic (Vectorized)
  nMajor.saved <- nMajor

  # Validation logic for negative nMinor
  neg_idx <- which(nMinor < 0)
  if (length(neg_idx) > 0) {
    b_req_sub <- BAF_req[neg_idx]

    # Identify which ones are BAF_req == 1
    is_one <- abs(b_req_sub - 1) < 1e-9

    # Case 1: BAF == 1 -> Major = 1000
    nMajor[neg_idx[is_one]] <- 1000

    # Case 2: BAF != 1 -> Recalculate Major
    not_one <- neg_idx[!is_one]
    if (length(not_one) > 0) {
      val <- nMajor[not_one] + BAF_req[not_one] * (0.01 - nMinor[not_one]) / (1 - BAF_req[not_one])
      # Clamp to 1000 if negative
      val[val < 0] <- 1000
      nMajor[not_one] <- val
    }

    nMinor[neg_idx] <- 0.01
  }

  # prioritizeCopyNumbers is now vectorized (assumed - we will update it next)
  all.edges <- prioritizeCopyNumbers(
    rho = rho, psi = psi, BAF_req = BAF_req,
    nMajor = nMajor, nMinor = nMinor, full = TRUE
  )

  # Columns: 1=nM1, 2=nm1, 3=nM2, 4=nm2
  nMaj.test <- all.edges[, c(1, 3), drop = FALSE]
  nMin.test <- all.edges[, c(2, 4), drop = FALSE]

  # Calculate levels for both options (Option 1 and Option 2)
  calc_baf <- function(nM, nm) {
    num <- 1 - rho + rho * nM
    den <- 2 - 2 * rho + rho * (nM + nm)
    lev <- num / den
    lev[nM == 0 & nm == 0] <- 0.5
    lev
  }

  lev1 <- calc_baf(nMaj.test[, 1], nMin.test[, 1])
  lev2 <- calc_baf(nMaj.test[, 2], nMin.test[, 2])

  dist1 <- abs(lev1 - BAF_req)
  dist2 <- abs(lev2 - BAF_req)

  # Vectorized choice of best index
  choose_2 <- dist2 < dist1

  best_nMaj <- ifelse(choose_2, nMaj.test[, 2], nMaj.test[, 1])
  best_nMin <- ifelse(choose_2, nMin.test[, 2], nMin.test[, 1])
  best_level <- ifelse(choose_2, lev2, lev1)

  # P-value calculation
  # Handle BAF_sd == 0 case
  pval <- numeric(length(BAF_req))
  valid_sd <- BAF_sd > 0

  if (any(valid_sd)) {
    # Assuming calc_Pvalue_t_twotailed is vectorized
    pval[valid_sd] <- calc_Pvalue_t_twotailed(
      BAF_size[valid_sd], BAF_req[valid_sd],
      BAF_sd[valid_sd], best_level[valid_sd], maxdist_BAF
    )
  }
  # SD == 0 stays 0

  balanced <- (best_nMaj == best_nMin)

  # Clonal decision
  is_clonal <- (pval > siglevel_BAF)

  # Stability check (Vectorized)
  unstable <- (nMajor - nMajor.saved) >= 1
  is_clonal[unstable] <- FALSE

  return(list(
    is_clonal = is_clonal,
    balanced = balanced,
    nMaj = best_nMaj,
    nMin = best_nMin
  ))
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
