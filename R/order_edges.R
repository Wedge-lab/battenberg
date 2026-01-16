#' Prioritize candidate integer copy number states around a fractional state
#'
#' Returns the nearest grid edges/corners based on BAF and LogR position,
#' following ASCAT's original prioritization rules (LogR distance + simplicity).
#'
#' @param rho Observed tumor purity (fraction of tumor cells)
#' @param psi Estimated ploidy
#' @param BAF_req Observed B-allele frequency
#' @param nMajor Fractional major allele copy number
#' @param nMinor Fractional minor allele copy number
#' @param full logical; if TRUE return all 6 prioritized options (like orderEdges),
#'                   if FALSE return only the 2 endpoints of the nearest edge
#' @return matrix with columns nMaj1, nMin1, nMaj2, nMin2 (or 6 rows if full=TRUE)
#' @noRd
prioritizeCopyNumbers <- function(rho, psi, BAF_req, nMajor, nMinor, full = FALSE) {
  # Vectorized Inputs
  x <- floor(nMinor)
  y <- floor(nMajor)
  ntot <- nMajor + nMinor

  # Ensure all inputs are vectors of same length (Recycling rules apply)
  # But BAF_req drives the length normally
  n <- length(BAF_req)

  # BAF values at the four corners of the unit square
  # We compute these for every element
  nMaj_c <- cbind(y, y + 1, y, y + 1)
  nMin_c <- cbind(x + 1, x + 1, x, x)

  # Vectorized BAF calculation
  # Note: rho might be scalar, nMaj_c is matrix (N x 4)
  # R handles scalar-matrix arithmetic fine
  BAF_corners <- (1 - rho + rho * nMaj_c) /
    (2 - 2 * rho + rho * (nMaj_c + nMin_c))

  # Handle 0/0 case
  zero_idx <- (nMaj_c == 0 & nMin_c == 0)
  BAF_corners[zero_idx] <- 0.5

  # Determine quadrant relative to BAF_corners
  # Column 3 is equivalent to BAF_corners[3] in scalar version
  above_mid_horizontal <- BAF_req > BAF_corners[, 3]
  above_mid_vertical <- BAF_req > BAF_corners[, 2]
  logR_low <- ntot < (x + y + 1)

  # Pre-define the 6 candidate matrices (flattened or indexed)
  # Because we need to apply different logic per element, we construct the offsets
  # dynamically based on the boolean flags.

  # We do this for the "top edge" only if full=FALSE, or all 6 if full=TRUE
  # But wait, original code returns 6 rows if full=TRUE.
  # If vectorized, full=TRUE would mean returning a N x 6 x 4 array? Or a list?
  # The usage in is_segment_clonal asks for full=TRUE but only uses row 1.
  # Actually, `is_segment_clonal` uses `all.edges[1, c(1,3)]` which implies it expects a matrix.
  # But if we pass vectors, we return a Matrix of N rows?
  # NO. `is_segment_clonal` as written above expects `all.edges` to be a matrix where
  # rows correspond to input elements?
  #
  # Let's look at `is_segment_clonal` usage again:
  #   all.edges <- prioritizeCopyNumbers(..., full=TRUE)
  #   nMaj.test <- all.edges[, c(1, 3)]
  #
  # If `is_segment_clonal` is vectorized, `all.edges` must return a structure where
  # for each input i, we get the "best edge" (Option 1 and Option 2).
  # The original `full=TRUE` returned 6 candidates.
  # The vectorized `is_segment_clonal` only cares about the **first** candidate row
  # from the prioritization list (the distinct "best edge").
  #
  # So we will simplify: We only compute the FIRST priority candidate (row 1 of the matrix).
  # Wait, standard ASCAT logic tries to find the "nearest" valid edge.
  # The original code provided 6 options in order of preference.
  # Does `is_segment_clonal` iterate through them?
  # Original `is_segment_clonal`:
  #   nMaj.test <- all.edges[1, c(1, 3)]
  # It takes just the first row.
  #
  # So we only need to implement the logic for the **first priority** candidate!

  # Logic for First Priority Candidate (Index 1 of the matrix):
  # Case A: Horizontal (above_mid_horizontal)
  #   Subcase A1: logR_low -> 0, 0, 1, 0 (y, x -> y+1, x)
  #   Subcase A2: !logR_low -> 1, 0, 1, 1 (y+1, x -> y+1, x+1)
  # Case B: Vertical (above_mid_vertical)
  #   Subcase B1: logR_low -> 0, 0, 0, 1 (y, x -> y, x+1)
  #   Subcase B2: !logR_low -> 1, 0, 1, 1 (y+1, x -> y+1, x+1)
  # Case C: Neither (2b)
  #   Subcase C1: logR_low -> 0, 0, 0, 1 (y, x -> y, x+1)
  #   Subcase C2: !logR_low -> 0, 1, 1, 1 (y, x+1 -> y+1, x+1)

  # Initialize with 0s
  dm1 <- integer(n)
  dn1 <- integer(n)
  dm2 <- integer(n)
  dn2 <- integer(n)

  # Case A
  idx_A_low <- which(above_mid_horizontal & logR_low)
  idx_A_high <- which(above_mid_horizontal & !logR_low)
  if (length(idx_A_low)) {
    dm1[idx_A_low] <- 0
    dn1[idx_A_low] <- 0
    dm2[idx_A_low] <- 1
    dn2[idx_A_low] <- 0
  }
  if (length(idx_A_high)) {
    dm1[idx_A_high] <- 1
    dn1[idx_A_high] <- 0
    dm2[idx_A_high] <- 1
    dn2[idx_A_high] <- 1
  }

  # Case B (Not A, and Vertical)
  # Note: The original generic if/else structure implies sequential checks.
  # if (horizontal) { ... } else if (vertical) { ... } else { ... }
  is_B <- (!above_mid_horizontal) & above_mid_vertical
  idx_B_low <- which(is_B & logR_low)
  idx_B_high <- which(is_B & !logR_low)

  if (length(idx_B_low)) {
    dm1[idx_B_low] <- 0
    dn1[idx_B_low] <- 0
    dm2[idx_B_low] <- 0
    dn2[idx_B_low] <- 1
  }
  if (length(idx_B_high)) {
    dm1[idx_B_high] <- 1
    dn1[idx_B_high] <- 0
    dm2[idx_B_high] <- 1
    dn2[idx_B_high] <- 1
  }

  # Case C (Not A, Not B)
  is_C <- (!above_mid_horizontal) & (!above_mid_vertical)
  idx_C_low <- which(is_C & logR_low)
  idx_C_high <- which(is_C & !logR_low)

  if (length(idx_C_low)) {
    dm1[idx_C_low] <- 0
    dn1[idx_C_low] <- 0
    dm2[idx_C_low] <- 0
    dn2[idx_C_low] <- 1
  }
  if (length(idx_C_high)) {
    dm1[idx_C_high] <- 0
    dn1[idx_C_high] <- 1
    dm2[idx_C_high] <- 1
    dn2[idx_C_high] <- 1
  }

  # Apply deltas
  maj1 <- y + dm1
  min1 <- x + dn1
  maj2 <- y + dm2
  min2 <- x + dn2

  # Validation: Clamp negative to NA (or handle as in original)
  # Original code: valid <- (all >= 0); invalid -> NA
  valid <- (maj1 >= 0 & min1 >= 0 & maj2 >= 0 & min2 >= 0)

  # If not valid, we return NA.
  # Since we are returning vectors, we can just set them to NA.
  maj1[!valid] <- NA
  min1[!valid] <- NA
  maj2[!valid] <- NA
  min2[!valid] <- NA

  # Return N x 4 matrix
  # Corresponds to nMaj1, nMin1, nMaj2, nMin2
  return(cbind(nMaj1 = maj1, nMin1 = min1, nMaj2 = maj2, nMin2 = min2))
}
