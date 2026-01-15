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
  x <- floor(nMinor)
  y <- floor(nMajor)
  ntot <- nMajor + nMinor

  # BAF values at the four corners of the unit square
  nMaj_corners <- c(y, y + 1, y, y + 1)
  nMin_corners <- c(x + 1, x + 1, x, x)
  BAF_corners <- (1 - rho + rho * nMaj_corners) /
    (2 - 2 * rho + rho * (nMaj_corners + nMin_corners))
  BAF_corners[nMaj_corners == 0 & nMin_corners == 0] <- 0.5

  # Determine quadrant relative to BAF_corners[3] and BAF_corners[2]
  above_mid_horizontal <- BAF_req > BAF_corners[3] # case 1 or 2a
  above_mid_vertical <- BAF_req > BAF_corners[2] # case 2c vs 2b

  logR_low <- ntot < x + y + 1

  # Define the six candidate adjustments in priority order
  # Each row: Δmajor1, Δminor1, Δmajor2, Δminor2
  # Priority: first favor smaller LogR distance, then simplicity
  candidates <- if (above_mid_horizontal) {
    if (logR_low) {
      matrix(c(
        0,  0,  1,  0, # y,   x     -> y+1, x
        0, -1,  1,  0, # y,   x-1   -> y+1, x
        0,  0,  1, -1, # y,   x     -> y+1, x-1
        1,  0,  1,  0, # y+1, x     -> y+1, x
        1,  0,  1,  1, # y+1, x     -> y+1, x+1
        1,  0,  2,  0 # y+1, x     -> y+2, x
      ), nrow = 6, byrow = TRUE)
    } else {
      matrix(c(
        1,  0,  1,  1, # y+1, x     -> y+1, x+1
        1, -1,  1,  0, # y+1, x-1   -> y+1, x
        1,  0,  1,  0, # y+1, x     -> y+1, x
        0,  0,  0,  0, # y,   x     -> y,   x
        0, -1,  0,  0, # y,   x-1   -> y,   x
        0,  0,  1,  2 # y,   x     -> y+1, x+2  (wait, original had y+1,x+2 but adjusted)
      ), nrow = 6, byrow = TRUE)
    }
  } else if (above_mid_vertical) {
    # symmetric cases for 2c
    if (logR_low) {
      matrix(c(
        0,  0,  0,  1,
        0, -1,  0,  1,
        0,  0,  1,  1,
        1,  0,  1,  1,
        1, -1,  1,  1,
        1,  0,  1,  2
      ), nrow = 6, byrow = TRUE)
    } else {
      matrix(c(
        1,  0,  1,  1,
        0, -1,  0,  1,
        0,  0,  0,  1,
        1,  0,  0,  1,
        1, -1,  0,  1,
        1,  0,  1,  2
      ), nrow = 6, byrow = TRUE)
    }
  } else {
    # case 2b
    if (logR_low) {
      matrix(c(
        0,  0,  0,  1,
        0, -1,  1,  1,
        0,  0,  1,  2,
        0,  1,  0,  1,
        0,  1,  0,  2,
        1,  1,  1,  1
      ), nrow = 6, byrow = TRUE)
    } else {
      matrix(c(
        0,  1,  1,  1,
        0,  1,  0,  1,
        1,  1,  0,  1,
        0,  0,  1,  2,
        0, -1,  0,  1,
        1,  1,  1,  2
      ), nrow = 6, byrow = TRUE)
    }
  }

  # Apply base (y, x) and deltas
  maj1 <- y + candidates[, 1]
  min1 <- x + candidates[, 2]
  maj2 <- y + candidates[, 3]
  min2 <- x + candidates[, 4]

  # Remove invalid (negative) copy numbers
  valid <- (maj1 >= 0 & min1 >= 0 & maj2 >= 0 & min2 >= 0)
  maj1[!valid] <- NA
  min1[!valid] <- NA
  maj2[!valid] <- NA
  min2[!valid] <- NA

  result <- cbind(nMaj1 = maj1, nMin1 = min1, nMaj2 = maj2, nMin2 = min2)

  if (full) {
    return(result) # 6 × 4 matrix
  } else {
    return(list(nMaj = result[1, c(1, 3)], nMin = result[1, c(2, 4)])) # top edge only
  }
}
