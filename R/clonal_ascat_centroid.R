####################################################################################################
#' This function is an alternative procedure for finding the optimum (psi, rho) pair.
#' This function first finds all the find all the global optima,
#' and then finds the centroid of this set of globla optima.
#' Then we find the global optimum which is nearest to the centroid.
#' (When the set of global optima is convex, we expect the selected optimum to be at the centroid.)
#' @param d A distance matrix
#' @param ref_seg_matrix The corresponding ref seg matrix that belongs to d
#' @param ref_major The corresponding major allele values with d
#' @param ref_minor The corresponding minor allele values with d
#' @param s A segmented BAF/LogR data.frame from \code{get_segment_info}
#' @param dist_choice Some distance metrics require adaptation of the data (i.e. log transform)
#' @param minimise Boolean whether we're minimising or maximising
#' @param new_bounds The rho/psi boundaries between we are searching for a solution. This is a named list with values psi_min, psi_max, rho_min, rho_max
#' @param distancepng String where the sunrise distance plot will be saved
#' @param gamma_param The platform gamma
#' @param siglevel_BAF The level at which BAF becomes significant TODO: this option is no longer used
#' @param maxdist_BAF TODO: this option is no longer used
#' @param siglevel_LogR The p-value at which logR becomes significant when establishing whether a segment should be subclonal
#' @param maxdist_LogR The maximum distance allowed as slack when establishing the significance. This allows for the case when a breakpoint is missed, the segment would then not automatically become subclonal
#' @param allow100percent Boolean whether to allow for a 100"\%" cellularity solution
#' @param uninformative_baf_threshold The threshold above which BAF becomes uninformative
#' @param read_depth TODO: this option is no longer used
#' @return A list with fields optima_info_without_ref and optima_info
#' @export
find_centroid_of_global_minima <- function(
  d, ref_seg_matrix,
  ref_major, ref_minor,
  s, dist_choice, minimise,
  new_bounds, distancepng,
  gamma_param, siglevel_BAF,
  maxdist_BAF, siglevel_LogR,
  maxdist_LogR, allow100percent,
  uninformative_baf_threshold,
  read_depth
) {
  if (!minimise) d <- -d

  # Get global minimum value and grid indices
  gmin <- collapse::fmin(d)
  optima_indices <- which(d == gmin, arr.ind = TRUE)
  nropt <- nrow(optima_indices)

  # Pre-extract numeric grid values from row/col names
  psi_grid <- as.numeric(rownames(d))
  rho_grid <- as.numeric(colnames(d))

  # Map indices to specific psi and rho values for all global optima
  psis <- psi_grid[optima_indices[, 1]]
  rhos <- rho_grid[optima_indices[, 2]]

  # Pre-calculate segment-level constants
  s_length <- s[, "length"]
  s_r <- s[, "r"]
  total_len <- sum(s_length)

  # Calculate the segment-specific term: 2^(r / gamma)
  s_term <- 2^(s_r / gamma_param)

  weighted_s_term <- collapse::fsum(s_term, w = s_length, na.rm = FALSE)
  sum_s_length <- sum(s_length)

  # Calculate the specific ploidy for every global optimum in one vectorized step
  rho_psi_term <- ((1 - rhos) * 2) + (rhos * psis)
  ploidy_vector <- ((2 * rhos - 2) * sum_s_length + (weighted_s_term * rho_psi_term)) / (rhos * total_len)

  # Using collapse::fmedian for C-based speed on the indices
  centre <- c(
    collapse::fmedian(optima_indices[, 1]),
    collapse::fmedian(optima_indices[, 2])
  )

  # Calculate Euclidean distance to the centroid for all points
  row_diffs <- optima_indices[, 1] - centre[1]
  col_diffs <- optima_indices[, 2] - centre[2]
  dists <- (row_diffs^2) + (col_diffs^2)

  best_idx <- which.min(dists)

  # Extract final optimized values
  grid_x <- optima_indices[best_idx, 1]
  grid_y <- optima_indices[best_idx, 2]

  # Format return values
  psi_opt1 <- psi_grid[grid_x]
  rho_opt1 <- min(rho_grid[grid_y], 1)
  ploidy_opt1 <- ploidy_vector[best_idx]
  # Retrieve the reference segment index for the selected grid point
  goodness_of_fit_opt1 <- if (minimise) gmin else -gmin

  ref_seg <- ref_seg_matrix[grid_x, grid_y]
  optima_info_without_ref <- list(
    nropt = nropt,
    psi_opt1 = psi_opt1,
    rho_opt1 = rho_opt1,
    ploidy_opt1 = ploidy_opt1,
    ref_seg = ref_seg,
    goodness_of_fit_opt1 = goodness_of_fit_opt1
  )

  # Handle the logic for determining the final psi/rho based on reference segments
  if (ref_seg == 0) {
    psi_opt1 <- 2
    rho_opt1 <- 1
    ploidy_opt1 <- 2
    goodness_of_fit_opt1 <- 1
  } else {
    ref_segment_info <- get_psi_rho_from_ref_seg(
      ref_seg, s,
      ref_major[grid_x, grid_y],
      ref_minor[grid_x, grid_y],
      gamma_param
    )

    psi_opt1 <- ref_segment_info$psi
    rho_opt1 <- ref_segment_info$rho
    ploidy_opt1 <- ref_segment_info$ploidy

    # Recalculate goodness of fit if a valid rho was found
    if (!is.na(rho_opt1)) {
      distance_info <- calc_distance_clonal(
        s, dist_choice, rho_opt1, psi_opt1, gamma_param,
        read_depth, siglevel_BAF, maxdist_BAF,
        siglevel_LogR, maxdist_LogR, uninformative_baf_threshold
      )
      goodness_of_fit_opt1 <- distance_info$distance_value
    } else {
      goodness_of_fit_opt1 <- Inf
    }
  }

  # Generate the diagnostic sunrise plot if a file path is provided
  if (!is.na(distancepng)) {
    grDevices::png(filename = distancepng, width = 1000, height = 1000, res = 1000 / 7, type = "cairo")
    clonal_findcentroid_plot(minimise, dist_choice, -d, c(psi_opt1), c(rho_opt1), new_bounds)
    grDevices::dev.off()
  }

  # Return the structured results containing both raw and reference-adjusted optima
  return(list(
    optima_info_without_ref = optima_info_without_ref,
    optima_info = list(
      nropt = nropt,
      psi_opt1 = psi_opt1,
      rho_opt1 = rho_opt1,
      ploidy_opt1 = ploidy_opt1,
      ref_seg = ref_seg,
      goodness_of_fit_opt1 = goodness_of_fit_opt1
    )
  ))
}
