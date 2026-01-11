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

#' A modified ASCAT main function to fit Battenberg
#'
#' This function returns an initial rho and psi estimate for a clonal copy number fit. It uses an internal distance metric to create a distance matrix.
#' Using that matrix it will search for a rho and psi combination that yields the least heavy penalty.
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param baf (unsegmented) B Allele Frequency, in genomic sequence (all probes), with probe IDs
#' @param lrrsegmented log R, segmented, in genomic sequence (all probes), with probe IDs
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param chromosomes a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#' @param dist_choice The distance metric to be used internally to penalise a copy number solution
#' @param distancepng if NA: distance is plotted, if filename is given, the plot is written to a .png file (Default NA)
#' @param copynumberprofilespng if NA: possible copy number profiles are plotted, if filename is given, the plot is written to a .png file (Default NA)
#' @param nonroundedprofilepng if NA: copy number profile before rounding is plotted (total copy number as well as the copy number of the minor allele), if filename is given, the plot is written to a .png file (Default NA)
#' @param cnaStatusFile File where the copy number profile status is written to. This contains either the message "No suitable copy number solution found" or "X copy number solutions found" (Default copynumber_solution_status.txt)
#' @param gamma technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 "\%" aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays) (Default 0.55)
#' @param allow100percent A boolean whether to allow a 100"\%" cellularity solution
#' @param reliabilityFile String to where fit reliabilty information should be written. This file contains backtransformed BAF and LogR values for segments using the fitted copy number profile (Default NA)
#' @param min_ploidy The minimum ploidy to consider (Default 1.6)
#' @param max_ploidy The maximum ploidy to consider (Default 4.8)
#' @param min_rho The minimum cellularity to consider (Default 0.1)
#' @param max_rho The maximum cellularity to consider (Default 1.0)
#' @param min_goodness The minimum goodness of fit for a solution to have to be considered (Default 63)
#' @param uninformative_baf_threshold The threshold beyond which BAF becomes uninformative (Default 0.51)
#' @param chr_names A vector with chromosome names used for plotting
#' @param analysis A String representing the type of analysis to be run, this determines whether the distance figure is produced (Default paired)
#' @return A list with fields psi, rho and ploidy
#' @export
# the limit on rho is lenient and may lead to spurious solutions
runASCAT <- function(
  lrr, baf, lrrsegmented,
  bafsegmented, chromosomes,
  dist_choice, distancepng = NA,
  copynumberprofilespng = NA,
  nonroundedprofilepng = NA,
  cnaStatusFile = "copynumber_solution_status.txt",
  gamma = 0.55, allow100percent,
  reliabilityFile = NA, min_ploidy = 1.6,
  max_ploidy = 4.8, min_rho = 0.1,
  max_rho = 1.0, min_goodness = 63,
  uninformative_baf_threshold = 0.51,
  chr_names, analysis = "paired"
) {
  # Setup inputs and segments
  ch <- chromosomes
  b <- bafsegmented
  r <- lrrsegmented[names(bafsegmented)]

  # Adapt the rho/psi boundaries
  dist_min_psi <- max(min_ploidy - 0.6, 0)
  dist_max_psi <- max_ploidy + 0.6
  dist_min_rho <- max(min_rho - 0.03, 0.05)
  dist_max_rho <- max_rho + 0.03

  s <- ASCAT::make_segments(r, b)
  dist_matrix_info <- create_distance_matrix(
    s, dist_choice, gamma,
    uninformative_baf_threshold = uninformative_baf_threshold,
    min_psi = dist_min_psi,
    max_psi = dist_max_psi,
    min_rho = dist_min_rho,
    max_rho = dist_max_rho
  )
  d <- dist_matrix_info$distance_matrix
  minimise <- dist_matrix_info$minimise

  # Calculate theoretical max distance for goodness of fit
  TheoretMaxdist <- sum(rep(0.25, dim(s)[1]) * s[, "length"], na.rm = TRUE)
  total_len <- sum(s[, "length"])

  # Ensure we are always searching for a minimum
  if (!minimise) d <- -d

  # VECTORIZED LOCAL MINIMA SEARCH (Pixel-perfect replacement for 7x7 loop)
  nr <- nrow(d)
  nc <- ncol(d)
  is_local_min <- matrix(TRUE, nrow = nr, ncol = nc)

  # Constrain search to the interior to match 4:(dim-3) logic
  row_range <- 4:(nr - 3)
  col_range <- 4:(nc - 3)

  # Check every neighbor in the 7x7 window (48 neighbors)
  for (dx in -3:3) {
    for (dy in -3:3) {
      if (dx == 0 && dy == 0) next
      is_local_min[row_range, col_range] <- is_local_min[row_range, col_range] &
        (d[row_range, col_range] < d[row_range + dx, col_range + dy])
    }
  }

  # Zero out the margins to match original loop boundaries
  is_local_min[-row_range, ] <- FALSE
  is_local_min[, -col_range] <- FALSE

  # Extraction helper to process candidates
  evaluate_candidates <- function(indices, current_d) {
    if (nrow(indices) == 0) {
      return(NULL)
    }

    # Pre-calculate segment masks for efficiency
    is_not_balanced <- s[, "b"] != 0.5
    weight_unbalanced <- sum(s[, "length"] * is_not_balanced)

    results <- apply(indices, 1, function(idx) {
      i <- idx[1]
      j <- idx[2]
      m <- current_d[i, j]
      psi <- as.numeric(rownames(current_d)[i])
      rho <- as.numeric(colnames(current_d)[j])

      # Copy number algebra
      common_term <- 2^(s[, "r"] / gamma) * ((1 - rho) * 2 + rho * psi)
      nA <- (rho - 1 - (s[, "b"] - 1) * common_term) / rho
      nB <- (rho - 1 + s[, "b"] * common_term) / rho

      ploidy <- sum((nA + nB) * s[, "length"]) / total_len

      # Biological viability checks
      is_nA_zero <- round(nA) == 0
      is_nB_zero <- round(nB) == 0
      percentzero <- (sum(is_nA_zero * s[, "length"]) + sum(is_nB_zero * s[, "length"])) / total_len
      perczeroAbb <- (sum(is_nA_zero * s[, "length"] * is_not_balanced) + sum(is_nB_zero * s[, "length"] * is_not_balanced)) / weight_unbalanced
      if (is.na(perczeroAbb)) perczeroAbb <- 0

      # Goodness of fit calculation
      fit <- if (minimise) (1 - m / TheoretMaxdist) * 100 else -m / TheoretMaxdist * 100

      # Return data if it meets primary constraints (percentzero checks applied later if allow100percent is used)
      return(list(m = m, i = i, j = j, ploidy = ploidy, fit = fit, pz = percentzero, pza = perczeroAbb, rho = rho, psi = psi))
    })
    return(results)
  }

  # First pass: find optima meeting the percentzero conditions
  opt_indices <- which(is_local_min, arr.ind = TRUE)
  candidates <- evaluate_candidates(opt_indices, d)

  # Filtering based on standard Battenberg criteria
  valid_optima <- Filter(function(x) {
    x$ploidy >= min_ploidy && x$ploidy <= max_ploidy &&
      x$rho >= min_rho && x$fit >= min_goodness &&
      (x$pz > 0.01 || x$pza > 0.1)
  }, candidates)

  # Second pass: If allow100percent is TRUE and no solutions found, relax constraints
  if (allow100percent && length(valid_optima) == 0) {
    # Penalize cellularity > 1 as per original code
    cold_idx <- which(as.numeric(colnames(d)) > 1)
    d[, cold_idx] <- 1e20

    # Re-evaluate all local minima with relaxed biological constraints
    valid_optima <- Filter(function(x) {
      x$ploidy > min_ploidy && x$ploidy < max_ploidy &&
        x$rho >= min_rho && x$fit >= min_goodness
    }, candidates)
  }

  # Process the winning solution
  nropt <- length(valid_optima)
  psi_opt1_plot <- vector(mode = "numeric")
  rho_opt1_plot <- vector(mode = "numeric")

  if (nropt > 0) {
    data.table::fwrite(paste(nropt, " copy number solutions found", sep = ""), file = cnaStatusFile, quote = FALSE, col_names = FALSE, row.names = FALSE)

    # Find the global minimum among the local optima
    all_m <- sapply(valid_optima, function(x) x$m)
    optlim <- min(all_m)

    # Extract ties for plotting and set the final result
    for (opt in valid_optima) {
      if (opt$m == optlim) {
        psi_opt1 <- opt$psi
        rho_opt1 <- min(opt$rho, 1)
        ploidy_opt1 <- opt$ploidy
        goodness_of_fit_opt1 <- opt$fit

        psi_opt1_plot <- c(psi_opt1_plot, psi_opt1)
        rho_opt1_plot <- c(rho_opt1_plot, rho_opt1)
      }
    }
  } else {
    data.table::fwrite("no copy number solutions found", file = cnaStatusFile, quote = FALSE, col_names = FALSE, row.names = FALSE)
    print("No suitable copy number solution found")
    psi <- ploidy <- rho <- NA
    psi_opt1_plot <- rho_opt1_plot <- -1
  }

  # Plotting Sunrise (if paired)
  if (analysis == "paired") {
    if (!is.na(distancepng)) {
      grDevices::png(filename = distancepng, width = 1000, height = 1000, res = 1000 / 7, type = "cairo")
      ASCAT::ascat.plotSunrise(-d, psi_opt1_plot, rho_opt1_plot, minimise)
      grDevices::dev.off()
    }
  }

  # Final calculations for the best solution
  if (nropt > 0) {
    rho <- rho_opt1
    psi <- psi_opt1
    ploidy <- ploidy_opt1

    # Full genomic fit
    nAfull <- (rho - 1 - (b - 1) * 2^(r / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
    nBfull <- (rho - 1 + b * 2^(r / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
    nA <- pmax(round(nAfull), 0)
    nB <- pmax(round(nBfull), 0)

    # Reliability and back-transformation
    rBT <- gamma * log((rho * (nA + nB) + (1 - rho) * 2) / ((1 - rho) * 2 + rho * psi), 2)
    bBT <- (1 - rho + rho * nB) / (2 - 2 * rho + rho * (nA + nB))

    if (!is.na(reliabilityFile)) {
      data.table::fwrite(data.frame(segmentedBAF = b, backTransformedBAF = bBT, segmentedR = r, backTransformedR = rBT, nA = nA, nB = nB, nAfull = nAfull, nBfull = nBfull), reliabilityFile, sep = ",", row.names = FALSE)
    }

    # Generate Profile Plots
    if (!is.na(copynumberprofilespng)) {
      grDevices::png(
        filename = copynumberprofilespng,
        width = 2000, height = 500,
        res = 200, type = "cairo"
      )
      ASCAT::ascat.plotAscatProfile(
        n1all = nA, n2all = nB, heteroprobes = TRUE,
        ploidy = ploidy_opt1, rho = rho_opt1,
        goodness_of_fit = goodness_of_fit_opt1,
        nonaberrant = FALSE, ch = ch,
        lrr = lrr, bafsegmented = bafsegmented,
        chrs = chr_names
      )
      grDevices::dev.off()
    }

    if (!is.na(nonroundedprofilepng)) {
      grDevices::png(
        filename = nonroundedprofilepng,
        width = 2000, height = 500,
        res = 200, type = "cairo"
      )
      ASCAT::ascat.plotNonRounded(
        ploidy = ploidy_opt1, rho = rho_opt1,
        goodness_of_fit = goodness_of_fit_opt1,
        nonaberrant = FALSE, nAfull = nAfull,
        nBfull = nBfull, bafsegmented = bafsegmented,
        ch = ch, lrr = lrr, chrs = chr_names
      )
      grDevices::dev.off()
    }
  }

  return(list(psi = psi, rho = rho, ploidy = ploidy))
}

####################################################################################################
#' ASCAT like function to obtain a clonal copy number profile
#'
#' This function takes an initial optimum rho/psi pair and uses
#' an internal distance metric to calculate a score for each rho/psi pair allowed.
#' The solution with the best score is then taken to obtain a global copy number
#' profile. This function performs both a grid search and tries to find a reference
#' segment, but the grid search result is always used for now.
#' @param lrr (unsegmented) log R, in genomic sequence (all probes), with probe IDs
#' @param baf (unsegmented) B Allele Frequency, in genomic sequence (all probes), with probe IDs
#' @param lrrsegmented log R, segmented, in genomic sequence (all probes), with probe IDs
#' @param bafsegmented B Allele Frequency, segmented, in genomic sequence (only probes heterozygous in germline), with probe IDs
#' @param chromosomes a list containing c vectors, where c is the number of chromosomes and every vector contains all probe numbers per chromosome
#' @param segBAF_table Segmented BAF data.frame from \code{get_segment_info}
#' @param input_optimum_pair A list containing fields for rho, psi and ploidy, as is output from \code{runASCAT}
#' @param dist_choice The distance metric to be used internally to penalise a copy number solution
#' @param distancepng if NA: distance is plotted, if filename is given, the plot is written to a .png file (Default NA)
#' @param copynumberprofilespng if NA: possible copy number profiles are plotted, if filename is given, the plot is written to a .png file (Default NA)
#' @param nonroundedprofilepng if NA: copy number profile before rounding is plotted (total copy number as well as the copy number of the minor allele), if filename is given, the plot is written to a .png file (Default NA)
#' @param gamma_param technology parameter, compaction of Log R profiles (expected decrease in case of deletion in diploid sample, 100 "\%" aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays) (Default 0.55)
#' @param read_depth TODO: unused parameter that should be removed
#' @param uninformative_baf_threshold The threshold beyond which BAF becomes uninformative
#' @param allow100percent A boolean whether to allow a 100"\%" cellularity solution
#' @param reliabilityFile String to where fit reliabilty information should be written. This file contains backtransformed BAF and LogR values for segments using the fitted copy number profile (Default NA)
#' @param psi_min_initial Minimum psi value to be considered (Default: 1.0)
#' @param psi_max_initial Maximum psi value to be considered (Default: 5.4)
#' @param rho_min_initial Minimum rho value to be considered (Default: 0.1)
#' @param rho_max_initial Maximum rho value to be considered (Default: 1.05)
#' @param chr_names A vector with chromosome names used for plotting
#' @return A list with fields output_optimum_pair, output_optimum_pair_without_ref, distance, distance_without_ref, minimise and is_ref_better
#' @export
run_clonal_ASCAT <- function(
  lrr, baf, lrrsegmented,
  bafsegmented, chromosomes,
  segBAF_table, input_optimum_pair,
  dist_choice, distancepng = NA,
  copynumberprofilespng = NA,
  nonroundedprofilepng = NA,
  gamma_param, read_depth,
  uninformative_baf_threshold,
  allow100percent,
  reliabilityFile = NA,
  psi_min_initial = 1.0,
  psi_max_initial = 5.4,
  rho_min_initial = 0.1,
  rho_max_initial = 1.05,
  chr_names
) {
  siglevel_BAF <- 0.05
  maxdist_BAF <- 0.01

  # DCW 160314 - much more lenient logR thresholds (allow anything!)
  #  # TODO: This parameter is pushed down to is_segment_clonal but not used there (maybe not used at all?)
  siglevel_LogR <- -0.01
  maxdist_LogR <- 1

  ininitial_bounds <- list(psi_min = psi_min_initial, psi_max = psi_max_initial, rho_min = rho_min_initial, rho_max = rho_max_initial)

  new_bounds <- get_new_bounds(input_optimum_pair, ininitial_bounds)


  ch <- chromosomes
  b <- bafsegmented
  r <- lrrsegmented[names(bafsegmented)]

  s <- get_segment_info(lrrsegmented, segBAF_table)
  # Make sure no segment of length 1 remains - TODO: this should not occur and needs to be prevented upstream
  s <- s[s[, 3] > 1, ]
  dist_matrix_info <- create_distance_matrix_clonal(s, dist_choice, gamma_param, read_depth, siglevel_BAF, maxdist_BAF, siglevel_LogR, maxdist_LogR, uninformative_baf_threshold, new_bounds) # kjd 10-2-2013

  d <- dist_matrix_info$distance_matrix
  minimise <- dist_matrix_info$minimise

  # DCW 210314
  if (minimise) {
    best.distance <- min(d)
  } else {
    best.distance <- max(d)
  }

  ref_seg_matrix <- dist_matrix_info$ref_seg_matrix

  ref_major <- dist_matrix_info$ref_major
  ref_minor <- dist_matrix_info$ref_minor

  #########################################################

  ret <- find_centroid_of_global_minima(
    d, ref_seg_matrix, ref_major,
    ref_minor, s, dist_choice, minimise,
    new_bounds, distancepng, gamma_param,
    siglevel_BAF, maxdist_BAF, siglevel_LogR,
    maxdist_LogR, allow100percent,
    uninformative_baf_threshold, read_depth
  )
  optima_info_without_ref <- ret$optima_info_without_ref
  optima_info <- ret$optima_info

  nropt <- optima_info$nropt
  psi_opt1 <- optima_info$psi_opt1
  rho_opt1 <- optima_info$rho_opt1
  ploidy_opt1 <- optima_info$ploidy_opt1
  goodness_of_fit_opt1 <- optima_info$goodness_of_fit_opt1

  distance.from.ref.seg <- goodness_of_fit_opt1

  is_ref_better <- FALSE
  if (is.na(rho_opt1)) {
    print("reference segment did not provide a possible solution")
  } else if (psi_opt1 >= psi_min_initial && psi_opt1 <= psi_max_initial && rho_opt1 >= rho_min_initial && rho_opt1 <= rho_max_initial && ((minimise && distance.from.ref.seg < best.distance) || (!minimise && distance.from.ref.seg > best.distance))) {
    is_ref_better <- T
    print("reference segment gives better results than grid search")
  } else {
    print("reference segment gives no better results than grid search. Reverting to grid search solution")
  }

  psi_without_ref <- optima_info_without_ref$psi_opt1
  rho_without_ref <- optima_info_without_ref$rho_opt1
  ploidy_without_ref <- optima_info_without_ref$ploidy_opt1
  goodness_of_fit_without_ref <- optima_info_without_ref$goodness_of_fit_opt1

  #########################################################

  if (nropt > 0) {
    rho <- rho_without_ref
    psi <- psi_without_ref
    ploidy <- ploidy_without_ref
    goodness_of_fit <- goodness_of_fit_without_ref * 100
    nAfull <- (rho - 1 - (b - 1) * 2^(r / gamma_param) * ((1 - rho) * 2 + rho * psi)) / rho
    nBfull <- (rho - 1 + b * 2^(r / gamma_param) * ((1 - rho) * 2 + rho * psi)) / rho
    nA <- pmax(round(nAfull), 0)
    nB <- pmax(round(nBfull), 0)

    rBacktransform <- gamma_param * log((rho * (nA + nB) + (1 - rho) * 2) / ((1 - rho) * 2 + rho * psi), 2)
    bBacktransform <- (1 - rho + rho * nB) / (2 - 2 * rho + rho * (nA + nB))
    rConf <- ifelse(abs(rBacktransform) > 0.15, pmin(100, pmax(0, 100 * (1 - abs(rBacktransform - r) / abs(r)))), NA)
    bConf <- ifelse(bBacktransform != 0.5, pmin(100, pmax(0, ifelse(b == 0.5, 100, 100 * (1 - abs(bBacktransform - b) / abs(b - 0.5))))), NA)
    # DCW 150711 - get deviations from expected values
    if (!is.na(reliabilityFile)) {
      data.table::fwrite(data.frame(segmentedBAF = b, backTransformedBAF = bBacktransform, confidenceBAF = bConf, segmentedR = r, backTransformedR = rBacktransform, confidenceR = rConf, nA = nA, nB = nB, nAfull = nAfull, nBfull = nBfull), reliabilityFile, sep = ",", row.names = FALSE)
    }

    # Make plots
    if (!is.na(copynumberprofilespng)) {
      grDevices::png(
        filename = copynumberprofilespng,
        width = 2000, height = 500,
        res = 200, type = "cairo"
      )
    }
    ASCAT::ascat.plotAscatProfile(
      n1all = nA, n2all = nB,
      heteroprobes = TRUE,
      ploidy = ploidy, rho = rho,
      goodness_of_fit = goodness_of_fit, nonaberrant = FALSE,
      ch = ch, lrr = lrr,
      bafsegmented = bafsegmented,
      chrs = chr_names
    )
    if (!is.na(copynumberprofilespng)) {
      grDevices::dev.off()
    }

    # separated plotting from logic: create nonrounded copy number profile plot here
    if (!is.na(nonroundedprofilepng)) {
      grDevices::png(
        filename = nonroundedprofilepng,
        width = 2000, height = 500,
        res = 200, type = "cairo"
      )
    }
    ASCAT::ascat.plotNonRounded(
      ploidy = ploidy, rho = rho,
      goodness_of_fit = goodness_of_fit,
      nonaberrant = FALSE, nAfull = nAfull,
      nBfull = nBfull, bafsegmented = bafsegmented,
      ch = ch, lrr = lrr, chrs = chr_names
    )
    if (!is.na(nonroundedprofilepng)) {
      grDevices::dev.off()
    }
  }

  # Recalculate the psi_t for this rho using only clonal segments
  psi_t <- recalc_psi_t(psi_without_ref, rho_without_ref, gamma_param, lrrsegmented, segBAF_table, siglevel_BAF, maxdist_BAF, include_subcl_segments = FALSE)

  # If there aren't any clonally fit segments, the above yields NA. In this case, revert to the original grid search psi_t
  if (is.na(psi_t)) {
    print("Recalculated psi_t was NA, reverting to grid search solution. This occurs when no segment could be fit with a clonal state, check sample for contamination")
    psi_t <- psi_without_ref
  }

  output_optimum_pair <- list(psi = psi_opt1, rho = rho_opt1, ploidy = ploidy_opt1)
  # output_optimum_pair_without_ref = list(psi = psi_without_ref, rho = rho_without_ref, ploidy = ploidy_without_ref)
  # Use the recalculated psi_t from the clonal segments as our final estimate of psi_t which is data driven with rho fixed
  output_optimum_pair_without_ref <- list(psi = psi_t, rho = rho_without_ref, ploidy = ploidy_without_ref)
  return(list(output_optimum_pair = output_optimum_pair, output_optimum_pair_without_ref = output_optimum_pair_without_ref, distance = distance.from.ref.seg, distance_without_ref = best.distance, minimise = minimise, is_ref_better = is_ref_better)) # kjd 20-2-2014, adapted by DCW 140314
}

#' Function extends the ASCAT \code{make_segments} function to make segments
#' of constant BAF and LogR. This function returns a matrix with for each
#' segment the LogR, BAF, the length of the segment (twice), and the mean and
#' standard deviation of the BAF values
#' @noRd
get_segment_info <- function(segLogR, segBAF_table) {
  # Column 5: Segmented BAF (b), Column 4: Phased BAF (BAFke)
  b_raw <- segBAF_table[, 5]
  b_phased <- segBAF_table[, 4]

  # Match original make_segments(r, b) call
  pcf_segments <- make_segments(segLogR, b_raw)

  # To match 'which(segBAF_table[, 5] == BAF_req)' exactly:
  # We group by the BAF value itself, not the segment position.
  # collapse::GRP is extremely fast for this.
  val_g <- collapse::GRP(b_raw)

  # Calculate stats for every unique BAF value once (O(N))
  all_means <- as.numeric(collapse::fmean(b_phased, val_g))
  all_sds <- as.numeric(collapse::fsd(b_phased, val_g))
  all_sizes <- as.numeric(collapse::fnobs(b_phased, val_g))

  # Map the calculated stats to each segment by matching the segment's BAF
  # value back to the group values.
  match_idx <- match(pcf_segments[, "b"], val_g$groups)

  # Build final matrix
  segs <- cbind(
    pcf_segments,
    size = all_sizes[match_idx],
    mean = all_means[match_idx],
    sd   = all_sds[match_idx]
  )

  return(segs)
}


#' Optimized Segment Maker
make_segments <- function(r, b) {
  keep <- !is.na(r) & !is.na(b)
  if (!any(keep)) {
    return(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c("r", "b", "length"))))
  }

  r_clean <- r[keep]
  b_clean <- b[keep]
  ids <- data.table::rleid(r_clean, b_clean)

  # To get 'r' and 'b' for each segment (the values at the start of each group):
  # which(!duplicated(ids)) finds the index of the first row of every new segment.
  first_idx <- which(!duplicated(ids))

  # To get 'length' (the count of rows in each group):
  # collapse::fnobs counts observations per group ID extremely quickly.
  # we cast to numeric to match the original matrix type perfectly.
  res_len <- as.numeric(collapse::fnobs(r_clean, g = ids))

  # Creating the matrix via cbind on atomic vectors is nearly instantaneous.
  # This avoids the 'as.matrix' call that slows down data.frame-based approaches.
  pcf_segments <- cbind(
    r      = r_clean[first_idx],
    b      = b_clean[first_idx],
    length = res_len
  )

  return(pcf_segments)
}
