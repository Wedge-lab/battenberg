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
    log_info("reference segment did not provide a possible solution")
  } else if (psi_opt1 >= psi_min_initial && psi_opt1 <= psi_max_initial && rho_opt1 >= rho_min_initial && rho_opt1 <= rho_max_initial && ((minimise && distance.from.ref.seg < best.distance) || (!minimise && distance.from.ref.seg > best.distance))) {
    is_ref_better <- T
    log_info("reference segment gives better results than grid search")
  } else {
    log_info("reference segment gives no better results than grid search. Reverting to grid search solution")
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
    log_info("Recalculated psi_t was NA, reverting to grid search solution. This occurs when no segment could be fit with a clonal state, check sample for contamination")
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
  # Fast removal of NAs
  keep <- which(!is.na(r) & !is.na(b))

  if (length(keep) == 0) {
    return(matrix(
      nrow = 0, ncol = 6,
      dimnames = list(NULL, c("r", "b", "length", "size", "mean", "sd"))
    ))
  }

  r_clean <- r[keep]
  b_clean <- b[keep]

  # 1. Robust Grouping
  # We round to 8 decimal places to avoid floating point noise breaking segments
  ids <- data.table::rleid(round(r_clean, 8), round(b_clean, 8))

  # 2. Ultra-fast Aggregation using collapse
  # We use ffirst to get the segment values and fnobs/fmean/fsd for the stats
  # g = ids tells collapse to perform these operations by group in C

  # pre-allocate matrix for speed
  n_seg <- ids[length(ids)]
  pcf_segments <- matrix(nrow = n_seg, ncol = 6)
  colnames(pcf_segments) <- c("r", "b", "length", "size", "mean", "sd")

  # Populate columns
  pcf_segments[, "r"] <- collapse::ffirst(r_clean, g = ids)
  pcf_segments[, "b"] <- collapse::ffirst(b_clean, g = ids)
  pcf_segments[, "length"] <- as.numeric(collapse::fnobs(r_clean, g = ids))
  pcf_segments[, "size"] <- pcf_segments[, "length"]
  pcf_segments[, "mean"] <- as.numeric(collapse::fmean(b_clean, g = ids))

  # Standard deviation requires a safety check for single-probe segments
  sds <- collapse::fsd(b_clean, g = ids)
  pcf_segments[, "sd"] <- ifelse(is.na(sds), 0, as.numeric(sds))

  return(pcf_segments)
}
