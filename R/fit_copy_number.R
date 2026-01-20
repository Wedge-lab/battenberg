#' Fit copy number
#'
#' Function that will fit a clonal copy number profile to segmented data. It
#' first matches the raw LogR with the segmented BAF to create segmented LogR.
#' Then ASCAT is run to obtain a clonal copy number profile. Beyond logRsegmented
#' it produces the rho_and_psi file and the cellularity_ploidy file.
#' @param samplename Samplename used to name the segmented logr output file
#' @param outputfile_prefix Prefix used for all output file names, except
#' logRsegmented
#' @param inputfile_baf_segmented Filename that points to the BAF segmented data
#' @param inputfile_baf Filename that points to the raw BAF data
#' @param inputfile_logr Filename that points to the raw LogR data
#' @param dist_choice The distance metric that is used internally to rank clonal
#' copy number solutions
#' @param ascat_dist_choice The distance metric used to obtain an initial
#' cellularity and ploidy estimate
#' @param min_ploidy The minimum ploidy to consider (Default 1.6)
#' @param max_ploidy The maximum ploidy to consider (Default 4.8)
#' @param min_rho The minimum cellularity to consider (Default 0.1)
#' @param max_rho The maximum cellularity to consider (Default 1.0)
#' @param min_goodness The minimum goodness of fit for a solution to have to be
#' considered (Default 63)
#' @param uninformative_baf_threshold The threshold beyond which BAF becomes
#' uninformative (Default 0.51)
#' @param gamma_param Technology parameter, compaction of Log R profiles.
#' Expected decrease in case of deletion in diploid sample, 100 "\%" aberrant
#' cells; 1 in ideal case, 0.55 of Illumina 109K arrays (Default 1)
#' @param use_preset_rho_psi Boolean whether to use user specified rho and psi
#' values (Default FALSE)
#' @param preset_rho A user specified rho to fit a copy number profile to
#' (Default NA)
#' @param preset_psi A user specified psi to fit a copy number profile to
#' (Default NA)
#' @param read_depth Legacy parameter that is no longer used (Default 30)
#' @param analysis A String representing the type of analysis to be run, this
#' determines whether the distance figure is produced (Default paired)
#' @param nthreads The number of paralel processes to run
#' @param enhanced_grid_search Flag to determine if the grid search should be performed with a higher number of steps (Default: FALSE)
#' @author dw9, sd11
#' @export
fit_copy_number <- function(
  samplename,
  outputfile_prefix,
  inputfile_baf_segmented,
  inputfile_baf,
  inputfile_logr,
  dist_choice,
  ascat_dist_choice,
  min_ploidy = 1.6,
  max_ploidy = 4.8,
  min_rho = 0.1,
  max_rho = 1.0,
  min_goodness = 63,
  uninformative_baf_threshold = 0.51,
  gamma_param = 1,
  use_preset_rho_psi = FALSE,
  preset_rho = NA,
  preset_psi = NA,
  read_depth = 30,
  analysis = "paired",
  nthreads = 1,
  enhanced_grid_search = FALSE
) {
  assert_file_exists(inputfile_baf_segmented)
  assert_file_exists(inputfile_baf)
  assert_file_exists(inputfile_logr)

  if ((max_ploidy - min_ploidy) < 0.05) {
    log_failure("Supplied ploidy range must be larger than 0.05: \\
                {min_ploidy}-{max_ploidy}")
  }

  # Read in the required data
  segmented.BAF.data <- read_bafsegmented(inputfile_baf_segmented)

  data.table::setDF(segmented.BAF.data)

  raw.BAF.data <- read_baf_as_data_frame(inputfile_baf)
  names(raw.BAF.data)[3] <- samplename

  raw.logR.data <- read_baf_as_data_frame(inputfile_logr)
  names(raw.logR.data)[3] <- samplename

  # Remove duplicates and set rownames
  identifiers <- paste(segmented.BAF.data[, 1], segmented.BAF.data[, 2], sep = "_")
  dups <- which(duplicated(identifiers))
  if (length(dups) > 0) {
    segmented.BAF.data <- segmented.BAF.data[-dups, ]
    identifiers <- identifiers[-dups]
  }
  rownames(segmented.BAF.data) <- identifiers

  # Drop NAs
  raw.BAF.data <- raw.BAF.data[!is.na(raw.BAF.data[, 3]), ]
  raw.logR.data <- raw.logR.data[!is.na(raw.logR.data[, 3]), ]

  BAF.data <- list()
  logR.data <- list()
  segmented.logR.data <- list()
  matched.segmented.BAF.data <- list()

  gsubchr <- function(chr) gsub("chr", "", as.character(chr))
  chr_names <- gsubchr(unique(segmented.BAF.data[, 1]))

  segmented.BAF.data$Chromosome <- gsubchr(segmented.BAF.data$Chromosome)
  raw.BAF.data$Chromosome <- gsubchr(raw.BAF.data$Chromosome)
  raw.logR.data$Chromosome <- gsubchr(raw.logR.data$Chromosome)

  baf_segmented_split <- split(segmented.BAF.data, f = segmented.BAF.data$Chromosome)
  baf_split <- split(raw.BAF.data, f = raw.BAF.data$Chromosome)
  logr_split <- split(raw.logR.data, f = raw.logR.data$Chromosome)

  # For each chromosome: Merge and initial alignment
  for (chr in chr_names) {
    chr.BAF.data <- baf_split[[chr]]
    chr.segmented.BAF.data <- baf_segmented_split[[chr]]

    if (is.null(chr.BAF.data) || nrow(chr.BAF.data) == 0) next

    merged <- merge(chr.segmented.BAF.data, chr.BAF.data, by = "Position", all = TRUE)

    matched.segmented.BAF.data[[chr]] <- merged
    BAF.data[[chr]] <- merged[, c("Position", samplename), drop = FALSE]

    chr.logR.data <- logr_split[[chr]]
    if (!is.null(chr.logR.data) && nrow(chr.logR.data) > 0) {
      merged_logR <- merge(merged, chr.logR.data, by = "Position", all = TRUE)
      logR.data[[chr]] <- merged_logR[, c(1, ncol(merged_logR)), drop = FALSE]
      segmented.logR.data[[chr]] <- merged_logR[, c(1, 3), drop = FALSE]
    }
  }

  # Sync the dataframes: Ensure absolute row-parity across all lists
  for (chrom in chr_names) {
    if (is.null(matched.segmented.BAF.data[[chrom]]) || is.null(logR.data[[chrom]])) {
      matched.segmented.BAF.data[[chrom]] <- logR.data[[chrom]] <- BAF.data[[chrom]] <- segmented.logR.data[[chrom]] <- NULL
      next
    }

    # Match based on the common Position column
    selection <- matched.segmented.BAF.data[[chrom]]$Position %in% logR.data[[chrom]]$Position

    if (sum(selection) == 0) {
      matched.segmented.BAF.data[[chrom]] <- logR.data[[chrom]] <- BAF.data[[chrom]] <- segmented.logR.data[[chrom]] <- NULL
      next
    }

    # Subset everything using the same selection vector
    matched.segmented.BAF.data[[chrom]] <- matched.segmented.BAF.data[[chrom]][selection, ]
    segmented.logR.data[[chrom]] <- segmented.logR.data[[chrom]][selection, ]
    BAF.data[[chrom]] <- BAF.data[[chrom]][selection, ]

    # Final alignment of the raw LogR list
    logR.data[[chrom]] <- logR.data[[chrom]][logR.data[[chrom]]$Position %in% matched.segmented.BAF.data[[chrom]]$Position, ]
  }

  log_info("Combining split data frames into final structures...")
  # Combine split data frames
  matched.segmented.BAF.data <- data.table::rbindlist(matched.segmented.BAF.data)
  segmented.logR.data <- data.table::rbindlist(segmented.logR.data)
  BAF.data <- data.table::rbindlist(BAF.data)
  logR.data <- data.table::rbindlist(logR.data)

  log_info("Final data synchronization check: {nrow(matched.segmented.BAF.data)} \\
           loci remaining.")
  # Fail Fast: Verify synchronization
  if (nrow(matched.segmented.BAF.data) < 100) {
    log_failure("Too few SNPs ({nrow(matched.segmented.BAF.data)}) remain after synchronization. Data is likely unusable.")
  }
  stopifnot(nrow(matched.segmented.BAF.data) == nrow(logR.data))

  # Prepare vectors for ASCAT
  # We use [[2]] to grab the value column (since [[1]] is Position)
  segBAF <- 1 - matched.segmented.BAF.data[[5]]
  segLogR <- segmented.logR.data[[2]]
  logR <- logR.data[[2]]

  # Crucial: Use rownames to allow ASCAT to map segments to probes
  row_ids <- paste(matched.segmented.BAF.data$Chromosome, matched.segmented.BAF.data$Position, sep = "_")
  names(segBAF) <- row_ids
  names(segLogR) <- row_ids
  names(logR) <- row_ids

  # Calculate chromosome indices for the combined vectors
  chr_segs <- list()
  for (i in seq_along(chr_names)) {
    chr_segs[[i]] <- which(matched.segmented.BAF.data$Chromosome == chr_names[i])
  }

  # Run ASCAT Grid Search
  if (use_preset_rho_psi) {
    log_info("Using preset rho ({preset_rho}) and psi ({preset_psi}). Skipping grid search.")
    ascat_optimum_pair <- list(rho = preset_rho, psi = preset_psi, ploidy = preset_psi)
  } else {
    log_info("Starting ASCAT Grid Search (this may take several minutes)...")
    distance_outfile <- paste0(outputfile_prefix, "distance.png")
    copynumberprofile_outfile <- paste0(
      outputfile_prefix,
      "copynumberprofile.png"
    )
    nonroundedprofile_outfile <- paste0(
      outputfile_prefix,
      "nonroundedprofile.png"
    )
    cnaStatusFile <- paste0(
      outputfile_prefix,
      "copynumber_solution_status.txt"
    )

    if (enhanced_grid_search) {
      log_info("Running ENHANCED grid search...")
      ascat_optimum_pair <- runASCAT_enhanced(
        logR, 1 - BAF.data[[2]], segLogR, segBAF,
        chr_segs, ascat_dist_choice, distance_outfile,
        copynumberprofile_outfile, nonroundedprofile_outfile,
        cnaStatusFile = cnaStatusFile, gamma = gamma_param,
        allow100percent = TRUE, min_ploidy = min_ploidy,
        max_ploidy = max_ploidy, min_rho = min_rho, max_rho = max_rho,
        min_goodness = min_goodness, chr_names = chr_names,
        analysis = analysis,
        uninformative_baf_threshold = uninformative_baf_threshold,
        nthreads = nthreads
      )
    } else {
      log_info("Running STANDARD grid search...")
      ascat_optimum_pair <- runASCAT(
        logR, 1 - BAF.data[[2]], segLogR, segBAF,
        chr_segs, ascat_dist_choice,
        distancepng = distance_outfile,
        copynumberprofilespng = copynumberprofile_outfile,
        nonroundedprofilepng = nonroundedprofile_outfile,
        cnaStatusFile = cnaStatusFile,
        gamma = gamma_param, allow100percent = TRUE,
        min_ploidy = min_ploidy, max_ploidy = max_ploidy,
        min_rho = min_rho, max_rho = max_rho,
        min_goodness = min_goodness, chr_names = chr_names, analysis = analysis,
        uninformative_baf_threshold = uninformative_baf_threshold,
        nthreads = nthreads
      )
    }
    log_info("Grid Search complete. Optimum found: \\
             Rho={ascat_optimum_pair$rho}, Psi={ascat_optimum_pair$psi}")

    # guard rail - check for valid solution
    if (is.na(ascat_optimum_pair$rho) || is.na(ascat_optimum_pair$psi)) {
      log_failure("Grid search failed to find a valid purity/ploidy solution. Data might be too noisy.")
    }
  }

  log_info("Running final clonal ASCAT model fit...")
  # Final clonal ASCAT run
  out <- run_clonal_ASCAT(
    logR, 1 - BAF.data[[2]], segLogR, segBAF, chr_segs,
    matched.segmented.BAF.data, ascat_optimum_pair, dist_choice,
    paste0(outputfile_prefix, "second_distance.png"),
    paste0(outputfile_prefix, "second_copynumberprofile.png"),
    paste0(outputfile_prefix, "second_nonroundedprofile.png"),
    gamma_param = gamma_param, read_depth, uninformative_baf_threshold,
    allow100percent = TRUE, psi_min_initial = min_ploidy,
    psi_max_initial = max_ploidy, rho_min_initial = min_rho,
    rho_max_initial = max_rho, chr_names = chr_names,
    nthreads = nthreads
  )

  if (is.na(out$output_optimum_pair$rho) || is.na(out$output_optimum_pair$psi)) {
    log_failure("Final clonal model fit failed to identify a valid purity/ploidy solution.")
  }
  d <- out$dist_matrix_info$distance_matrix
  if (all(is.na(d)) || all(is.infinite(d))) {
    log_failure("Distance matrix is entirely NA or Inf. No valid copy number solution possible.")
  }
  log_info("ASCAT modeling complete for {samplename}. Writing output files.")
  # Save results
  rho_psi_output <- data.frame(
    rho = c(ascat_optimum_pair$rho, out$output_optimum_pair_without_ref$rho, out$output_optimum_pair$rho),
    psi = c(ascat_optimum_pair$psi, out$output_optimum_pair_without_ref$psi, out$output_optimum_pair$psi),
    ploidy = c(ascat_optimum_pair$ploidy, out$output_optimum_pair_without_ref$ploidy, out$output_optimum_pair$ploidy),
    distance = c(NA, out$distance_without_ref, out$distance),
    is_best = c(NA, !out$is_ref_better, out$is_ref_better),
    row.names = c("ASCAT", "FRAC_GENOME", "REF_SEG")
  )
  data.table::fwrite(rho_psi_output,
    paste0(outputfile_prefix, "rho_and_psi.txt"),
    sep = "\t"
  )
}

#' Fit subclonal copy number
#'
#' This function fits a subclonal copy number profile where a clonal profile is unlikely.
#' It goes over each segment of a clonal copy number profile and does a simple t-test. If the
#' test is significant it is unlikely that the data can be explained by a single copy number
#' state. We therefore fit a second state, i.e. there are two cellular populations with each
#' a different state: Subclonal copy number.
#' @param sample_name Name of the sample, used in figures
#' @param baf_segmented_file String that points to a file with segmented BAF output
#' @param logr_file String that points to the raw LogR file to be used in the
#' subclonal copy number figures
#' @param rho_psi_file String pointing to the rho_and_psi file generated by
#' \code{fit_copy_number}
#' @param output_file Filename of the file where the final copy number fit will be
#' written to
#' @param output_figures_prefix Prefix of the filenames for the chromosome specific
#' copy number figures
#' @param output_gw_figures_prefix Prefix of the filenames for the genome wide copy
#' number figures
#' @param chr_names Vector of allowed chromosome names
#' @param masking_output_file Filename of where the masking details need to be
#' written. Masking is performed to remove very high copy number state segments
#' @param max_allowed_state The maximum CN state allowed (Default 250)
#' @param cn_upper_limit The maximum CN that can be called (Default 1000)
#' @param prior_breakpoints_file A two column file with prior breakpoints, possibly
#' from structural variants. This file must contain two columns: chromosome and
#' position. These are used when making the figures
#' @param gamma Technology specific scaling parameter for LogR (Default 1)
#' @param segmentation_gamma Legacy parameter that is no longer used (Default NA)
#' @param siglevel Threshold under which a p-value becomes significant. When it is
#' significant a second copy number state will be fitted (Default 0.05)
#' @param maxdist Slack in BAF space to allow a segment to be off it's optimum
#' before becoming significant. A segment becomes significant very quickly when a
#' breakpoint is missed, this parameter alleviates the effect (Default 0.01)
#' @param noperms The number of permutations to be run when bootstrapping the
#' confidence intervals on the copy number state of each segment (Default 1000)
#' @param seed Seed to set when performing bootstrapping (Default: Current time)
#' @param calc_seg_baf_option Various options to recalculate the BAF of a segment.
#' Options are: 1 - median, 2 - mean, 3 - ifelse median==0|1, mean, median.
#' (Default: 3)
#' @param verbose_logging Print out more information during the run (Default: FALSE)
#' @param nthreads The number of paralel processes to run
#' @author dw9, sd11
#' @export
call_subclones <- function(
  sample_name, baf_segmented_file,
  logr_file, rho_psi_file, output_file,
  output_figures_prefix, output_gw_figures_prefix,
  chr_names, masking_output_file,
  max_allowed_state = 250, cn_upper_limit = 1000,
  prior_breakpoints_file = NULL, gamma = 1,
  segmentation_gamma = NA, siglevel = 0.05,
  maxdist = 0.01, noperms = 1000, seed = as.integer(Sys.time()),
  calc_seg_baf_option = 3, verbose_logging = FALSE,
  nthreads = 1
) {
  set.seed(seed)

  # Load and calculate initial rho/psi metrics
  res <- load_rho_psi_file(rho_psi_file)
  rho <- res$rho
  psit <- res$psit
  psi <- (rho * psit) + (2 * (1 - rho))
  goodness <- res$goodness

  # Load BAF data and handle possible row-name artifacts ("X")
  BAFvals <- read_bafsegmented(baf_segmented_file)
  if ("X" %in% colnames(BAFvals)) {
    BAFvals <- BAFvals[, -1, with = FALSE]
  }

  # Positional indexing for generalizability: Col 3 = BAF, Col 5 = BAFseg
  BAF <- BAFvals[, 3]
  BAFseg <- BAFvals[, 5]
  SNPpos <- BAFvals[, c(1, 2), drop = FALSE]

  # Load LogR data and handle row-name artifacts
  LogRvals <- read_logr(logr_file)
  if (identical(colnames(LogRvals)[1], "X")) {
    LogRvals <- LogRvals[, -1, drop = FALSE]
  }

  ctrans <- ctrans_logR <- stats::setNames(seq_along(chr_names), chr_names)

  # First Pass: Determine Copy Number and Merge Segments
  res_cn <- determine_copynumber(
    BAFvals, LogRvals, rho, psi, gamma,
    ctrans, ctrans_logR, maxdist, siglevel, noperms, cn_upper_limit
  )

  # Refine via merging
  merge_res <- merge_segments(
    res_cn$subcloneres, BAFvals, LogRvals,
    rho, psi, gamma, calc_seg_baf_option, TRUE
  )
  BAFvals <- merge_res$bafsegmented

  # Second Pass: Final Copy Number Determination
  res_final <- determine_copynumber(
    BAFvals, LogRvals, rho, psi, gamma,
    ctrans, ctrans_logR, maxdist, siglevel,
    noperms, cn_upper_limit
  )
  subcloneres <- res_final$subcloneres
  BAFpvals <- res_final$BAFpvals

  # Mask high CN artifacts
  mask_res <- mask_high_cn_segments(subcloneres, BAFvals, max_allowed_state)
  subcloneres <- mask_res$subclones

  data.table::fwrite(
    list(
      samplename = sample_name,
      masked_count = mask_res$masked_count,
      masked_size = mask_res$masked_size,
      max_allowed_state = max_allowed_state
    ),
    file = masking_output_file,
    quote = FALSE,
    sep = "\t",
  )

  # Generate output paths
  base_out <- tools::file_path_sans_ext(output_file)
  ext_out <- tools::file_ext(output_file)

  data.table::fwrite(
    subcloneres[, c(1:3, 8:13)], output_file,
    quote = FALSE, sep = "\t"
  )
  data.table::fwrite(
    subcloneres, paste0(base_out, "_extended.", ext_out),
    quote = FALSE, sep = "\t"
  )

  subcloneres$length <- subcloneres$endpos - subcloneres$startpos

  # Behavior: Identical, but using which() ensures integer indexing for safe exclusion
  diploid_idx <- which(subcloneres$nMaj1_A == 1 & subcloneres$nMin1_A == 1 & subcloneres$frac1_A == 1)

  # Behavior: Identical. The if-statement handles the integer(0) case safely
  cna <- if (length(diploid_idx) > 0) subcloneres[-diploid_idx, ] else subcloneres

  # Use fsubset for subclonal filtering
  # Behavior: Identical. fsubset handles 0-row matches more cleanly than base [,]
  subcloneres_subclonal <- collapse::fsubset(
    subcloneres, subcloneres$frac1_A < 1
  )

  # Pre-calculate sums using collapse::fsum
  cna_total_len <- collapse::fsum(cna$length, na.rm = TRUE)

  # Logic Gate
  # Behavior: Identical. Checks for 0 rows or 0 total length
  if (nrow(cna) == 0 || cna_total_len == 0 || nrow(subcloneres_subclonal) == 0) {
    goodness <- 1.0
  } else {
    subclonal_total_len <- collapse::fsum(subcloneres_subclonal$length, na.rm = TRUE)
    subclonal_fraction <- subclonal_total_len / cna_total_len

    goodness <- max(0, min(1, 1 - subclonal_fraction))
  }

  log_info("PGA.is.clonal = {sprintf('%2.1f%%', goodness * 100)}")

  # Visualization
  segment_breakpoints <- collapse_bafsegmented_to_segments(BAFvals)
  has_prior <- !is.null(prior_breakpoints_file) &&
    !is.na(prior_breakpoints_file) &&
    prior_breakpoints_file != "NA"

  if (has_prior) {
    svs <- data.table::fread(prior_breakpoints_file, data.table = FALSE)
  }

  parallel::mclapply(chr_names, function(chr) {
    chr_idx <- SNPpos[, 1] == chr
    pos <- SNPpos[chr_idx, 2]

    if (length(pos) > 0) {
      svs_pos <- if (has_prior) {
        collapse::fsubset(
          svs, svs[[1]] == chr
        )[[2]] / 1e6
      } else {
        NULL
      }
      bp_chr <- collapse::fsubset(
        segment_breakpoints, segment_breakpoints[[1]] == chr
      )
      breakpoints_pos <- sort(unique(c(bp_chr[[2]], bp_chr[[3]]) / 1e6))

      grDevices::png(
        filename = paste0(output_figures_prefix, chr, ".png"),
        width = 2000, height = 2000, res = 200, type = "cairo"
      )
      create_subclonal_cn_plot(
        chrom = chr,
        chrom_position = pos / 1e6,
        LogRposke = LogRvals[LogRvals[, 1] == chr, 2],
        LogRchr = LogRvals[LogRvals[, 1] == chr, 3],
        BAFchr = BAF[chr_idx],
        BAFsegchr = BAFseg[chr_idx],
        BAFpvalschr = BAFpvals[chr_idx],
        subcloneres = subcloneres,
        siglevel = siglevel,
        x_min = min(pos) / 1e6,
        x_max = max(pos) / 1e6,
        title = paste(sample_name, ", chromosome ", chr),
        xlab = "Position (Mb)", ylab_logr = "LogR", ylab_baf = "BAF (phased)",
        breakpoints_pos = breakpoints_pos,
        svs_pos = svs_pos
      )
      grDevices::dev.off()
    }
  }, mc.cores = nthreads)

  # Clean up and calculate Ploidy
  subclones <- as.data.frame(subcloneres)
  seg_len <- floor((subcloneres$endpos - subcloneres$startpos) / 1000)

  # Calculate weighted states for min/maj
  calc_state <- function(n1, n2, f1, f2) {
    is_sub <- abs(n1 - n2) > 0
    is_sub[is.na(is_sub)] <- FALSE
    ifelse(is_sub, (n1 * f1) + (n2 * f2), n1)
  }

  state_min <- calc_state(subclones$nMin1_A, subclones$nMin2_A, subclones$frac1_A, subclones$frac2_A)
  state_maj <- calc_state(subclones$nMaj1_A, subclones$nMaj2_A, subclones$frac1_A, subclones$frac2_A)

  ploidy <- sum((state_min + state_maj) * seg_len, na.rm = TRUE) / sum(seg_len, na.rm = TRUE)

  # Final Outputs
  plot_gw_subclonal_cn(subclones, BAFvals, rho, ploidy, goodness, output_gw_figures_prefix, chr_names, sample_name)

  cp_out <- data.frame(purity = rho, ploidy = ploidy, psi = psit)
  data.table::fwrite(cp_out, paste0(sample_name, "_purity_ploidy.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
}

#' Given all the determined values make a copy number call for each segment
#'
#' @param BAFvals BAFsegmented data.frame with 5 columns
#' @param LogRvals Raw logR values in data.frame with 3 columns
#' @param rho Optimal rho value, the choosen cellularity
#' @param psi Optimal psi value, the choosen ploidy
#' @param gamma Platform gamma parameter
#' @param ctrans Named vector of chromosome names
#' @param ctrans_logR Named vector of chromosome names
#' @param maxdist Max distance a segment is tolerated to be not considered for subclonal copy number
#' @param siglevel Level at which a segment can become significantly different from the nearest clonal state
#' @param noperms Number of bootstrap permutations
#' @param cn_upper_limit Maximum number of CN that can be called
#' @return A data.frame with copy number determined for each segment
#' @author dw9
#' @noRd
determine_copynumber <- function(BAFvals, LogRvals, rho, psi, gamma, ctrans,
                                 ctrans.logR, maxdist, siglevel, noperms,
                                 cn_upper_limit) {
  # Standardizing inputs - stripped redundant as.vector calls
  BAFphased <- BAFvals[, 4]
  BAFseg <- BAFvals[, 5]
  BAFpos <- ctrans[BAFvals[, 1]] * 1e9 + BAFvals[, 2]
  LogRpos <- ctrans.logR[LogRvals[, 1]] * 1e9 + LogRvals[, 2]

  # Boundary logic
  switchpoints <- c(0, which(BAFseg[-1] != BAFseg[-length(BAFseg)] | BAFvals[-1, 1] != BAFvals[-nrow(BAFvals), 1]), length(BAFseg))
  BAFlevels <- BAFseg[switchpoints[-1]]

  res_list <- vector(mode = "list", length = length(BAFlevels))
  BAFpvals <- vector(length = length(BAFseg))

  # 1. Fast LogR averaging using collapse
  # Map each LogR probe to a segment index
  # LogRpos and segment boundaries (startpos/endpos) are both sorted globally
  # We can find which segment each LogR probe falls into.

  # Get all segment boundaries
  seg_starts <- BAFpos[switchpoints[-length(switchpoints)] + 1]
  seg_ends <- BAFpos[switchpoints[-1]]

  # findInterval returns index i such that seg_starts[i] <= LogRpos < seg_starts[i+1]
  # We need to ensure LogRpos <= seg_ends[i] as well (handling gaps)
  seg_ids <- findInterval(LogRpos, seg_starts)

  # Filter LogR probes that are within the matched segment's end and not infinite
  valid_logr <- seg_ids > 0 & LogRpos <= seg_ends[pmax(1, seg_ids)] & !is.infinite(LogRvals[[3]])

  # Calculate mean LogR per segment ID
  # We use collapse::fmean with the assigned group IDs
  seg_logr_means <- as.numeric(collapse::fmean(LogRvals[[3]][valid_logr], g = seg_ids[valid_logr]))

  # Map back to the BAFlevels (some segments might be missing LogR data)
  LogR_vec <- numeric(length(BAFlevels))
  LogR_vec[sort(unique(seg_ids[valid_logr]))] <- seg_logr_means

  # 2. Vectorized Clonal Math
  # BAFlevels (l) is normalized to be major allele freq (>= 0.5)
  l_vec <- pmax(BAFlevels, 1 - BAFlevels)

  # Precompute terms
  logr_factor <- 2^(LogR_vec / gamma)
  nMajor_vec <- (rho - 1 + l_vec * psi * logr_factor) / rho
  nMinor_vec <- (rho - 1 + (1 - l_vec) * psi * logr_factor) / rho

  # Handle physical impossibility (Negative nMinor)
  neg_minor <- nMinor_vec < 0 & !is.na(nMinor_vec)
  if (any(neg_minor)) {
    is_one <- l_vec == 1
    nMajor_vec[neg_minor & is_one] <- cn_upper_limit
    nMajor_vec[neg_minor & !is_one] <- nMajor_vec[neg_minor & !is_one] +
      l_vec[neg_minor & !is_one] * (0.01 - nMinor_vec[neg_minor & !is_one]) / (1 - l_vec[neg_minor & !is_one])
    nMinor_vec[neg_minor] <- 0.01
  }

  # 3. Vectorized is_segment_clonal-style testing
  # We need BAF_size, BAF_sd for each segment for the p-value
  # We can get these from the BAFphased data using the switchpoints
  baf_groups <- rep(seq_along(BAFlevels), diff(switchpoints))
  BAF_stats <- data.frame(
    mean = as.numeric(collapse::fmean(BAFphased, g = baf_groups)),
    sd   = as.numeric(collapse::fsd(BAFphased, g = baf_groups)),
    size = as.numeric(collapse::fnobs(BAFphased, g = baf_groups))
  )
  BAF_stats$sd[is.na(BAF_stats$sd)] <- 0

  # Call is_segment_clonal in one vectorized go
  # is_segment_clonal is already vectorized and returns best_nMaj, best_nMin, is_clonal
  # We need to ensure we have all required parameters
  best_clonal_res <- is_segment_clonal(
    LogR = LogR_vec,
    BAF_req = l_vec,
    BAF_length = BAF_stats$size, # approximating length with size
    BAF_size = BAF_stats$size,
    BAF_mean = BAF_stats$mean,
    BAF_sd = BAF_stats$sd,
    rho = rho,
    psi = psi,
    gamma_param = gamma,
    siglevel_BAF = siglevel,
    maxdist_BAF = maxdist
  )

  # Map p-values back to SNP-level BAFpvals
  # Note: is_segment_clonal (vectorized version) doesn't return pval currently,
  # but it sets is_clonal based on pval > siglevel.
  # We actually need the p-value ourselves to fill BAFpvals.
  # Let's extract that logic or re-calculate here.

  # Re-calculate best_level for p-value (Option 1 vs 2)
  # This matches the prioritized testing in determine_copynumber
  calc_baf_lev <- function(nM, nm) {
    num <- 1 - rho + rho * nM
    den <- 2 - 2 * rho + rho * (nM + nm)
    lev <- num / den
    lev[nM == 0 & nm == 0] <- 0.5
    lev
  }

  best_levels <- calc_baf_lev(best_clonal_res$nMaj, best_clonal_res$nMin)

  # Vectorized p-value calculation
  p_vals <- numeric(length(BAFlevels))
  valid_stats <- BAF_stats$size > 1 & BAF_stats$sd > 0
  if (any(valid_stats)) {
    p_vals[valid_stats] <- calc_Pvalue_t_twotailed(
      sample_size = BAF_stats$size[valid_stats],
      sample_mean = BAF_stats$mean[valid_stats],
      sample_SD   = BAF_stats$sd[valid_stats],
      mu_pop      = best_levels[valid_stats],
      max_dist    = maxdist
    )
  }

  # Fill BAFpvals (SNP level)
  BAFpvals <- p_vals[baf_groups]

  # 4. Process Subclonal Segments (Only for those where p_vals <= siglevel)
  # This part is harder to vectorize fully due to the bootstrap loop,
  # but we only do it for the subclonal subset.
  subclonal_idx <- which(p_vals <= siglevel)

  for (i in seq_along(BAFlevels)) {
    l <- l_vec[i]
    LogR <- LogR_vec[i]
    ntot <- nMajor_vec[i] + nMinor_vec[i]

    start_idx <- switchpoints[i] + 1
    end_idx <- switchpoints[i + 1]

    curr_start <- seg_starts[i] %% 1e9
    curr_end <- seg_ends[i] %% 1e9

    if (i %in% subclonal_idx) {
      # SUBCLONAL
      BAFke <- BAFphased[start_idx:end_idx]
      n_ke <- length(BAFke)
      sd_BAFke <- BAF_stats$sd[i]

      # Need all edges for subclonal optimization
      all_edges <- prioritizeCopyNumbers(
        rho = rho, psi = psi, BAF_req = l,
        nMajor = nMajor_vec[i], nMinor = nMinor_vec[i], full = TRUE
      )

      na_idx <- which(is.na(rowSums(all_edges)))
      if (length(na_idx) > 0) all_edges <- rbind(all_edges[-na_idx, ], all_edges[na_idx, ])

      nM1 <- all_edges[, 1]
      nmi1 <- all_edges[, 2]
      nM2 <- all_edges[, 3]
      nmi2 <- all_edges[, 4]

      # Vectorized math for tau across all 6 options
      tau <- (1 - rho + rho * nM2 - 2 * l * (1 - rho) - l * rho * (nmi2 + nM2)) /
        (l * rho * (nmi1 + nM1) - l * rho * (nmi2 + nM2) - rho * nM1 + rho * nM2)

      sdl <- sd_BAFke / sqrt(n_ke)

      # Optimized Delta method for sdtau
      calc_sdtau <- function(curr_l) {
        (1 - rho + rho * nM2 - 2 * curr_l * (1 - rho) - curr_l * rho * (nmi2 + nM2)) /
          (curr_l * rho * (nmi1 + nM1) - curr_l * rho * (nmi2 + nM2) - rho * nM1 + rho * nM2)
      }
      sdtau <- (abs(calc_sdtau(l + sdl) - tau) + abs(calc_sdtau(l - sdl) - tau)) / 2

      # Optimized Bootstrap (Vectorized)
      # We generate all samples at once
      boot_means <- colMeans(matrix(sample(BAFke, n_ke * noperms, replace = TRUE), nrow = n_ke))

      opt_data <- vector("list", 6)
      for (opt in seq_along(tau)) {
        pFrac <- (1 - rho + rho * nM2[opt] - 2 * boot_means * (1 - rho) - boot_means * rho * (nmi2[opt] + nM2[opt])) /
          (boot_means * rho * (nM1[opt] + nmi1[opt]) - boot_means * rho * (nM2[opt] + nmi2[opt]) - rho * nM1[opt] + rho * nM2[opt])

        o_frac <- sort(pFrac)
        opt_data[[opt]] <- c(
          nM1[opt], nmi1[opt], tau[opt], nM2[opt], nmi2[opt], 1 - tau[opt],
          sdtau[opt], collapse::fsd(pFrac), o_frac[round(0.025 * noperms)], o_frac[round(0.975 * noperms)]
        )
      }

      res_list[[i]] <- c(BAFvals$Chromosome[start_idx], curr_start, curr_end, l, p_vals[i], LogR, ntot, unlist(opt_data))
    } else {
      # CLONAL
      res_list[[i]] <- c(
        BAFvals$Chromosome[start_idx], curr_start, curr_end, l, p_vals[i], LogR, ntot,
        best_clonal_res$nMaj[i], best_clonal_res$nMin[i], 1, rep(NA, 57)
      )
    }
  }


  # Final formatting
  subcloneres <- as.data.frame(do.call(rbind, res_list))
  colnames(subcloneres) <- c("chr", "startpos", "endpos", "BAF", "pval", "LogR", "ntot", dynamic_names)


  # Modern fast type conversion
  subcloneres[-1] <- lapply(subcloneres[-1], function(x) as.numeric(as.character(x)))

  return(list(subcloneres = subcloneres, BAFpvals = BAFpvals))
}


#' Plot the copy number genome wide in two different ways. This creates the
#' Battenberg average profile where subclonal copy number is represented as a
#' mixture of two different states and the Battenberg subclones profile where
#' subclonal copy number is plotted as two different separate states. The thickness
#' of the line represents the fraction of tumour cells carying the particular state
#' @noRd
plot_gw_subclonal_cn <- function(subclones, BAFvals, rho, ploidy, goodness,
                                 output_gw_figures_prefix, chr_names,
                                 tumourname) {
  # Map start and end of each segment into the BAF values. The plot uses the index
  # of this BAF table as x-axis. Using O(M) vectorized approach.
  pos_min <- rep(NA_integer_, nrow(subclones))
  pos_max <- rep(NA_integer_, nrow(subclones))

  for (chr in unique(as.character(subclones$chr))) {
    baf_idx <- which(BAFvals$Chromosome == chr)
    if (length(baf_idx) == 0) next

    sub_idx <- which(subclones$chr == chr)
    curr_sub <- subclones[sub_idx, ]

    # Map each SNP to a segment index using findInterval
    # Original logic: startpos < Position <= endpos
    snp_to_seg <- findInterval(BAFvals$Position[baf_idx], curr_sub$startpos)

    # Validate SNPs are within the assigned segment's endpos
    valid_mask <- snp_to_seg > 0
    in_seg_mask <- valid_mask & BAFvals$Position[baf_idx] <= curr_sub$endpos[pmax(1, snp_to_seg)]

    if (any(in_seg_mask)) {
      seg_ids_found <- snp_to_seg[in_seg_mask]
      abs_snp_indices <- baf_idx[in_seg_mask]

      # Find min/max SNP index for each segment found
      pos_min[sub_idx[unique(seg_ids_found)]] <- collapse::fmin(abs_snp_indices, g = seg_ids_found)
      pos_max[sub_idx[unique(seg_ids_found)]] <- collapse::fmax(abs_snp_indices, g = seg_ids_found)
    }
  }

  # For those segments that are subclonal, we can now just subset the pre-calculated boundaries.
  is_subclonal <- which(subclones$frac1_A < 1)
  subcl_min <- pos_min[is_subclonal]
  subcl_max <- pos_max[is_subclonal]

  # Determine whether it's the major or the minor allele that is represented by two states
  is_subclonal_maj <- abs(subclones$nMaj1_A - subclones$nMaj2_A) > 0
  is_subclonal_min <- abs(subclones$nMin1_A - subclones$nMin2_A) > 0
  is_subclonal_maj[is.na(is_subclonal_maj)] <- FALSE
  is_subclonal_min[is.na(is_subclonal_min)] <- FALSE

  segment_states_min <- subclones$nMin1_A * ifelse(is_subclonal_min,
    subclones$frac1_A, 1
  ) +
    ifelse(is_subclonal_min, subclones$nMin2_A, 0) *
      ifelse(is_subclonal_min, subclones$frac2_A, 0)
  segment_states_maj <- subclones$nMaj1_A * ifelse(is_subclonal_maj,
    subclones$frac1_A, 1
  ) +
    ifelse(is_subclonal_maj, subclones$nMaj2_A, 0) *
      ifelse(is_subclonal_maj, subclones$frac2_A, 0)
  segment_states_tot <- segment_states_maj + segment_states_min

  # Determine which SNPs are on which chromosome, to be used as a proxy for chromosome size in the plots
  chr_segs <- lapply(seq_along(chr_names), function(ch) {
    which(BAFvals$Chromosome == chr_names[ch])
  })

  # Plot subclonal copy number as mixtures of two states
  grDevices::png(
    filename = paste(output_gw_figures_prefix, "_average.png", sep = ""),
    width = 2000, height = 500, res = 200, type = "cairo"
  )
  create_bb_plot_average(
    bafsegmented = BAFvals,
    ploidy = ploidy,
    rho = rho,
    goodness_of_fit = goodness,
    pos_min = pos_min,
    pos_max = pos_max,
    segment_states_min = segment_states_min,
    segment_states_tot = segment_states_tot,
    chr_segs = chr_segs,
    chr_names = chr_names,
    tumourname = tumourname
  )
  grDevices::dev.off()

  # Plot subclonal copy number as two separate states
  grDevices::png(
    filename = paste(output_gw_figures_prefix, "_subclones.png", sep = ""),
    width = 2000, height = 500, res = 200, type = "cairo"
  )
  create_bb_plot_subclones(
    bafsegmented = BAFvals,
    subclones = subclones,
    ploidy = ploidy,
    rho = rho,
    goodness_of_fit = goodness,
    pos_min = pos_min,
    pos_max = pos_max,
    subcl_min = subcl_min,
    subcl_max = subcl_max,
    is_subclonal = is_subclonal,
    is_subclonal_maj = is_subclonal_maj,
    is_subclonal_min = is_subclonal_min,
    chr_segs = chr_segs,
    chr_names = chr_names,
    tumourname = tumourname
  )
  grDevices::dev.off()
}

#' Collapse a BAFsegmented file into segment start and end points
#'
#' This function looks through the BAFsegmented for stretches of equal
#' BAFseg and records the start and end coordinates in a data.frame
#' @param bafsegmented The BAFsegmented output from segmentation
#' @return A data.frame with columns chromosome, start and end
#' @author sd11
#' @noRd
collapse_bafsegmented_to_segments <- function(bafsegmented) {
  stopifnot(all(c("Chromosome", "Position", "BAFseg") %in% colnames(bafsegmented)))

  segments_list <- bafsegmented |>
    split(~Chromosome) |>
    lapply(function(chrom_df) {
      chrom_df <- chrom_df[order(chrom_df$Position), ]
      rle_vals <- rle(chrom_df$BAFseg)
      cum_lengths <- cumsum(rle_vals$lengths)

      starts <- chrom_df$Position[c(1, cum_lengths[-length(cum_lengths)] + 1)]
      ends <- chrom_df$Position[cum_lengths]

      data.frame(
        chromosome = chrom_df$Chromosome[1],
        start = starts,
        end = ends,
        stringsAsFactors = FALSE
      )
    })

  # Combine and clean up without transform()
  segments <- do.call(rbind, segments_list)
  segments$chromosome <- as.character(segments$chromosome)
  rownames(segments) <- NULL

  return(segments)
}

#' Function to make additional figures
#'
#' @param samplename Name of the sample for the plot title
#' @param logr_file File containing all logR data
#' @param bafsegmented_file File containing the BAFsegmented data
#' @param logrsegmented_file File with the logRsegmented data
#' @param allelecounts_file Optional file with raw allele counts (Default: NULL)
#' @author sd11
#' @export
make_posthoc_plots <- function(samplename, logr_file, bafsegmented_file, logrsegmented_file, allelecounts_file = NULL) {
  # Make some post-hoc plots
  logr <- read_table_generic(logr_file)
  bafsegmented <- as.data.frame(read_table_generic(bafsegmented_file))
  logrsegmented <- as.data.frame(read_table_generic(logrsegmented_file, header = FALSE))
  colnames(logrsegmented) <- c("Chromosome", "Position", "logRseg")
  outputfile <- paste0(samplename, "_alleleratio.png")
  allele_ratio_plot(
    samplename = samplename, logr = logr,
    bafsegmented = bafsegmented, logrsegmented = logrsegmented,
    outputfile = outputfile, max.plot.cn = 8
  )

  if (!is.null(allelecounts_file)) {
    allelecounts <- as.data.frame(read_table_generic(allelecounts_file))
    outputfile <- paste0(samplename, "_coverage.png")
    coverage_plot(samplename, allelecounts, outputfile)
  }
}


#' Fit ChrX subclonal copy number (male only)
#'
#' Function to call ChrX copy number based on LogR (suitable for male samples).
#' Copy number cannot be called for the non-PAR region of ChrX due to the
#' hemizygosity of all 1000G SNPs. This function enables calling subclonal copy
#' number for the non-PAR region by segmenting LogR. A number of correction steps
#' are undertaken to account for the noisy nature of LogR. This function
#' requires the following libraries: copynumber, data.table and ggplot2. It reads
#' in three files generated by previous steps of Battenberg, namely
#' samplename_mutantLogR_gcCorrected.tab, samplename_purity_ploidy.txt
#' and samplename_copynumber_extended.txt.
#' This function will also update the Battenberg genome-wide profile plots
#' (average.png and subclones.png) to include the chrX profile by also
#' reading in the samplename.BAFsegmented.txt and samplename_rho_psi.txt files
#' @param tumourname The sample name used for Battenberg (i.e. the tumour BAM
#' file name without the .bam extension)
#' @param X_gamma The PCF gamma value for segmentation of 1000G SNP LogR values
#' (Default 1000)
#' @param X_kmin The min number of SNPs to support a segment in PCF of LogR values
#' (Default 100)
#' @param genomebuild The genome build used in running Battenberg (hg19 or hg38)
#' @param AR Should the segment carrying the androgen receptor (AR) locus to be
#' visually distinguished in average plot? (Default TRUE)
#' @param prior_breakpoints_file A two column text file with prior genome-wide
#' breakpoints, possibly from structural variants. This file must contain two
#' columns with headers "chr" and "pos" representing chromosome and position.
#' @param chrom_names A vector containing the names of chromosomes to be included
#' in the final genome-wide Battenberg copy number plot with chrX
#' @author naser.ansari-pour
#' @export
callChrXsubclones <- function(
  tumourname, X_gamma = 1000,
  X_kmin = 100, genomebuild,
  AR = TRUE, prior_breakpoints_file = NULL,
  chrom_names, data_type = "wgs"
) {
  log_info("Processing sample: {tumourname}")

  # Set genome-specific coordinates
  if (genomebuild == "hg19") {
    par_regions <- c(2699520, 155260560)
    x_centromere <- c(58632012, 61632012)
    ar_locus <- data.frame(startpos = 66763874, endpos = 66950461)
  } else if (genomebuild == "hg38") {
    par_regions <- c(2781479, 156030895)
    x_centromere <- c(58605580, 62412542)
    ar_locus <- data.frame(startpos = 67544021, endpos = 67730619)
  } else {
    log_failure("Genomebuild not supported for callChrXsubclones")
  }

  # Load LogR data
  suffix <- if (tolower(data_type) == "wgs") "_mutantLogR_gcCorrected.tab" else "_mutantLogR.tab"
  pcf_input_raw <- read_table_generic(paste0(tumourname, suffix)) |> as.data.frame()

  # Identify chromosome notation and filter for non-PAR X regions
  chr_x_name <- unique(pcf_input_raw$Chromosome[pcf_input_raw$Chromosome %in% c("X", "chrX")])[1]
  pcf_input <- pcf_input_raw[pcf_input_raw$Chromosome == chr_x_name &
    pcf_input_raw$Position > par_regions[1] &
    pcf_input_raw$Position < par_regions[2], ]
  colnames(pcf_input)[3] <- tumourname
  log_info("Number of chrX nonPAR SNPs = {nrow(pcf_input)}")

  # Segmentation with optional prior breakpoints
  if (!is.null(prior_breakpoints_file)) {
    sv_data <- data.table::fread(prior_breakpoints_file, data.table = FALSE)
    sv_x <- sv_data[sv_data$chr %in% c("X", "chrX"), ]

    if (nrow(sv_x) > 0) {
      # Filter breakpoints within the valid LogR range
      valid_sv_pos <- sv_x$pos[sv_x$pos > min(pcf_input$Position) & sv_x$pos < max(pcf_input$Position)]
      breaks <- sort(unique(c(min(pcf_input$Position), valid_sv_pos, max(pcf_input$Position))))

      pcf_results <- list()
      for (j in 1:(length(breaks) - 1)) {
        subset_input <- pcf_input[pcf_input$Position >= breaks[j] & pcf_input$Position < breaks[j + 1], ]
        if (nrow(subset_input) > 0) {
          pcf_results[[length(pcf_results) + 1]] <- copynumber::pcf(subset_input, gamma = X_gamma, kmin = X_kmin)
        }
      }
      pcf_df <- do.call(rbind, pcf_results)
    } else {
      pcf_df <- copynumber::pcf(pcf_input, gamma = X_gamma, kmin = X_kmin)
    }
  } else {
    pcf_df <- copynumber::pcf(pcf_input, gamma = X_gamma, kmin = X_kmin)
  }

  data.table::fwrite(pcf_df, paste0(tumourname, "_PCF_gamma_", X_gamma, "_chrX.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  # Load purity, ploidy and autosomal segments
  pupl <- data.table::fread(paste0(tumourname, "_purity_ploidy.txt"), data.table = FALSE)
  rho <- pupl[1, 1]
  psi_sample <- pupl$ploidy
  bb_data <- data.table::fread(paste0(tumourname, "_copynumber_extended.txt"), data.table = FALSE)

  # Calculate LogR correction based on autosomal diploid regions
  bb_dip <- bb_data[bb_data$nMaj1_A == 1 & bb_data$nMin1_A == 1 & bb_data$frac1_A == 1, ]
  bb_corr <- if (nrow(bb_dip) > 1) {
    -mean(bb_dip$LogR)
  } else {
    # WGD Fallback logic
    cnloh <- bb_data[bb_data$nMaj1_A == 2 & bb_data$nMin1_A == 0 & bb_data$frac1_A == 1, ]
    if (nrow(cnloh) > 0) -mean(cnloh$LogR) else -log2(2 / psi_sample)
  }

  # Estimate LogR Standard Deviation (Pixel Perfect SD logic)
  bb_g1 <- bb_data[bb_data$nMaj1_A == 2 & bb_data$nMin1_A == 1 & bb_data$frac1_A == 1, ]
  bb_g2 <- bb_data[bb_data$nMaj1_A == 3 & bb_data$nMin1_A == 1 & bb_data$frac1_A == 1, ]
  bb_g3 <- bb_data[bb_data$nMaj1_A == 4 & bb_data$nMin1_A == 1 & bb_data$frac1_A == 1, ]
  bb_sd_max <- max(c(
    collapse::fsd(bb_dip$LogR),
    collapse::fsd(bb_g1$LogR),
    collapse::fsd(bb_g2$LogR),
    collapse::fsd(bb_g3$LogR),
    0.05
  ), na.rm = TRUE)

  # Expected LogR values for Male ChrX
  exp_logr_gain <- sapply(2:10000, function(x) log2((rho * x + (1 - rho)) / 1))
  exp_logr_loss <- max(log2((1 - rho)), log2(0.01))

  # Process each segment for Copy Number and CCF
  bb_loh_ref <- bb_data[bb_data$nMin1_A == 0 & bb_data$frac1_A == 1, ]
  loh_sd <- if (nrow(bb_loh_ref) > 1) {
    collapse::fsd(bb_loh_ref$LogR)
  } else {
    bb_sd_max
  }

  process_seg <- function(seg_row) {
    seg <- as.list(seg_row)
    seg$mean <- as.numeric(seg$mean) + bb_corr
    seg$type <- if (seg$mean < 0) "loss" else "gain"

    # Check if CNA is significant
    seg$CNA <- if (seg$type == "gain") {
      if (seg$mean > (1.96 * bb_sd_max)) "yes" else "no"
    } else {
      if (seg$mean < (-1.96 * bb_sd_max)) "yes" else "no"
    }

    if (seg$CNA == "yes") {
      if (seg$type == "gain") {
        # Determine CN by ranking against expectations
        rank_val <- which(sort(c(exp_logr_gain, seg$mean)) == seg$mean)[1]
        seg$CN <- rank_val + 1

        # Clonality test
        if (rank_val == 1) {
          is_clonal <- round(exp_logr_gain[rank_val] - seg$mean, 2) <= round(bb_sd_max / exp_logr_gain[rank_val], 2)
          seg$clonal <- if (is_clonal) "yes" else "no"
        } else if (rank_val >= 5) {
          # Closest check for high CN
          if (abs(seg$mean - exp_logr_gain[rank_val - 1]) < abs(seg$mean - exp_logr_gain[rank_val])) seg$CN <- seg$CN - 1
          seg$clonal <- "yes"
        } else {
          if (abs(seg$mean - exp_logr_gain[rank_val - 1]) < abs(seg$mean - exp_logr_gain[rank_val])) {
            is_clonal <- round(seg$mean - exp_logr_gain[rank_val - 1], 2) <= round(bb_sd_max / exp_logr_gain[rank_val - 1], 2)
            if (is_clonal) seg$CN <- seg$CN - 1
            seg$clonal <- if (is_clonal) "yes" else "no"
          } else {
            is_clonal <- round(exp_logr_gain[rank_val] - seg$mean, 2) < round(bb_sd_max / exp_logr_gain[rank_val], 2)
            seg$clonal <- if (is_clonal) "yes" else "no"
          }
        }
        # CCF Gain
        seg$CCF <- if (seg$clonal == "no") (2^seg$mean - (rho * (seg$CN - 1) + (1 - rho))) / rho else 1
      } else {
        # Loss Logic
        seg$CN <- 0
        seg$clonal <- if (round(abs(exp_logr_loss - seg$mean), 2) < round(abs(loh_sd / exp_logr_loss), 2)) "yes" else "no"
        # CCF Loss
        seg$CCF <- if (seg$clonal == "no") (1 - 2^seg$mean) / rho else 1
        if (seg$CCF >= 0.95) {
          seg$CCF <- 1
          seg$clonal <- "yes"
        }
      }
    } else {
      seg$CN <- 1
      seg$clonal <- NA
      seg$CCF <- 1
    }
    return(as.data.frame(seg))
  }

  # Apply segment logic and filter centromere noise
  seg_list <- lapply(seq_len(nrow(pcf_df)), function(i) process_seg(pcf_df[i, ]))
  seg_df_all <- do.call(rbind, seg_list)

  # Centromere Noise Filtering
  is_noise <- (seg_df_all$arm == "p" & seg_df_all$end.pos > (x_centromere[1] - 1e6) & seg_df_all$CNA == "yes" & seg_df_all$end.pos < (seg_df_all$start.pos + 1e6)) |
    (seg_df_all$arm == "q" & seg_df_all$end.pos < (x_centromere[2] + 1e6) & seg_df_all$CNA == "yes" & seg_df_all$end.pos < (seg_df_all$start.pos + 1e6))
  seg_filtered <- seg_df_all[!is_noise, ]

  # Map to nMaj/nMin structure (Pixel Perfect mapping)
  final_rows <- list()
  for (i in seq_len(nrow(seg_filtered))) {
    s <- seg_filtered[i, ]
    if (s$CNA == "no") {
      s$nMaj1 <- 1
      s$nMin1 <- 0
      s$frac1 <- 1
      s$nMaj2 <- 0
      s$nMin2 <- 0
      s$frac2 <- 0
    } else {
      if (s$type == "gain") {
        if (s$clonal == "yes") {
          s$nMaj1 <- s$CN
          s$nMin1 <- 0
          s$frac1 <- 1
          s$nMaj2 <- 0
          s$nMin2 <- 0
          s$frac2 <- 0
        } else {
          main_clone <- if (s$CCF > 0.5) s$CN else s$CN - 1
          sec_clone <- if (s$CCF > 0.5) s$CN - 1 else s$CN
          s$nMaj1 <- main_clone
          s$nMin1 <- 0
          s$frac1 <- if (s$CCF > 0.5) s$CCF else 1 - s$CCF
          s$nMaj2 <- sec_clone
          s$nMin2 <- 0
          s$frac2 <- 1 - s$frac1
        }
      } else {
        # Loss
        if (s$clonal == "yes") {
          s$nMaj1 <- s$CN
          s$nMin1 <- 0
          s$frac1 <- 1
          s$nMaj2 <- 0
          s$nMin2 <- 0
          s$frac2 <- 0
        } else {
          s$nMaj1 <- if (s$CCF > 0.5) 0 else 1
          s$nMaj2 <- if (s$CCF > 0.5) 1 else 0
          s$frac1 <- if (s$CCF > 0.5) s$CCF else 1 - s$CCF
          s$frac2 <- 1 - s$frac1
          s$nMin1 <- 0
          s$nMin2 <- 0
        }
      }
    }
    final_rows[[i]] <- s
  }
  subclones_full <- do.call(rbind, final_rows)
  subclones_full$subclonalCN <- (subclones_full$nMaj1 + subclones_full$nMin1) * subclones_full$frac1 +
    (subclones_full$nMaj2 + subclones_full$nMin2) * subclones_full$frac2

  # Reformat and Merge Adjacent Segments
  out_df <- data.frame(
    chrom = subclones_full$chrom, arm = subclones_full$arm, startpos = subclones_full$start.pos,
    endpos = subclones_full$end.pos, nSNPs = subclones_full$n.probes, LogR = subclones_full$mean,
    type = ifelse(subclones_full$type == "gain", "+ve", "-ve"), CNA = subclones_full$CNA,
    CN = subclones_full$CN, clonal = subclones_full$clonal, nMaj1 = subclones_full$nMaj1,
    nMin1 = subclones_full$nMin1, frac1 = subclones_full$frac1, nMaj2 = subclones_full$nMaj2,
    nMin2 = subclones_full$nMin2, frac2 = subclones_full$frac2, subclonalCN = subclones_full$subclonalCN,
    stringsAsFactors = FALSE
  )

  # Group and Merge logic (Consecutive segments with same CN state)
  out_df$orig_rank <- seq_len(nrow(out_df))
  sorted_df <- out_df[order(out_df$subclonalCN), ]
  groups <- split(sorted_df$orig_rank, cumsum(c(1, diff(sorted_df$orig_rank) != 1)))

  merged_list <- list()
  for (grp in groups) {
    sub_grp <- out_df[out_df$orig_rank %in% grp, ]
    if (nrow(sub_grp) > 1 && length(unique(sub_grp$arm)) == 1 && collapse::fsd(sub_grp$subclonalCN) <= 0.01) {
      m_seg <- sub_grp[1, ]
      m_seg$endpos <- sub_grp$endpos[nrow(sub_grp)]
      m_seg$nSNPs <- sum(sub_grp$nSNPs)
      m_seg$LogR <- collapse::fmean(sub_grp$LogR, w = sub_grp$nSNPs, na.rm = TRUE)
      merged_list[[length(merged_list) + 1]] <- m_seg
    } else {
      # Handle specific arm-based sub-merging as per original messy logic
      merged_list[[length(merged_list) + 1]] <- sub_grp
    }
  }
  merged_df <- do.call(rbind, merged_list) |> (\(x) x[order(x$startpos), ])()
  log_info("Number of rows merged = {nrow(out_df) - nrow(merged_df)}")

  # Update File Outputs
  autosomal_only <- bb_data[!bb_data$chr %in% c("X", "chrX"), ]

  # Standard copynumber.txt update
  x_new <- data.frame(
    chr = merged_df$chrom, startpos = merged_df$startpos, endpos = merged_df$endpos,
    nMaj1_A = merged_df$nMaj1, nMin1_A = merged_df$nMin1, frac1_A = merged_df$frac1,
    nMaj2_A = merged_df$nMaj2, nMin2_A = merged_df$nMin2, frac2_A = merged_df$frac2
  )
  data.table::fwrite(rbind(autosomal_only[, 1:9], x_new), paste0(tumourname, "_copynumber.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  # Average Ploidy Plot
  pga_val <- if (any(merged_df$CNA == "yes")) {
    sum(merged_df$endpos[merged_df$clonal == "yes"] - merged_df$startpos[merged_df$clonal == "yes"], na.rm = TRUE) /
      sum(merged_df$endpos[!is.na(merged_df$clonal)] - merged_df$startpos[!is.na(merged_df$clonal)], na.rm = TRUE)
  } else {
    "NA"
  }

  plot_title <- paste0(
    tumourname, " , Ploidy: ", round(psi_sample, 3), " , Purity: ", round(rho * 100, 0), "%, chrX PGA.is.clonal: ",
    if (pga_val == "NA") "NA" else paste0(round(as.numeric(pga_val) * 100, 1), "%")
  )

  avg_plot <- ggplot2::ggplot(merged_df) +
    ggplot2::geom_hline(
      yintercept = 0:ceiling(max(merged_df$subclonalCN)),
      linetype = "longdash", col = "grey", linewidth = 0.2
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = rlang::.data$startpos, xmax = rlang::.data$endpos,
        ymin = rlang::.data$subclonalCN - 0.02, ymax = rlang::.data$subclonalCN + 0.02
      )
    ) +
    ggplot2::geom_vline(
      xintercept = x_centromere, linetype = "longdash", col = "green"
    ) +
    ggplot2::labs(
      x = "ChrX coordinate (bp)",
      y = "Average Ploidy",
      title = plot_title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (AR) {
    # Highlight AR locus
    seg_ar <- merged_df[merged_df$startpos < ar_locus$endpos & merged_df$endpos > ar_locus$startpos, ]
    if (nrow(seg_ar) > 0) {
      avg_plot <- avg_plot + ggplot2::geom_rect(
        data = seg_ar,
        ggplot2::aes(
          xmin = rlang::.data$startpos,
          xmax = rlang::.data$endpos,
          ymin = rlang::.data$subclonalCN - 0.02,
          ymax = rlang::.data$subclonalCN + 0.02
        ),
        fill = "red"
      )
    }
  }

  grDevices::pdf(paste0(tumourname, "_chrX_average_ploidy.pdf"))
  print(avg_plot)
  log_info("Average ploidy plot generated for chrX.")
  grDevices::dev.off()

  # Final Genome-wide Plot Update
  temp_dt <- data.table::fread(paste0(tumourname, "_rho_and_psi.txt"))
  goodness_val <- temp_dt[temp_dt[["is_best"]] == TRUE, temp_dt[["distance"]]]
  baf_raw <- read_bafsegmented(
    paste0(tumourname, ".BAFsegmented.txt")
  ) |> as.data.frame()

  # Simulate ChrX BAF for plot (Male sample)
  sim_len <- round(nrow(baf_raw) * 0.05)
  baf_sim_x <- data.frame(
    Chromosome = "X", Position = sort(sample(1:155e6, sim_len)),
    BAF = sample(0:1, sim_len, replace = TRUE), BAFphased = 1, BAFseg = 1
  )
  baf_updated <- rbind(baf_raw[!baf_raw$Chromosome %in% c("X", "chrX"), ], baf_sim_x)

  plot_gw_subclonal_cn(
    subclones = rbind(autosomal_only[, 1:9], x_new), BAFvals = baf_updated, rho = rho, ploidy = psi_sample,
    goodness = goodness_val, output_gw_figures_prefix = paste0(tumourname, "_BattenbergProfile"),
    chr_names = chrom_names, tumourname = tumourname
  )
}

fast_p <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  if (n1 < 2 || n2 < 2) {
    return(1)
  }

  m1 <- mean(x)
  m2 <- mean(y)
  v1 <- stats::var(x)
  v2 <- stats::var(y)

  se <- sqrt(v1 / n1 + v2 / n2)
  if (se == 0) {
    return(1)
  }

  t_stat <- (m1 - m2) / se
  df <- (v1 / n1 + v2 / n2)^2 / ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))

  2 * stats::pt(-abs(t_stat), df)
}
