#' Fit copy number
#'
#' Function that will fit a clonal copy number profile to segmented data. It first
#' matches the raw LogR with the segmented BAF to create segmented LogR. Then ASCAT
#' is run to obtain a clonal copy number profile. Beyond logRsegmented it produces
#' the rho_and_psi file and the cellularity_ploidy file.
#' @param samplename Samplename used to name the segmented logr output file
#' @param outputfile_prefix Prefix used for all output file names, except logRsegmented
#' @param inputfile_baf_segmented Filename that points to the BAF segmented data
#' @param inputfile_baf Filename that points to the raw BAF data
#' @param inputfile_logr Filename that points to the raw LogR data
#' @param dist_choice The distance metric that is used internally to rank clonal copy number solutions
#' @param ascat_dist_choice The distance metric used to obtain an initial cellularity and ploidy estimate
#' @param min_ploidy The minimum ploidy to consider (Default 1.6)
#' @param max_ploidy The maximum ploidy to consider (Default 4.8)
#' @param min_rho The minimum cellularity to consider (Default 0.1)
#' @param max_rho The maximum cellularity to consider (Default 1.0)
#' @param min_goodness The minimum goodness of fit for a solution to have to be considered (Default 63)
#' @param uninformative_baf_threshold The threshold beyond which BAF becomes uninformative (Default 0.51)
#' @param gamma_param Technology parameter, compaction of Log R profiles. Expected decrease in case of deletion in diploid sample, 100 "\%" aberrant cells; 1 in ideal case, 0.55 of Illumina 109K arrays (Default 1)
#' @param use_preset_rho_psi Boolean whether to use user specified rho and psi values (Default FALSE)
#' @param preset_rho A user specified rho to fit a copy number profile to (Default NA)
#' @param preset_psi A user specified psi to fit a copy number profile to (Default NA)
#' @param read_depth Legacy parameter that is no longer used (Default 30)
#' @param analysis A String representing the type of analysis to be run, this determines whether the distance figure is produced (Default paired)
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
  nthreads,
  enhanced_grid_search = FALSE
) {
  assert_file_exists(inputfile_baf_segmented)
  assert_file_exists(inputfile_baf)
  assert_file_exists(inputfile_logr)
  # Check for enough options supplied for rho and psi
  if ((max_ploidy - min_ploidy) < 0.05) {
    log_failure("Supplied ploidy range must be larger than 0.05: {min_ploidy}-{max_ploidy}")
  }
  if ((max_rho - min_rho) < 0.01) {
    log_failure("Supplied rho range must be larger than 0.01: {min_rho}-{max_rho}")
  }

  # Read in the required data
  segmented.BAF.data <- read_bafsegmented(inputfile_baf_segmented)
  data.table::setDF(segmented.BAF.data)

  raw.BAF.data <- read_baf_as_data_frame(inputfile_baf)
  raw.logR.data <- read_baf_as_data_frame(inputfile_logr)

  # Assign rownames as those are required by various clonal_ascat.R functions
  # If there are duplicates (possible with old versions of BB) then remove those
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

  # For each chromosome
  for (chr in chr_names) {
    chr.BAF.data <- baf_split[[chr]]

    # Skip the rest if there is no data for this chromosome
    if (is.null(chr.BAF.data) || nrow(chr.BAF.data) == 0) {
      next
    }
    # Match segments with chromosome position
    chr.segmented.BAF.data <- baf_segmented_split[[chr]]
    indices <- match(chr.segmented.BAF.data[, 2], chr.BAF.data$Position)

    if (sum(is.na(indices)) == length(indices) || length(indices) == 0) {
      next
    }

    # Drop NAs here too
    chr.segmented.BAF.data <- chr.segmented.BAF.data[!is.na(indices), ]

    # Append the segmented data
    matched.segmented.BAF.data[[chr]] <- chr.segmented.BAF.data
    BAF.data[[chr]] <- chr.BAF.data[indices[!is.na(indices)], ]

    # Append raw LogR
    chr.logR.data <- logr_split[[chr]]
    indices <- match(chr.segmented.BAF.data[, 2], chr.logR.data$Position)
    logR.data[[chr]] <- chr.logR.data[indices[!is.na(indices)], ]
    chr.segmented.logR.data <- chr.logR.data[indices[!is.na(indices)], ]

    # Append segmented LogR
    segs <- rle(chr.segmented.BAF.data[, 5])$lengths
    cum.segs <- c(0, cumsum(segs))
    for (s in seq_along(segs)) {
      chr.segmented.logR.data[(cum.segs[s] + 1):cum.segs[s + 1], 3] <- mean(chr.segmented.logR.data[(cum.segs[s] + 1):cum.segs[s + 1], 3], na.rm = TRUE)
    }
    segmented.logR.data[[chr]] <- chr.segmented.logR.data
  }

  # Sync the dataframes
  selection <- c()
  for (chrom in chr_names) {
    matched.segmented.BAF.data.chr <- matched.segmented.BAF.data[[chrom]] # matched.segmented.BAF.data[matched.segmented.BAF.data[,1]==chrom,]
    logR.data.chr <- logR.data[[chrom]] # logR.data[logR.data[,1]==chrom,]

    selection <- matched.segmented.BAF.data.chr[, 2] %in% logR.data.chr[, 2]
    matched.segmented.BAF.data[[chrom]] <- matched.segmented.BAF.data.chr[selection, ]
    segmented.logR.data[[chrom]] <- segmented.logR.data[[chrom]][selection, ]
  }

  # Combine the split data frames into a single for the subsequent steps
  matched.segmented.BAF.data <- do.call(rbind, matched.segmented.BAF.data)
  matched.segmented.BAF.data <- data.table::rbindlist(matched.segmented.BAF.data)
  segmented.logR.data <- do.call(rbind, segmented.logR.data)
  BAF.data <- do.call(rbind, BAF.data)
  logR.data <- do.call(rbind, logR.data)
  names(matched.segmented.BAF.data)[5] <- samplename

  # write out the segmented logR data
  row.names(segmented.logR.data) <- row.names(matched.segmented.BAF.data)
  row.names(logR.data) <- row.names(matched.segmented.BAF.data)
  data.table::fwrite(segmented.logR.data, paste(samplename, ".logRsegmented.txt", sep = ""), sep = "\t", quote = FALSE, col_names = FALSE, row.names = FALSE)

  # Prepare the data for going into the runASCAT functions
  segBAF <- 1 - matched.segmented.BAF.data[, 5]
  segLogR <- segmented.logR.data[, 3]
  logR <- logR.data[, 3]
  names(segBAF) <- rownames(matched.segmented.BAF.data)
  names(segLogR) <- rownames(matched.segmented.BAF.data)
  names(logR) <- rownames(matched.segmented.BAF.data)

  chr_segs <- NULL
  for (ch in seq_along(chr_names)) {
    chr_segs[[ch]] <- which(logR.data[, 1] == chr_names[ch])
  }

  if (use_preset_rho_psi) {
    ascat_optimum_pair <- list(rho = preset_rho, psi = preset_psi, ploidy = preset_psi)
  } else {
    distance_outfile <- paste(outputfile_prefix, "distance.png", sep = "", collapse = "")
    copynumberprofile_outfile <- paste(outputfile_prefix, "copynumberprofile.png", sep = "", collapse = "")
    nonroundedprofile_outfile <- paste(outputfile_prefix, "nonroundedprofile.png", sep = "", collapse = "")
    cnaStatusFile <- paste(outputfile_prefix, "copynumber_solution_status.txt", sep = "", collapse = "")

    if (enhanced_grid_search) {
      ascat_optimum_pair <- runASCAT_enhanced(
        logR, 1 - BAF.data[, 3], segLogR, segBAF,
        chr_segs, ascat_dist_choice, distance_outfile,
        copynumberprofile_outfile, nonroundedprofile_outfile,
        cnaStatusFile = cnaStatusFile, gamma = gamma_param,
        allow100percent = TRUE, reliabilityFile = NA, min_ploidy = min_ploidy,
        max_ploidy = max_ploidy, min_rho = min_rho, max_rho = max_rho,
        min_goodness = min_goodness, chr_names = chr_names, analysis = analysis,
        uninformative_baf_threshold = uninformative_baf_threshold,
        verbose = TRUE
      )
    } else {
      ascat_optimum_pair <- runASCAT(
        logR, 1 - BAF.data[, 3], segLogR, segBAF,
        chr_segs, ascat_dist_choice,
        distancepng = distance_outfile,
        copynumberprofilespng = copynumberprofile_outfile,
        nonroundedprofilepng = nonroundedprofile_outfile,
        cnaStatusFile = cnaStatusFile,
        gamma = gamma_param, allow100percent = TRUE,
        reliabilityFile = NA, min_ploidy = min_ploidy,
        max_ploidy = max_ploidy, min_rho = min_rho, max_rho = max_rho,
        min_goodness = min_goodness, chr_names = chr_names, analysis = analysis,
        uninformative_baf_threshold = uninformative_baf_threshold
      )
    }
  }

  distance_outfile <- paste(outputfile_prefix, "second_distance.png", sep = "", collapse = "")
  copynumberprofile_outfile <- paste(outputfile_prefix, "second_copynumberprofile.png", sep = "", collapse = "")
  nonroundedprofile_outfile <- paste(outputfile_prefix, "second_nonroundedprofile.png", sep = "", collapse = "")

  # All is set up, now run ASCAT to obtain a clonal copynumber profile
  out <- run_clonal_ASCAT(
    logR, 1 - BAF.data[, 3], segLogR, segBAF, chr_segs,
    matched.segmented.BAF.data, ascat_optimum_pair, dist_choice,
    distance_outfile, copynumberprofile_outfile, nonroundedprofile_outfile,
    gamma_param = gamma_param, read_depth, uninformative_baf_threshold,
    allow100percent = TRUE, reliabilityFile = NA, psi_min_initial = min_ploidy,
    psi_max_initial = max_ploidy, rho_min_initial = min_rho,
    rho_max_initial = max_rho, chr_names = chr_names
  )

  ascat_optimum_pair_fraction_of_genome <- out$output_optimum_pair_without_ref
  ascat_optimum_pair_ref_seg <- out$output_optimum_pair
  is_ref_better <- out$is_ref_better

  # Save rho, psi and ploidy for future reference
  rho_psi_output <- data.frame(
    rho = c(ascat_optimum_pair$rho, ascat_optimum_pair_fraction_of_genome$rho, ascat_optimum_pair_ref_seg$rho),
    psi = c(ascat_optimum_pair$psi, ascat_optimum_pair_fraction_of_genome$psi, ascat_optimum_pair_ref_seg$psi),
    ploidy = c(ascat_optimum_pair$ploidy, ascat_optimum_pair_fraction_of_genome$ploidy, ascat_optimum_pair_ref_seg$ploidy),
    distance = c(NA, out$distance_without_ref, out$distance),
    is_best = c(NA, !is_ref_better, is_ref_better),
    row.names = c("ASCAT", "FRAC_GENOME", "REF_SEG")
  )
  data.table::fwrite(rho_psi_output, paste(outputfile_prefix, "rho_and_psi.txt", sep = ""), quote = FALSE, sep = "\t")
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
#' @param logr_file String that points to the raw LogR file to be used in the subclonal copy number figures
#' @param rho_psi_file String pointing to the rho_and_psi file generated by \code{fit_copy_number}
#' @param output_file Filename of the file where the final copy number fit will be written to
#' @param output_figures_prefix Prefix of the filenames for the chromosome specific copy number figures
#' @param output_gw_figures_prefix Prefix of the filenames for the genome wide copy number figures
#' @param chr_names Vector of allowed chromosome names
#' @param masking_output_file Filename of where the masking details need to be written. Masking is performed to remove very high copy number state segments
#' @param max_allowed_state The maximum CN state allowed (Default 250)
#' @param cn_upper_limit The maximum CN that can be called (Default 1000)
#' @param prior_breakpoints_file A two column file with prior breakpoints, possibly from structural variants. This file must contain two columns: chromosome and position. These are used when making the figures
#' @param gamma Technology specific scaling parameter for LogR (Default 1)
#' @param segmentation_gamma Legacy parameter that is no longer used (Default NA)
#' @param siglevel Threshold under which a p-value becomes significant. When it is significant a second copy number state will be fitted (Default 0.05)
#' @param maxdist Slack in BAF space to allow a segment to be off it's optimum before becoming significant. A segment becomes significant very quickly when a breakpoint is missed, this parameter alleviates the effect (Default 0.01)
#' @param noperms The number of permutations to be run when bootstrapping the confidence intervals on the copy number state of each segment (Default 1000)
#' @param seed Seed to set when performing bootstrapping (Default: Current time)
#' @param calc_seg_baf_option Various options to recalculate the BAF of a segment. Options are: 1 - median, 2 - mean, 3 - ifelse median==0|1, mean, median. (Default: 3)
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
  calc_seg_baf_option = 3, verbose_logging = FALSE
) {
  set.seed(seed)

  # Load and calculate initial rho/psi metrics
  res <- load_rho_psi_file(rho_psi_file)
  rho <- res$rho
  psit <- res$psit
  psi <- (rho * psit) + (2 * (1 - rho))
  goodness <- res$goodness

  # Load BAF data and handle possible row-name artifacts ("X")
  BAFvals <- read_bafsegmented(baf_segmented_file) |> as.data.frame()
  if ("X" %in% colnames(BAFvals)) {
    BAFvals <- BAFvals[, -1, with = FALSE]
  }

  # Positional indexing for generalizability: Col 3 = BAF, Col 5 = BAFseg
  BAF <- BAFvals[, 3]
  BAFseg <- BAFvals[, 5]
  SNPpos <- BAFvals[, c(1, 2), drop = FALSE]

  # Load LogR data and handle row-name artifacts
  LogRvals <- read_logr(logr_file) |> as.data.frame()
  if (identical(colnames(LogRvals)[1], "X")) {
    LogRvals <- LogRvals[, -1, drop = FALSE]
  }

  # Create named index vectors for chromosomes
  ctrans <- setNames(seq_along(chr_names), chr_names)
  ctrans.logR <- setNames(seq_along(chr_names), chr_names)

  # First Pass: Determine Copy Number and Merge Segments
  res_cn <- determine_copynumber(
    BAFvals, LogRvals, rho, psi, gamma,
    ctrans, ctrans.logR, maxdist, siglevel, noperms, cn_upper_limit
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
    ctrans, ctrans.logR, maxdist, siglevel,
    noperms, cn_upper_limit
  )
  subcloneres <- res_final$subcloneres
  BAFpvals <- res_final$BAFpvals

  # Mask high CN artifacts
  mask_res <- mask_high_cn_segments(subcloneres, BAFvals, max_allowed_state)
  subcloneres <- mask_res$subclones

  # Output Masking Details
  masking_details <- data.frame(
    samplename = sample_name,
    masked_count = mask_res$masked_count,
    masked_size = mask_res$masked_size,
    max_allowed_state = max_allowed_state
  )
  data.table::fwrite(
    masking_details,
    file = masking_output_file,
    quote = FALSE, sep = "\t", row.names = FALSE
  )

  # Generate output paths
  base_out <- tools::file_path_sans_ext(output_file)
  ext_out <- tools::file_ext(output_file)

  data.table::fwrite(
    subcloneres[, c(1:3, 8:13)], output_file,
    quote = FALSE, sep = "\t", row.names = FALSE
  )
  data.table::fwrite(
    subcloneres, paste0(base_out, "_extended.", ext_out),
    quote = FALSE, sep = "\t", row.names = FALSE
  )

  # Calculate Clonal PGA (Percent Genome Altered)
  subcloneres$length <- subcloneres$endpos - subcloneres$startpos
  diploid_idx <- which(subcloneres$nMaj1_A == 1 & subcloneres$nMin1_A == 1 & subcloneres$frac1_A == 1)

  cna <- if (length(diploid_idx) > 0) subcloneres[-diploid_idx, ] else subcloneres
  subcloneres_subclonal <- subcloneres[subcloneres$frac1_A < 1, ]

  if (nrow(cna) == 0 || sum(cna$length, na.rm = TRUE) == 0 || nrow(subcloneres_subclonal) == 0) {
    goodness <- 1.0
  } else {
    subclonal_fraction <- sum(subcloneres_subclonal$length) / sum(cna$length)
    goodness <- max(0, min(1, 1 - subclonal_fraction))
  }

  message(sprintf("PGA.is.clonal = %2.1f%%", goodness * 100))

  # Visualization
  segment_breakpoints <- collapse_bafsegmented_to_segments(BAFvals)
  has_prior <- !is.null(prior_breakpoints_file) &&
    !is.na(prior_breakpoints_file) &&
    prior_breakpoints_file != "NA"

  if (has_prior) {
    svs <- data.table::fread(prior_breakpoints_file, data.table = FALSE)
  }

  for (chr in chr_names) {
    chr_idx <- SNPpos[, 1] == chr
    pos <- SNPpos[chr_idx, 2]

    if (length(pos) == 0) next

    # Using positional indexing for svs (Col 1: Chr, Col 2: Pos)
    svs_pos <- if (has_prior) svs[svs[, 1] == chr, 2] / 1e6 else NULL

    # Identify breakpoints (Col 1: Chr, Col 2: Start, Col 3: End)
    bp_chr <- segment_breakpoints[segment_breakpoints[, 1] == chr, ]
    breakpoints_pos <- sort(unique(c(bp_chr[, 2], bp_chr[, 3]) / 1e6))

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

  # Clean up and calculate Ploidy
  subclones <- as.data.frame(subcloneres)
  num_cols <- 2:ncol(subclones)
  subclones[num_cols] <- lapply(subclones[num_cols], function(x) as.numeric(as.character(x)))

  seg_len <- floor((subclones$endpos - subclones$startpos) / 1000)

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
#' @param ctrans.logR Named vector of chromosome names
#' @param maxdist Max distance a segment is tolerated to be not considered for subclonal copy number
#' @param siglevel Level at which a segment can become significantly different from the nearest clonal state
#' @param noperms Number of bootstrap permutations
#' @param cn_upper_limit Maximum number of CN that can be called
#' @return A data.frame with copy number determined for each segment
#' @author dw9
#' @noRd
determine_copynumber <- function(
  BAFvals,
  LogRvals,
  rho,
  psi,
  gamma,
  ctrans,
  ctrans.logR,
  maxdist,
  siglevel,
  noperms,
  cn_upper_limit
) {
  # Positional extraction to maintain generalizability
  BAFphased <- BAFvals[, 4]
  BAFseg <- BAFvals[, 5]

  # Large integer coordinate mapping to handle multiple chromosomes in one vector
  # Using 1e9 as a spacer between chromosomes
  BAFpos <- as.vector(ctrans[as.vector(BAFvals[, 1])] * 1e9 + BAFvals[, 2])
  LogRpos <- as.vector(ctrans.logR[as.vector(LogRvals[, 1])] * 1e9 + LogRvals[, 2])

  # Identify segment boundaries where BAF level or Chromosome changes
  # switchpoints identifies the indices of the end of each segment
  n_rows <- nrow(BAFvals)
  boundary_idx <- which(BAFseg[-1] != BAFseg[-n_rows] | BAFvals[-1, 1] != BAFvals[-n_rows, 1])
  switchpoints <- c(0, boundary_idx, length(BAFseg))
  BAFlevels <- BAFseg[switchpoints[-1]]

  # Pre-allocate output containers
  BAFpvals <- vector(mode = "numeric", length = length(BAFseg))
  res_list <- vector(mode = "list", length = length(BAFlevels))

  for (i in seq_along(BAFlevels)) {
    l <- BAFlevels[i]

    # Ensure major/minor orientation (BAF >= 0.5)
    l <- max(l, 1 - l)

    # Segment slicing
    start_idx <- switchpoints[i] + 1
    end_idx <- switchpoints[i + 1]
    BAFke <- BAFphased[start_idx:end_idx]

    # Coordinate extraction
    s_pos <- BAFpos[start_idx:end_idx]
    startpos <- min(s_pos)
    endpos <- max(s_pos)

    # Extract chromosome name (Col 1 is Chromosome)
    chrom <- BAFvals[start_idx, 1]

    # LogR calculation: Mean of non-infinite values within segment range
    # LogRvals[, 3] is the LogR value
    logr_mask <- LogRpos >= startpos & LogRpos <= endpos & !is.infinite(LogRvals[, 3])
    LogR <- mean(LogRvals[logr_mask, 3], na.rm = TRUE)
    if (is.na(LogR)) LogR <- 0

    # Theoretical Copy Number calculation
    nMajor <- (rho - 1 + l * psi * 2^(LogR / gamma)) / rho
    nMinor <- (rho - 1 + (1 - l) * psi * 2^(LogR / gamma)) / rho

    if (is.na(nMinor)) next

    # Handle physical impossibility (negative copy number)
    if (nMinor < 0) {
      if (l == 1) {
        nMajor <- cn_upper_limit
      } else {
        nMajor <- nMajor + l * (0.01 - nMinor) / (1 - l)
      }
      nMinor <- 0.01
    }

    # Corner case testing
    nMaj_corners <- c(floor(nMajor), ceiling(nMajor), floor(nMajor), ceiling(nMajor))
    nMin_corners <- c(ceiling(nMinor), ceiling(nMinor), floor(nMinor), floor(nMinor))
    ntot <- nMajor + nMinor

    levels <- (1 - rho + rho * nMaj_corners) / (2 - 2 * rho + rho * (nMaj_corners + nMin_corners))
    levels[nMaj_corners == 0 & nMin_corners == 0] <- 0.5

    # Prioritize nearest clonal states
    all_edges <- prioritizeCopyNumbers(levels, l, ntot, floor(nMinor), floor(nMajor), full = TRUE)
    nMaj_test <- all_edges[1, c(1, 3)]
    nMin_test <- all_edges[1, c(2, 4)]
    test_levels <- (1 - rho + rho * nMaj_test) / (2 - 2 * rho + rho * (nMaj_test + nMin_test))
    best_idx <- which.min(abs(test_levels - l))

    # Significance testing
    p_val <- if (is.na(
      collapse::fsd(BAFke)
    ) || collapse::fsd(BAFke) == 0) {
      0
    } else {
      t.test(BAFke, mu = test_levels[best_idx])$p_value
    }
    if (abs(l - test_levels[best_idx]) < maxdist) p_val <- 1

    BAFpvals[start_idx:end_idx] <- p_val

    # Standardize result coordinates (remove the 1e9 multiplier)
    out_start <- startpos %% 1e9
    out_end <- endpos %% 1e9

    if (p_val <= siglevel) {
      # Subclonal fitting logic
      na_rows <- which(is.na(rowSums(all_edges)))
      if (length(na_rows) > 0) {
        all_edges <- rbind(all_edges[-na_rows, , drop = FALSE], all_edges[na_rows, , drop = FALSE])
      }

      nMaj1 <- all_edges[, 1]
      nMin1 <- all_edges[, 2]
      nMaj2 <- all_edges[, 3]
      nMin2 <- all_edges[, 4]

      # Fraction calculation
      tau <- (1 - rho + rho * nMaj2 - 2 * l * (1 - rho) - l * rho * (nMin2 + nMaj2)) /
        (l * rho * (nMin1 + nMaj1) - l * rho * (nMin2 + nMaj2) - rho * nMaj1 + rho * nMaj2)

      # Standard error estimation
      sdl <- collapse::fsd(BAFke, na.rm = TRUE) / sqrt(sum(!is.na(BAFke)))
      calc_tau <- function(curr_l) {
        (1 - rho + rho * nMaj2 - 2 * curr_l * (1 - rho) - curr_l * rho * (nMin2 + nMaj2)) /
          (curr_l * rho * (nMin1 + nMaj1) - curr_l * rho * (nMin2 + nMaj2) - rho * nMaj1 + rho * nMaj2)
      }
      sdtau <- (abs(calc_tau(l + sdl) - tau) + abs(calc_tau(l - sdl) - tau)) / 2

      # Bootstrap block
      sdtaubootstrap <- tau25 <- tau975 <- numeric(length(tau))
      for (opt in seq_along(tau)) {
        permFraction <- replicate(noperms, {
          pMean <- mean(sample(BAFke, length(BAFke), replace = TRUE))
          (1 - rho + rho * nMaj2[opt] - 2 * pMean * (1 - rho) - pMean * rho * (nMin2[opt] + nMaj2[opt])) /
            (pMean * rho * (nMaj1[opt] + nMin1[opt]) - pMean * rho * (nMaj2[opt] + nMin2[opt]) - rho * nMaj1[opt] + rho * nMaj2[opt])
        })
        ordered <- sort(permFraction)
        sdtaubootstrap[opt] <- collapse::fsd(permFraction)
        tau25[opt] <- ordered[25]
        tau975[opt] <- ordered[975]
      }

      # Compile 6-option subclonal result
      res_list[[i]] <- c(
        chrom, out_start, out_end, l, p_val, LogR, ntot,
        as.vector(t(cbind(nMaj1, nMin1, tau, nMaj2, nMin2, 1 - tau, sdtau, sdtaubootstrap, tau25, tau975)))[1:60]
      )
    } else {
      # Clonal result
      res_list[[i]] <- c(
        chrom, out_start, out_end, l, p_val, LogR, ntot,
        nMaj_test[best_idx], nMin_test[best_idx], 1, rep(NA, 57)
      )
    }
  }

  # Build dataframe efficiently from list
  subcloneres <- do.call(rbind, res_list) |> as.data.frame()

  # Standardized column names
  col_bases <- c("A", "B", "C", "D", "E", "F")
  col_suffixes <- c("nMaj1", "nMin1", "frac1", "nMaj2", "nMin2", "frac2", "SDfrac", "SDfrac_BS", "frac1_0.025", "frac1_0.975")
  dynamic_cols <- as.vector(t(outer(col_bases, col_suffixes, function(x, y) paste0(y, "_", x))))

  colnames(subcloneres) <- c("chr", "startpos", "endpos", "BAF", "pval", "LogR", "ntot", dynamic_cols)

  # Column coercion (avoiding loop-based factor conversion)
  subcloneres[-1] <- lapply(subcloneres[-1], function(x) as.numeric(as.character(x)))

  return(list(subcloneres = subcloneres, BAFpvals = BAFpvals))
}


#' Merge copy number segments
#'
#' Merges segments if there is not enough evidence for them to be separate. Two adjacent segments are merged
#' when they are either fit with the same clonal copy number state or when their BAF is not significantly different
#' and their logR puts them in the same square.
#' @param subclones A completely fit copy number profile in Battenberg output format
#' @param bafsegmented A BAFsegmented data.frame with the 5 columns that corresponds to the subclones file
#' @param logR The raw logR data
#' @param rho The rho estimate that the profile was fit with
#' @param psi the psi estimate that the profile was fit with
#' @param platform_gamma The gamma parameter for this platform
#' @param calc_seg_baf_option Various options to recalculate the BAF of a segment. Options are: 1 - median, 2 - mean, 3 - ifelse median== 0|1, mean, median. (Default: 3)
#' @param verbose A boolean to show merging operations (Default: FALSE)
#' @return A list with two fields: bafsegmented and subclones. The subclones field contains a data.frame in
#' Battenberg output format with the merged segments. The bafsegmented field contains the BAFsegmented data
#' corresponding to the provided subclones data.frame.
#' @author sd11, tl
#' @noRd
merge_segments <- function(
  subclones,
  bafsegmented,
  logR,
  rho,
  psi,
  platform_gamma,
  calc_seg_baf_option = 3,
  verbose_logging = FALSE
) {
  calc_nmin <- function(rho, psi, baf, logr, platform_gamma) {
    return((rho - 1 - (baf - 1) * 2^(logr / platform_gamma) * ((1 - rho) * 2 + rho * psi)) / rho)
  }
  calc_nmaj <- function(rho, psi, baf, logr, platform_gamma) {
    return((rho - 1 + baf * 2^(logr / platform_gamma) * ((1 - rho) * 2 + rho * psi)) / rho)
  }
  # Convert DF into GRanges objects
  df2gr <- function(DF, chr, pos1, pos2) {
    return(GenomicRanges::makeGRangesFromDataFrame(
      df = DF,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = NULL,
      seqnames.field = chr,
      start.field = pos1,
      end.field = pos2,
      starts.in.df.are.0based = FALSE
    ))
  }
  # Function called when two segments have not been merged so there is no need to recheck those again
  update_neighbour <- function(subclones, INDEX, INDEX_N) {
    if (INDEX_N > INDEX) {
      subclones$Next_checked[INDEX] <- TRUE
      subclones$Prev_checked[INDEX_N] <- TRUE
    } else {
      subclones$Prev_checked[INDEX] <- TRUE
      subclones$Next_checked[INDEX_N] <- TRUE
    }
    return(subclones)
  }
  # Function called when two segments have been merged so we need to recheck its two neighbours
  updateAround <- function(subclones, INDEX) {
    if (INDEX > 1) {
      subclones$Prev_checked[INDEX] <- FALSE
      subclones$Next_checked[INDEX - 1] <- FALSE
    } else {
      subclones$Prev_checked[INDEX] <- TRUE
    }
    if (INDEX < length(subclones)) {
      subclones$Next_checked[INDEX] <- FALSE
      subclones$Prev_checked[INDEX + 1] <- FALSE
    } else {
      subclones$Next_checked[INDEX] <- TRUE
    }
    return(subclones)
  }
  # Function called to test whether two segments must be checked
  check_status <- function(subclones, INDEX, INDEX_N) {
    if (INDEX_N > INDEX) {
      # Largest segment (INDEX_N) is after smallest one (INDEX)
      stopifnot(subclones$Next_checked[INDEX] == subclones$Prev_checked[INDEX_N])
      if (subclones$Next_checked[INDEX] && subclones$Prev_checked[INDEX_N]) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      # Largest segment (INDEX_N) is before smallest one (INDEX)
      stopifnot(subclones$Prev_checked[INDEX] == subclones$Next_checked[INDEX_N])
      if (subclones$Prev_checked[INDEX] && subclones$Next_checked[INDEX_N]) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }

  # Function to merge two segments
  merge_seg <- function(
    subclones, bafsegmented,
    logR, INDEX, INDEX_N,
    calc_seg_baf_option
  ) {
    # Standard GenomicRanges coordinate updates
    if (INDEX_N < INDEX) {
      GenomicRanges::end(
        subclones[INDEX_N]
      ) <- GenomicRanges::end(subclones[INDEX])
    } else {
      GenomicRanges::start(
        subclones[INDEX_N]
      ) <- GenomicRanges::start(subclones[INDEX])
    }

    # Remove the merged-from segment
    subclones <- subclones[-INDEX]
    if (INDEX_N < INDEX) INDEX <- INDEX - 1

    # Trigger local neighbor update logic
    subclones <- updateAround(subclones, INDEX)

    # Efficient overlap extraction
    # subjectHits is the linter-safe version of @to
    baf_idx <- S4Vectors::subjectHits(
      GenomicRanges::findOverlaps(subclones[INDEX], bafsegmented)
    )
    baf_vals <- bafsegmented$BAFphased[baf_idx]

    # Modernized BAF calculation with safety for NA values
    if (calc_seg_baf_option == 1) {
      NEW_BAF <- collapse::fmedian(baf_vals, na.rm = TRUE)
    } else if (calc_seg_baf_option == 2) {
      NEW_BAF <- collapse::fmean(baf_vals, na.rm = TRUE)
    } else if (calc_seg_baf_option == 3) {
      # Calculate both using high-performance C++ bindings
      m_baf <- collapse::fmedian(baf_vals, na.rm = TRUE)

      # Robust Logic: Only use the median if it's not NA
      # This avoids the "missing value where TRUE/FALSE needed" error
      if (!base::is.na(m_baf) && m_baf != 0 && m_baf != 1) {
        NEW_BAF <- m_baf
      } else {
        NEW_BAF <- collapse::fmean(baf_vals, na.rm = TRUE)
      }
    }

    # LogR update with safety for empty segments
    logr_idx <- S4Vectors::subjectHits(
      GenomicRanges::findOverlaps(subclones[INDEX], logR)
    )

    if (base::length(logr_idx) == 0) {
      subclones[INDEX]$LogR <- 0
    } else {
      subclones[INDEX]$LogR <- collapse::fmean(
        logR$logR[logr_idx],
        na.rm = TRUE
      )
    }

    # Update metadata on the S4 objects
    subclones[INDEX]$BAF <- NEW_BAF
    bafsegmented$BAFseg[baf_idx] <- NEW_BAF

    # Standard Evaluation sequence generation
    subclones$ID <- base::seq_along(subclones)

    list(subclones = subclones, bafsegmented = bafsegmented)
  }

  log_debug("Converting DFs into GRanges objects")

  subclones <- subclones |>
    df2gr("chr", "startpos", "endpos") |>
    GenomicRanges::sort()

  bafsegmented <- bafsegmented |>
    df2gr("Chromosome", "Position", "Position") |>
    GenomicRanges::sort()

  logR <- logR |>
    df2gr("Chromosome", "Position", "Position") |>
    GenomicRanges::sort()
  names(GenomicRanges::mcols(logR)) <- "logR"

  # Get unique chromosomes
  chr_names <- unique(as.character(GenomicRanges::seqnames(bafsegmented)))

  # Split by chromosome
  subclones <- split(subclones, GenomicRanges::seqnames(subclones))
  bafsegmented <- split(bafsegmented, GenomicRanges::seqnames(bafsegmented))
  logR <- split(logR, GenomicRanges::seqnames(logR))

  if (!all(chr_names %in% names(subclones)) || !all(chr_names %in% names(bafsegmented)) || !all(chr_names %in% names(logR))) {
    log_failure("Missing data for some chromosomes in one or more inputs")
  }

  # Process each chromosome
  for (CHR in chr_names) {
    log_debug("Merging segments within: {CHR}")

    subclones_chr <- subclones[[CHR]]
    bafsegmented_chr <- bafsegmented[[CHR]]
    logR_chr <- logR[[CHR]]

    # Initialize tracking columns
    subclones_chr$ID <- seq_along(subclones_chr)
    subclones_chr$prev_checked <- FALSE
    subclones_chr$next_checked <- FALSE
    subclones_chr$prev_checked[1] <- TRUE
    subclones_chr$next_checked[length(subclones_chr)] <- TRUE

    while (TRUE) {
      # Find segments needing checks
      unchecked <- which(!subclones_chr$prev_checked | !subclones_chr$next_checked)
      if (length(unchecked) == 0) break

      # Select smallest unchecked segment
      widths <- GenomicRanges::width(subclones_chr[unchecked])
      index <- unchecked[which.min(widths)]

      log_debug("Working on segment: {index} ({subclones_chr[index]})")

      # Determine possible neighbors
      n <- length(subclones_chr)
      neighbors <- integer(0)
      if (index > 1) neighbors <- c(neighbors, index - 1)
      if (index < n) neighbors <- c(neighbors, index + 1)

      if (length(neighbors) == 0) next

      # Sort neighbors by distance (closest first)
      dists <- GenomicRanges::distance(subclones_chr[index], subclones_chr[neighbors])
      sorted_neighbors <- neighbors[order(dists)]

      merged <- FALSE
      for (index_n in sorted_neighbors) {
        log_debug("Checking neighbour: {index_n} ({subclones_chr[index_n]}; distance={dists[which(neighbors == index_n)]})")

        # Skip if already checked
        if (check_status(subclones_chr, index, index_n)) {
          log_debug("Already checked")
          next
        }

        # Check distance threshold
        if (GenomicRanges::distance(subclones_chr[index], subclones_chr[index_n]) > 3e6) {
          log_debug("Distance > 3Mb - do not merge")
          subclones_chr <- update_neighbour(subclones_chr, index, index_n)
          next
        }

        # Check for identical clonal CN
        if (subclones_chr$nMaj1_A[index] == subclones_chr$nMaj1_A[index_n] &&
          subclones_chr$nMin1_A[index] == subclones_chr$nMin1_A[index_n] &&
          subclones_chr$frac1_A[index] == 1 &&
          subclones_chr$frac1_A[index_n] == 1) {
          log_debug("Same clonal CN solution - merge")
          res <- merge_seg(subclones_chr, bafsegmented_chr, logR_chr, index, index_n, calc_seg_baf_option)
          subclones_chr <- res$subclones
          bafsegmented_chr <- res$bafsegmented
          merged <- TRUE
          break
        }

        # Check for compatible CN via stats
        log_debug("Different CN solutions: check BAF and logR")
        nmin_curr <- round(calc_nmin(rho, psi, subclones_chr$BAF[index], subclones_chr$LogR[index], platform_gamma))
        nmaj_curr <- round(calc_nmaj(rho, psi, subclones_chr$BAF[index], subclones_chr$LogR[index], platform_gamma))
        nmin_other <- round(calc_nmin(rho, psi, subclones_chr$BAF[index_n], subclones_chr$LogR[index_n], platform_gamma))
        nmaj_other <- round(calc_nmaj(rho, psi, subclones_chr$BAF[index_n], subclones_chr$LogR[index_n], platform_gamma))

        if (nmin_curr == nmin_other || nmaj_curr == nmaj_other) {
          # Check sufficient data points
          logr_curr <- logR_chr$logR[GenomicRanges::findOverlaps(subclones_chr[index], logR_chr)@to]
          logr_other <- logR_chr$logR[GenomicRanges::findOverlaps(subclones_chr[index_n], logR_chr)@to]
          baf_curr <- bafsegmented_chr$BAFphased[GenomicRanges::findOverlaps(subclones_chr[index], bafsegmented_chr)@to]
          baf_other <- bafsegmented_chr$BAFphased[GenomicRanges::findOverlaps(subclones_chr[index_n], bafsegmented_chr)@to]

          if (sum(!is.na(logr_curr)) > 10 && sum(!is.na(logr_other)) > 10 &&
            sum(!is.na(baf_curr)) > 10 && sum(!is.na(baf_other)) > 10) {
            logr_p <- t.test(logr_curr, logr_other)$p_value
            baf_p <- t.test(baf_curr, baf_other)$p_value
            if (logr_p >= 0.05 && baf_p >= 0.05) {
              log_debug("No significant difference - merge")
              res <- merge_seg(subclones_chr, bafsegmented_chr, logR_chr, index, index_n, calc_seg_baf_option)
              subclones_chr <- res$subclones
              bafsegmented_chr <- res$bafsegmented
              merged <- TRUE
              break
            } else {
              log_debug("Significant difference - do not merge")
              subclones_chr <- update_neighbour(subclones_chr, index, index_n)
            }
          } else {
            log_debug("Too few values - do not merge")
            subclones_chr <- update_neighbour(subclones_chr, index, index_n)
          }
        } else {
          log_debug("Different squares - do not merge")
          subclones_chr <- update_neighbour(subclones_chr, index, index_n)
        }
      }
      if (merged) next # Continue while loop after merge
    }

    # Store back processed data
    subclones[[CHR]] <- subclones_chr
    bafsegmented[[CHR]] <- bafsegmented_chr
  }

  log_debug("Convert GRanges objects into DFs")

  # Combine and convert to data frames
  bafsegmented <- data.frame(Reduce(c, bafsegmented), stringsAsFactors = FALSE)[, -c(3:5)]
  bafsegmented$seqnames <- as.character(bafsegmented$seqnames)
  colnames(bafsegmented)[1:2] <- c("Chromosome", "Position")

  subclones <- data.frame(Reduce(c, subclones), stringsAsFactors = FALSE)[, -c(4:5)]
  subclones$seqnames <- as.character(subclones$seqnames)
  colnames(subclones)[1:3] <- c("chr", "startpos", "endpos")
  subclones$ID <- NULL
  subclones$prev_checked <- NULL
  subclones$next_checked <- NULL

  return(list(bafsegmented = bafsegmented, subclones = subclones))
}

#' Mask segments that have a too high CN state
#' @param subclones Subclones output data
#' @param bafsegmented BAFsegmented data
#' @param max_allowed_state The maximum state allowed before overruling takes place
#' @return A list with the masked subclones, bafsegmented and the number of segments masked and their total genome size
#' @author sd11
mask_high_cn_segments <- function(subclones, bafsegmented, max_allowed_state) {
  count <- 0
  masked_size <- 0
  for (i in seq_len(nrow(subclones))) {
    if (subclones$nMaj1_A[i] > max_allowed_state || subclones$nMin1_A[i] > max_allowed_state) {
      # Mask this segment
      subclones[i, "nMaj1_A"] <- NA
      subclones[i, "nMin1_A"] <- NA
      subclones[i, "nMaj2_A"] <- NA
      subclones[i, "nMin2_A"] <- NA
      # Mask the BAFsegmented
      bafsegmented[subclones$chr[i] == bafsegmented$Chromosome & subclones$startpos[i] < bafsegmented$Position & subclones$endpos[i] >= bafsegmented$Position, c("BAFseg")] <- NA
      count <- count + 1
      masked_size <- masked_size + (subclones$endpos[i] - subclones$startpos[i])
    }
  }
  return(list(subclones = subclones, bafsegmented = bafsegmented, masked_count = count, masked_size = masked_size))
}


#' Plot the copy number genome wide in two different ways. This creates the Battenberg average
#' profile where subclonal copy number is represented as a mixture of two different states and
#' the Battenberg subclones profile where subclonal copy number is plotted as two different
#' separate states. The thickness of the line represents the fraction of tumour cells carying
#' the particular state.
#' @noRd
plot_gw_subclonal_cn <- function(subclones, BAFvals, rho, ploidy, goodness, output_gw_figures_prefix, chr_names, tumourname) {
  # Map start and end of each segment into the BAF values. The plot uses the index of this BAF table as x-axis
  pos_min <- array(NA, nrow(subclones))
  pos_max <- array(NA, nrow(subclones))
  for (i in seq_len(nrow(subclones))) {
    segm_chr <- subclones$chr[i] == BAFvals$Chromosome & subclones$startpos[i] < BAFvals$Position & subclones$endpos[i] >= BAFvals$Position
    pos_min[i] <- min(which(segm_chr))
    pos_max[i] <- max(which(segm_chr))
  }

  # For those segments that are subclonal, Obtain the second state.
  is_subclonal <- which(subclones$frac1_A < 1)
  subcl_min <- array(NA, length(is_subclonal))
  subcl_max <- array(NA, length(is_subclonal))
  for (i in seq_along(is_subclonal)) {
    segment_index <- is_subclonal[i]
    segm_chr <- subclones$chr[segment_index] == BAFvals$Chromosome & subclones$startpos[segment_index] < BAFvals$Position & subclones$endpos[segment_index] >= BAFvals$Position
    subcl_min[i] <- min(which(segm_chr))
    subcl_max[i] <- max(which(segm_chr))
  }

  # Determine whether it's the major or the minor allele that is represented by two states
  is_subclonal_maj <- abs(subclones$nMaj1_A - subclones$nMaj2_A) > 0
  is_subclonal_min <- abs(subclones$nMin1_A - subclones$nMin2_A) > 0
  is_subclonal_maj[is.na(is_subclonal_maj)] <- FALSE
  is_subclonal_min[is.na(is_subclonal_min)] <- FALSE

  segment_states_min <- subclones$nMin1_A * ifelse(is_subclonal_min, subclones$frac1_A, 1) + ifelse(is_subclonal_min, subclones$nMin2_A, 0) * ifelse(is_subclonal_min, subclones$frac2_A, 0)
  segment_states_maj <- subclones$nMaj1_A * ifelse(is_subclonal_maj, subclones$frac1_A, 1) + ifelse(is_subclonal_maj, subclones$nMaj2_A, 0) * ifelse(is_subclonal_maj, subclones$frac2_A, 0)
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
  allele_ratio_plot(samplename = samplename, logr = logr, bafsegmented = bafsegmented, logrsegmented = logrsegmented, outputfile = outputfile, max.plot.cn = 8)

  if (!is.null(allelecounts_file)) {
    allelecounts <- as.data.frame(read_table_generic(allelecounts_file))
    outputfile <- paste0(samplename, "_coverage.png")
    coverage_plot(samplename, allelecounts, outputfile)
  }
}


#' Fit ChrX subclonal copy number (male only)
#'
#' Function to call ChrX copy number based on LogR (suitable for male samples). Copy number
#' cannot be called for the non-PAR region of ChrX due to the hemizygosity of all 1000G SNPs.
#' This function enables calling subclonal copy number for the non-PAR region by segmenting LogR.
#' A number of correction steps are undertaken to account for the noisy nature of LogR. This function
#' requires the following libraries: copynumber, data.table and ggplot2. It reads in three files generated
#' by previous steps of Battenberg, namely samplename_mutantLogR_gcCorrected.tab, samplename_purity_ploidy.txt
#' and samplename_copynumber_extended.txt.
#' This function will also update the Battenberg genome-wide profile plots (average.png and subclones.png) to include the chrX profile by also
#' reading in the samplename.BAFsegmented.txt and samplename_rho_psi.txt files
#' @param tumourname The sample name used for Battenberg (i.e. the tumour BAM file name without the .bam extension)
#' @param X_gamma The PCF gamma value for segmentation of 1000G SNP LogR values (Default 1000)
#' @param X_kmin The min number of SNPs to support a segment in PCF of LogR values (Default 100)
#' @param genomebuild The genome build used in running Battenberg (hg19 or hg38)
#' @param AR Should the segment carrying the androgen receptor (AR) locus to be visually distinguished in average plot? (Default TRUE)
#' @param prior_breakpoints_file A two column text file with prior genome-wide breakpoints, possibly from structural variants. This file must contain two columns with headers "chr" and "pos" representing chromosome and position.
#' @param chrom_names A vector containing the names of chromosomes to be included in the final genome-wide Battenberg copy number plot with chrX
#' @author naser.ansari-pour
#' @export
callChrXsubclones <- function(
  tumourname, X_gamma = 1000,
  X_kmin = 100, genomebuild,
  AR = TRUE, prior_breakpoints_file = NULL,
  chrom_names, data_type = "wgs"
) {
  message(paste("Processing sample:", tumourname))

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
  message(paste("Number of chrX nonPAR SNPs =", nrow(pcf_input)))

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
  message(paste("Number of rows merged =", nrow(out_df) - nrow(merged_df)))

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
        xmin = startpos, xmax = endpos,
        ymin = subclonalCN - 0.02, ymax = subclonalCN + 0.02
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
          xmin = startpos, xmax = endpos,
          ymin = subclonalCN - 0.02,
          ymax = subclonalCN + 0.02
        ),
        fill = "red"
      )
    }
  }

  grDevices::pdf(paste0(tumourname, "_chrX_average_ploidy.pdf"))
  print(avg_plot)
  grDevices::dev.off()

  # Final Genome-wide Plot Update
  temp_dt <- data.table::fread(paste0(tumourname, "_rho_and_psi.txt"))
  goodness_val <- temp_dt[temp_dt[["is_best"]] == TRUE, temp_dt[["distance"]]]
  baf_raw <- read_bafsegmented(
    paste0(tumourname, ".BAFsegmented.txt")
  ) |> as.data.frame()

  # Simulate ChrX BAF for plot (Male sample)
  sim_len <- round(nrow(baf_raw) * 0.05)
  baf_sim_x <- data.frame(Chromosome = "X", Position = sort(sample(1:155e6, sim_len)), BAF = sample(0:1, sim_len, replace = TRUE), BAFphased = 1, BAFseg = 1)
  baf_updated <- rbind(baf_raw[!baf_raw$Chromosome %in% c("X", "chrX"), ], baf_sim_x)

  plot_gw_subclonal_cn(
    subclones = rbind(autosomal_only[, 1:9], x_new), BAFvals = baf_updated, rho = rho, ploidy = psi_sample,
    goodness = goodness_val, output_gw_figures_prefix = paste0(tumourname, "_BattenbergProfile"),
    chr_names = chrom_names, tumourname = tumourname
  )
}
