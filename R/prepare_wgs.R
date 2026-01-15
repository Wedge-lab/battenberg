#' Obtain allele counts for 1000 Genomes loci through external program alleleCount
#'
#' @param bam.file A BAM alignment file on which the counter should be run.
#' @param output_file The file where output should go.
#' @param g1000.loci A file with 1000 Genomes SNP loci.
#' @param min.base.qual The minimum base quality required for it to be counted (optional, default=20).
#' @param min.map.qual The minimum mapping quality required for it to be counted (optional, default=35).
#' @param allelecounter.exe A pointer to where the alleleCounter executable can be found (optional, default points to $PATH).
#' @author sd11
#' @export
getAlleleCounts <- function(bam.file, output_file, g1000.loci, min.base.qual = 20, min.map.qual = 35, allelecounter.exe = "alleleCounter") {
  cmd <- paste(
    allelecounter.exe,
    "-b", bam.file,
    "-l", g1000.loci,
    "-o", output_file,
    "-m", min.base.qual,
    "-q", min.map.qual
  )


  # alleleCount >= v4.0.0 is sped up considerably on 1000G loci when run in dense-snp mode
  counter_version <- system(paste(allelecounter.exe, "--version"), intern = TRUE)
  if (as.integer(substr(x = counter_version, start = 1, stop = 1)) >= 4) {
    cmd <- paste(cmd, "--dense-snps")
  }

  exit_code <- system(cmd, wait = TRUE)
  stopifnot(exit_code == 0)
}


#' Obtain BAF and LogR from the allele counts (Optimized)
#' @export
getBAFsAndLogRs <- function(tumourAlleleCountsFile.prefix, normalAlleleCountsFile.prefix, figuresFile.prefix, BAFnormalFile, BAFmutantFile, logRnormalFile, logRmutantFile, combinedAlleleCountsFile, chr_names, g1000file.prefix, minCounts = NA, samplename = "sample1", seed = as.integer(Sys.time())) {
  set.seed(seed)

  # Fast data loading
  input_data <- concatenateAlleleCountFiles(tumourAlleleCountsFile.prefix, ".txt", chr_names)
  normal_input_data <- concatenateAlleleCountFiles(normalAlleleCountsFile.prefix, ".txt", chr_names)
  allele_data <- concatenateG1000SnpFiles(g1000file.prefix, ".txt", chr_names)

  log_info(
    "Data Loading Complete: Tumour {nrow(input_data)} rows, Normal {nrow(normal_input_data)} rows, G1000 Ref {nrow(allele_data)} rows",
  )

  # Efficient chr prefix stripping
  allele_data[[1]] <- gsub("chr", "", allele_data[[1]])
  normal_input_data[[1]] <- gsub("chr", "", normal_input_data[[1]])
  input_data[[1]] <- gsub("chr", "", input_data[[1]])

  # Fast Synchronisation: Using match/joins is faster than Reduce(intersect(paste))
  # To maintain pixel-perfect parity with the 'paste' key logic:
  key_allele <- paste0(allele_data[[1]], "_", allele_data[[2]])
  key_normal <- paste0(normal_input_data[[1]], "_", normal_input_data[[2]])
  key_tumour <- paste0(input_data[[1]], "_", input_data[[2]])

  # Find common keys
  common_keys <- intersect(intersect(key_allele, key_normal), key_tumour)

  # Filter data frames
  allele_data <- allele_data[collapse::fmatch(common_keys, key_allele), ]
  normal_input_data <- normal_input_data[collapse::fmatch(common_keys, key_normal), ]
  input_data <- input_data[collapse::fmatch(common_keys, key_tumour), ]

  log_info("Sync complete. Remaining SNPs: {nrow(input_data)}")

  rm(key_allele, key_normal, key_tumour, common_keys)

  names(input_data)[1] <- "CHR"
  names(normal_input_data)[1] <- "CHR"
  # Using matrix indexing for fast extraction
  norm_m <- as.matrix(normal_input_data[, 3:6])
  mut_m <- as.matrix(input_data[, 3:6])

  # Map alleles to counts
  len <- nrow(norm_m)

  idx_matrix <- cbind(seq_len(len), as.integer(allele_data[[3]]))
  idx_matrix2 <- cbind(seq_len(len), as.integer(allele_data[[4]]))

  # allele_data[,3] and [,4] contain the column indices for A and B alleles
  normCount1 <- norm_m[idx_matrix]
  normCount2 <- norm_m[idx_matrix2]
  mutCount1 <- mut_m[idx_matrix]
  mutCount2 <- mut_m[idx_matrix2]

  totalNormal <- normCount1 + normCount2
  totalMutant <- mutCount1 + mutCount2

  rm(norm_m, mut_m, allele_data, normal_input_data)

  # Apply coverage filters
  indices <- seq_len(nrow(input_data))
  if (!is.na(minCounts)) {
    indices <- which(totalNormal >= minCounts & totalMutant >= 1)
    totalNormal <- totalNormal[indices]
    totalMutant <- totalMutant[indices]
    normCount1 <- normCount1[indices]
    normCount2 <- normCount2[indices]
    mutCount1 <- mutCount1[indices]
    mutCount2 <- mutCount2[indices]
  }

  n <- length(indices)

  # Allele Randomization (Pixel-Perfect logic)
  # runif(n) generates values in [0,1], round() makes them 0 or 1
  selector <- round(stats::runif(n))
  is_zero <- selector == 0
  is_one <- !is_zero

  normalBAF <- numeric(n)
  mutantBAF <- numeric(n)

  normalBAF[is_zero] <- normCount1[is_zero] / totalNormal[is_zero]
  normalBAF[is_one] <- normCount2[is_one] / totalNormal[is_one]
  mutantBAF[is_zero] <- mutCount1[is_zero] / totalMutant[is_zero]
  mutantBAF[is_one] <- mutCount2[is_one] / totalMutant[is_one]

  # LogR Calculation
  # normalLogR is forced to integer 0 as per original script requirement
  normalLogR <- integer(n)
  mutantLogR_raw <- totalMutant / totalNormal
  tumorLogR_final <- log2(mutantLogR_raw / mean(mutantLogR_raw, na.rm = TRUE))

  # Prepare shared columns
  CHR_final <- input_data[[1]][indices]
  POS_final <- input_data[[2]][indices]

  baseDT <- data.table::data.table(
    Chromosome = CHR_final,
    Position   = POS_final
  )


  # Write Normal BAF
  baseDT[, (samplename) := normalBAF]
  data.table::fwrite(baseDT, file = BAFnormalFile, sep = "\t")
  log_info("Saved Normal BAF to: {normalizePath(BAFnormalFile, mustWork = FALSE)}")

  # Write Mutant BAF
  baseDT[, (samplename) := mutantBAF]
  data.table::fwrite(baseDT, file = BAFmutantFile, sep = "\t")
  log_info("Saved Mutant BAF to: {normalizePath(BAFmutantFile, mustWork = FALSE)}")

  # Write Normal LogR
  baseDT[, (samplename) := normalLogR]
  data.table::fwrite(baseDT, file = logRnormalFile, sep = "\t")
  log_info("Saved Normal LogR to: {normalizePath(logRnormalFile, mustWork = FALSE)}")


  # Write Mutant LogR
  baseDT[, (samplename) := tumorLogR_final]
  data.table::fwrite(baseDT, file = logRmutantFile, sep = "\t")
  log_info("Saved Mutant LogR to: {normalizePath(logRmutantFile, mustWork = FALSE)}")

  # Write Combined Allele Counts
  # We use a standard data.table definition here which is safe from list-bloat
  baseDT[, (samplename) := NULL] # Clean up the sample column before combining
  combinedDT <- cbind(baseDT, data.table::data.table(
    mutCountT1 = mutCount1,
    mutCountT2 = mutCount2,
    mutCountN1 = normCount1,
    mutCountN2 = normCount2
  ))

  data.table::fwrite(combinedDT, file = combinedAlleleCountsFile, sep = "\t")
  log_info("Saved combined Allele Counts to: {normalizePath(combinedAlleleCountsFile, mustWork = FALSE)}")

  # Plotting Setup
  # Re-using vectors to build the ASCAT list object without re-reading files
  SNPpos <- data.frame(
    Chromosome = CHR_final,
    Position = POS_final,
    stringsAsFactors = FALSE
  )

  # Optimized 'ch' list creation
  ch <- lapply(chr_names, function(x) {
    tmp <- which(SNPpos$Chromosome == x)
    if (length(tmp) == 0) {
      return(0)
    }
    return(tmp[1]:tmp[length(tmp)])
  })

  ascat_bc <- list(
    Tumor_LogR = data.frame(tumorLogR_final),
    Tumor_BAF = data.frame(mutantBAF),
    Germline_LogR = data.frame(normalLogR),
    Germline_BAF = data.frame(normalBAF),
    Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL,
    Tumor_counts = NULL, Germline_counts = NULL,
    SNPpos = SNPpos,
    chrs = chr_names,
    samples = samplename,
    chrom = split_genome(SNPpos),
    ch = ch
  )

  ASCAT::ascat.plotRawData(ascat_bc)
}

#' Prepare data for impute
#'
#' @param chrom The chromosome for which impute input should be generated.
#' @param tumour_allele_counts_file Output from the allele counter on the matched tumour for this chromosome.
#' @param normal_allele_counts_file Output from the allele counter on the matched normal for this chromosome.
#' @param output_file File where the impute input for this chromosome will be written.
#' @param imputeinfofile Info file with impute reference information.
#' @param is_male Boolean denoting whether this sample is male (TRUE), or female (FALSE).
#' @param problem_loci_file A file containing genomic locations that must be discarded (optional).
#' @param use_loci_file A file containing genomic locations that must be included (optional).
#' @param heterozygous_filter The cutoff where a SNP will be considered as heterozygous (default 0.1).
#' @author dw9, sd11
#' @export
generate_impute_input_wgs <- function(
  chrom, tumour_allele_counts_file, normal_allele_counts_file,
  output_file, imputeinfofile, is_male, problem_loci_file = NA,
  use_loci_file = NA, heterozygous_filter = 0.1
) {
  # Read in the reference file paths for the specified chrom
  impute_info <- parse_imputeinfofile(imputeinfofile, is_male, chrom = chrom)
  chrom_name <- chrom

  # Efficiently load and combine known SNP legend files
  # Replaces the for-loop/rbind pattern which is very slow in R
  known_SNPs <- lapply(impute_info$impute_legend, function(file) {
    data.table::fread(file, sep = " ", header = TRUE, data.table = FALSE)
  }) |>
    data.table::rbindlist() |>
    as.data.frame()

  # Filter out 'problem' SNPs (BAF streaks)
  if (!is.na(problem_loci_file) && problem_loci_file != "NA") {
    problem_snps_raw <- data.table::fread(problem_loci_file, header = TRUE, sep = "\t", data.table = FALSE)
    problem_positions <- problem_snps_raw$Pos[problem_snps_raw$Chr == chrom_name]
    known_SNPs <- known_SNPs[!(known_SNPs$position %in% problem_positions), ]
  }

  # Filter for 'good' SNPs (e.g., SNP6 positions)
  if (!is.na(use_loci_file) && use_loci_file != "NA") {
    good_snps_raw <- data.table::fread(use_loci_file, header = TRUE, sep = "\t", data.table = FALSE)
    good_positions <- good_snps_raw$pos[good_snps_raw$chr == chrom_name]
    known_SNPs <- known_SNPs[known_SNPs$position %in% good_positions, ]
  }

  # Load allele counts using fread (ignoring comments)
  # Tumour and Normal are combined column-wise to match legacy indexing
  snp_tumour <- data.table::fread(tumour_allele_counts_file, sep = "\t", header = FALSE, data.table = FALSE)
  snp_normal <- data.table::fread(normal_allele_counts_file, sep = "\t", header = FALSE, data.table = FALSE)

  # Combined data: [Tumour Cols 1-6] [Normal Cols 7-12]
  snp_combined <- cbind(snp_tumour, snp_normal)

  # Match known SNPs to the allele counter positions
  indices <- match(known_SNPs$position, snp_combined[, 2])
  mask <- !is.na(indices)
  found_snp_data <- snp_combined[indices[mask], ]
  valid_known_snps <- known_SNPs[mask, ]

  # Calculate BAF for the NORMAL sample to determine genotypes
  # Logic: Alt / (Alt + Ref).
  # Ref column index: match allele in col 3 + normal offset (ncol) + 2
  # Alt column index: match allele in col 4 + normal offset (ncol) + 2
  nucleotides <- c("A", "C", "G", "T")
  norm_col_count <- ncol(snp_normal)

  ref_cols <- match(valid_known_snps[, 3], nucleotides) + norm_col_count + 2
  alt_cols <- match(valid_known_snps[, 4], nucleotides) + norm_col_count + 2

  # Matrix indexing for high-speed extraction of specific allele counts
  row_idx <- seq_len(nrow(found_snp_data))
  alt_counts <- as.numeric(found_snp_data[cbind(row_idx, alt_cols)])
  ref_counts <- as.numeric(found_snp_data[cbind(row_idx, ref_cols)])

  bafs <- alt_counts / (alt_counts + ref_counts)
  bafs[is.nan(bafs)] <- 0

  # Determine genotypes for IMPUTE2 (1-hot encoded: HomRef, Het, HomAlt)
  min_baf <- min(heterozygous_filter, 1.0 - heterozygous_filter)
  max_baf <- max(heterozygous_filter, 1.0 - heterozygous_filter)

  genotypes <- matrix(0, nrow = nrow(found_snp_data), ncol = 3)
  genotypes[bafs <= min_baf, 1] <- 1
  genotypes[bafs > min_baf & bafs < max_baf, 2] <- 1
  genotypes[bafs >= max_baf, 3] <- 1

  # Create final output table
  # Format: [snpID] [Chr] [Pos] [Ref] [Alt] [G1] [G2] [G3]
  snp_names <- paste0("snp", seq_len(nrow(genotypes)))
  out_data <- cbind(snp_names, valid_known_snps[, 1:4], genotypes)

  # Write main output
  data.table::fwrite(out_data, file = output_file, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Legacy check: Write sample_g.txt if chrom_name is NA (usually for non-standard chrom processing)
  if (is.na(chrom_name)) {
    sample_g_file <- file.path(dirname(output_file), "sample_g.txt")
    sample_g_data <- data.frame(
      ID_1 = c(0, "INDIVI1"),
      ID_2 = c(0, "INDIVI1"),
      missing = c(0, 0),
      sex = c("D", 2)
    )
    data.table::fwrite(sample_g_data, file = sample_g_file, sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}

#' Function to correct LogR for waivyness that correlates with GC content
#' @param Tumour_LogR_file String pointing to the tumour LogR output
#' @param outfile String pointing to where the GC corrected LogR should be written
#' @param correlations_outfile File where correlations are to be saved
#' @param gc_content_file_prefix String pointing to where GC windows for this reference genome can be
#' found. These files should be split per chromosome and this prefix must contain the full path until
#' chr in its name. The .txt extension is automatically added.
#' @param replic_timing_file_prefix Like the gc_content_file_prefix, containing replication timing info (supply NULL if no replication timing correction is to be applied)
#' @param chrom_names A vector containing chromosome names to be considered
#' @param recalc_corr_afterwards Set to TRUE to recalculate correlations after correction
#' @author jdemeul, sd11
#' @export
gc_correct_wgs <- function(
  Tumour_LogR_file,
  outfile,
  correlations_outfile,
  gc_content_file_prefix,
  replic_timing_file_prefix,
  chrom_names,
  recalc_corr_afterwards = FALSE,
  debug = FALSE
) {
  # :: syntax used
  # Pure comments instead of numbering

  if (is.null(gc_content_file_prefix)) {
    stop("GC content reference files must be supplied")
  }

  Tumor_LogR <- read_logr(Tumour_LogR_file)

  # Efficiently load and combine GC data
  gc_files <- paste0(gc_content_file_prefix, chrom_names, ".txt.gz")
  GC_data <- do.call(rbind, lapply(gc_files, read_gccontent))

  # Clean up the GC_data headers
  # The first column is often a duplicate of the third; we remove it safely
  correct_headers <- colnames(GC_data)[2:ncol(GC_data)]
  GC_data <- GC_data[, -1]
  colnames(GC_data) <- trimws(correct_headers)
  data.table::setnames(GC_data, old = 1:2, new = c("Chromosome", "Position"))

  # Processing replication data if prefix is provided
  has_replic <- !is.null(replic_timing_file_prefix)
  if (has_replic) {
    replic_files <- paste0(replic_timing_file_prefix, chrom_names, ".txt.gz")
    replic_data <- do.call(rbind, lapply(replic_files, read_replication))
    colnames(replic_data) <- trimws(colnames(replic_data))
    if ("pos" %in% colnames(replic_data)) data.table::setnames(replic_data, "pos", "Position")
    if ("chr" %in% colnames(replic_data)) data.table::setnames(replic_data, "chr", "Chromosome")
  }

  # Fast Loci Matching
  logr_key <- paste0(Tumor_LogR$Chromosome, "_", Tumor_LogR$Position)
  gc_key <- paste0(GC_data$Chromosome, "_", GC_data$Position)
  locimatches <- match(logr_key, gc_key)

  num_matches <- sum(!is.na(locimatches))
  log_info("Alignment check: {num_matches} / {nrow(Tumor_LogR)} positions matched.")

  if (num_matches == 0) {
    log_failure("Zero overlap found! Check if LogR is hg19 while GC refs are hg38.")
  }


  valid_idx <- which(!is.na(locimatches))
  matched_gc <- locimatches[valid_idx]

  # Subsetting objects to matched rows
  Tumor_LogR <- Tumor_LogR[valid_idx, ]
  GC_data <- GC_data[matched_gc, ]
  if (has_replic) replic_data <- replic_data[matched_gc, ]

  # Clean up memory
  rm(logr_key, gc_key, locimatches)

  # Calculate correlations and identify best window sizes
  # We use collapse::pwcor for speed
  corr <- collapse::pwcor(GC_data[, 3:ncol(GC_data)], Tumor_LogR[[3]], use = "pairwise.complete.obs")
  corr <- abs(corr[, 1])

  # instead of capping it at 100kb go to the end of the frame
  index_2kb <- which(names(corr) == "2kb")
  maxGCcol_insert <- names(which.max(corr[1:index_2kb]))
  maxGCcol_amplic <- names(which.max(corr[(index_2kb + 1):length(corr)]))
  index_100kb <- which(names(corr) == "100kb")
  maxGCcol_amplic <- names(which.max(corr[(index_2kb + 2):index_100kb]))

  # Construct the design matrix for splines
  # We use intercept = TRUE for the first and FALSE for the others to avoid rank deficiency
  if (has_replic) {
    corr_rep <- collapse::pwcor(replic_data[, 3:ncol(replic_data)], Tumor_LogR[[3]], use = "pairwise.complete.obs")
    corr_rep <- abs(corr_rep[, 1])
    maxreplic <- names(which.max(corr_rep))

    X <- cbind(
      splines::ns(GC_data[[maxGCcol_insert]], df = 5, intercept = TRUE),
      splines::ns(GC_data[[maxGCcol_amplic]], df = 5, intercept = FALSE),
      splines::ns(replic_data[[maxreplic]], df = 5, intercept = FALSE)
    )
  } else {
    X <- cbind(
      splines::ns(GC_data[[maxGCcol_insert]], df = 5, intercept = TRUE),
      splines::ns(GC_data[[maxGCcol_amplic]], df = 5, intercept = FALSE)
    )
  }

  y <- as.numeric(Tumor_LogR[[3]])

  # Robust Linear Model fitting
  # We use stats::lm.fit directly for a balance of speed and numerical stability
  # It is faster than lm() but more stable than flm() for splines
  keep_idx <- stats::complete.cases(X) & !is.na(y)
  fit <- stats::lm.fit(x = as.matrix(X[keep_idx, ]), y = y[keep_idx])

  # Calculate residuals and cap them to remove outliers
  resids <- rep(NA, length(y))
  resids[keep_idx] <- fit$residuals
  resids <- pmax(pmin(resids, 5), -5)

  # Metrics for noise reduction
  sd_before <- stats::sd(y, na.rm = TRUE)
  sd_after <- stats::sd(resids, na.rm = TRUE)
  reduction <- ((sd_before - sd_after) / sd_before) * 100

  # Apply corrected LogR
  Tumor_LogR[[3]] <- resids

  # Log results
  message(paste0("Noise Reduction: ", round(reduction, 2), "%"))

  sd_before <- stats::sd(y, na.rm = TRUE)
  sd_after <- stats::sd(resids, na.rm = TRUE)
  reduction <- ((sd_before - sd_after) / sd_before) * 100
  # Post-correction correlation check
  corr_post_short <- abs(stats::cor(resids[keep_idx], GC_data[[maxGCcol_insert]][keep_idx], use = "complete.obs"))
  corr_post_long <- abs(stats::cor(resids[keep_idx], GC_data[[maxGCcol_amplic]][keep_idx], use = "complete.obs"))

  # Glue Log: Interpretation block
  log_info("Noise Reduction (SD): {round(reduction, 2)}%")
  log_info("Residual Correlation (Short): {round(corr_post_short, 4)} (Target: ~0)")
  log_info("Residual Correlation (Long): {round(corr_post_long, 4)} (Target: ~0)")
  log_info("LogR Mean Shift: {round(mean(resids, na.rm=TRUE), 6)} (Target: 0)")

  # Write corrected LogR
  data.table::fwrite(
    x = Tumor_LogR[!is.na(Tumor_LogR[[3]]), ],
    file = outfile,
    sep = "\t",
    quote = FALSE
  )
}

#' Prepare WGS data for haplotype construction
#'
#' This function performs part of the Battenberg WGS pipeline: Counting alleles, constructing BAF and logR
#' and performing GC content correction.
#'
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param tumourbam Full path to the tumour BAM file
#' @param normalbam Full path to the normal BAM file
#' @param tumourname Identifier to be used for tumour output files
#' @param normalname Identifier to be used for normal output files
#' @param g1000allelesprefix Prefix path to the 1000 Genomes alleles reference files
#' @param g1000prefix Prefix path to the 1000 Genomes SNP reference files
#' @param gccorrectprefix Prefix path to GC content reference data
#' @param repliccorrectprefix Prefix path to replication timing reference data (supply NULL if no replication timing correction is to be applied)
#' @param min_base_qual Minimum base quality required for a read to be counted
#' @param min_map_qual Minimum mapping quality required for a read to be counted
#' @param allelecounter_exe Path to the allele counter executable (can be found in $PATH)
#' @param min_normal_depth Minimum depth required in the normal for a SNP to be included
#' @param nthreads The number of paralel processes to run
#' @param skip_allele_counting Flag, set to TRUE if allele counting is already complete (files are expected in the working directory on disk)
#' @param skip_allele_counting_normal Flag, set to TRUE from the second sample onwards for multisample case (Default: FALSE)
#' @author sd11
#' @export
prepare_wgs <- function(
  chrom_names,
  tumourbam,
  normalbam,
  tumourname,
  normalname,
  g1000allelesprefix,
  g1000prefix,
  gccorrectprefix,
  repliccorrectprefix,
  min_base_qual,
  min_map_qual,
  allelecounter_exe,
  min_normal_depth,
  nthreads,
  skip_allele_counting,
  skip_allele_counting_normal = FALSE,
  debug = FALSE
) {
  if (!skip_allele_counting) {
    do_allele_counting <- function(i) {
      getAlleleCounts(
        bam.file = normalbam,
        output_file = paste(
          normalname,
          "_alleleFrequencies_chr",
          chrom_names[i], ".txt",
          sep = ""
        ),
        g1000.loci = paste(g1000prefix, chrom_names[i], ".txt", sep = ""),
        min.base.qual = min_base_qual,
        min.map.qual = min_map_qual,
        allelecounter.exe = allelecounter_exe
      )

      if (!skip_allele_counting_normal) {
        getAlleleCounts(
          bam.file = normalbam,
          output_file = paste(normalname,
            "_alleleFrequencies_chr",
            chrom_names[i], ".txt",
            sep = ""
          ),
          g1000.loci = paste(
            g1000prefix,
            chrom_names[i], ".txt",
            sep = ""
          ),
          min.base.qual = min_base_qual,
          min.map.qual = min_map_qual,
          allelecounter.exe = allelecounter_exe
        )
      }
    }
    run_parallel_or_serial(
      iterator = seq_along(chrom_names),
      func = do_allele_counting,
      debug = debug
    )
  }

  # Obtain BAF and LogR from the raw allele counts
  # getBAFsAndLogRs(
  #  tumourAlleleCountsFile.prefix = paste(tumourname, "_alleleFrequencies_chr", sep = ""),
  #  normalAlleleCountsFile.prefix = paste(normalname, "_alleleFrequencies_chr", sep = ""),
  #  figuresFile.prefix = paste(tumourname, "_", sep = ""),
  #  BAFnormalFile = paste(tumourname, "_normalBAF.tab", sep = ""),
  #  BAFmutantFile = paste(tumourname, "_mutantBAF.tab", sep = ""),
  #  logRnormalFile = paste(tumourname, "_normalLogR.tab", sep = ""),
  #  logRmutantFile = paste(tumourname, "_mutantLogR.tab", sep = ""),
  #  combinedAlleleCountsFile = paste(tumourname, "_alleleCounts.tab", sep = ""),
  #  chr_names = chrom_names,
  #  g1000file.prefix = g1000allelesprefix,
  #   minCounts = min_normal_depth,
  #   samplename = tumourname
  # )
  # Perform GC correction
  # gc_correct_wgs(
  #  Tumour_LogR_file = paste(tumourname, "_mutantLogR.tab", sep = ""),
  #  outfile = paste(tumourname, "_mutantLogR_gcCorrected.tab", sep = ""),
  #  correlations_outfile = paste(tumourname, "_GCwindowCorrelations.txt", sep = ""),
  #  gc_content_file_prefix = gccorrectprefix,
  #  replic_timing_file_prefix = repliccorrectprefix,
  #  chrom_names = chrom_names
  # )
}

#' A helper function to split the genome into parts
#' @param SNPpos A data.frame with a row for each SNP. First column is chromosome, second column position
#' @noRd
split_genome <- function(SNPpos) {
  # look for gaps of more than 1Mb and chromosome borders
  holesOver1Mb <- which(diff(SNPpos[, 2]) >= 1000000) + 1
  chrBorders <- which(diff(as.numeric(factor(SNPpos[, 1], levels = unique(SNPpos[, 1])))) != 0) + 1
  holes <- unique(sort(c(holesOver1Mb, chrBorders)))

  # find which segments are too small
  joincandidates <- which(diff(c(0, holes, dim(SNPpos)[1])) < 200)

  # if it's the first or last segment, just join to the one next to it, irrespective of chromosome and positions
  while (1 %in% joincandidates) {
    holes <- holes[-1]
    joincandidates <- which(diff(c(0, holes, dim(SNPpos)[1])) < 200)
  }
  while ((length(holes) + 1) %in% joincandidates) {
    holes <- holes[-length(holes)]
    joincandidates <- which(diff(c(0, holes, dim(SNPpos)[1])) < 200)
  }

  while (length(joincandidates) != 0) {
    # the while loop is because after joining, segments may still be too small..
    startseg <- c(1, holes)
    endseg <- c(holes - 1, dim(SNPpos)[1])

    # for each segment that is too short, see if it has the same chromosome as the segments before and after
    # the next always works because neither the first or the last segment is in joincandidates now
    previoussamechr <- SNPpos[endseg[joincandidates - 1], 1] == SNPpos[startseg[joincandidates], 1]
    nextsamechr <- SNPpos[endseg[joincandidates], 1] == SNPpos[startseg[joincandidates + 1], 1]

    distanceprevious <- SNPpos[startseg[joincandidates], 2] - SNPpos[endseg[joincandidates - 1], 2]
    distancenext <- SNPpos[startseg[joincandidates + 1], 2] - SNPpos[endseg[joincandidates], 2]

    # if both the same, decide based on distance, otherwise if one the same, take the other, if none, just take one.
    joins <- ifelse(previoussamechr & nextsamechr,
      ifelse(distanceprevious > distancenext, joincandidates, joincandidates - 1),
      ifelse(nextsamechr, joincandidates, joincandidates - 1)
    )

    holes <- holes[-joins]
    joincandidates <- which(diff(c(0, holes, dim(SNPpos)[1])) < 200)
  }
  # if two neighboring segments are selected, this may make bigger segments then absolutely necessary.
  startseg <- c(1, holes)
  endseg <- c(holes - 1, dim(SNPpos)[1])
  chr <- list()
  for (i in seq_along(startseg)) {
    chr[[i]] <- startseg[i]:endseg[i]
  }

  return(chr)
}
