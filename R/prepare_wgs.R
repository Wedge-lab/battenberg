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

  rm(key_allele, key_normal, key_tumour, common_keys)

  # Map alleles to counts
  len <- nrow(normal_input_data)
  # Using matrix indexing for fast extraction
  norm_m <- as.matrix(normal_input_data[, 3:6])
  mut_m <- as.matrix(input_data[, 3:6])

  # allele_data[,3] and [,4] contain the column indices for A and B alleles
  normCount1 <- norm_m[cbind(seq_len(len), allele_data[[3]])]
  normCount2 <- norm_m[cbind(seq_len(len), allele_data[[4]])]
  mutCount1 <- mut_m[cbind(seq_len(len), allele_data[[3]])]
  mutCount2 <- mut_m[cbind(seq_len(len), allele_data[[4]])]

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
  selector <- round(runif(n))
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

  # Fast File Saving (Direct List writing avoids data.frame overhead)
  data.table::fwrite(
    list(CHR_final, POS_final, normalBAF),
    file = BAFnormalFile,
    sep = "\t", col.names = c("Chromosome", "Position", samplename)
  )
  data.table::fwrite(
    list(CHR_final, POS_final, mutantBAF),
    file = BAFmutantFile,
    sep = "\t", col.names = c("Chromosome", "Position", samplename)
  )
  data.table::fwrite(
    list(CHR_final, POS_final, normalLogR),
    file = logRnormalFile,
    sep = "\t", col.names = c("Chromosome", "Position", samplename)
  )
  data.table::fwrite(
    list(CHR_final, POS_final, tumorLogR_final),
    file = logRmutantFile, sep = "\t",
    col.names = c("Chromosome", "Position", samplename)
  )
  data.table::fwrite(
    list(CHR_final, POS_final, mutCount1, mutCount2, normCount1, normCount2),
    file = combinedAlleleCountsFile, sep = "\t",
    col.names = c("Chromosome", "Position", "mutCountT1", "mutCountT2", "mutCountN1", "mutCountN2")
  )

  # Plotting Setup
  # Re-using vectors to build the ASCAT list object without re-reading files
  SNPpos <- data.frame(Chromosome = CHR_final, Position = POS_final, stringsAsFactors = FALSE)

  # Optimized 'ch' list creation
  ch <- lapply(chr_names, function(x) {
    tmp <- which(SNPpos$Chromosome == x)
    if (length(tmp) == 0) {
      return(0)
    }
    return(tmp[1]:tmp[length(tmp)])
  })

  ascat.bc <- list(
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

  ASCAT::ascat.plotRawData(ascat.bc)
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
  data.table::fwrite(out_data, file = output_file, sep = " ", row.names = FALSE, col_names = FALSE, quote = FALSE)

  # Legacy check: Write sample_g.txt if chrom_name is NA (usually for non-standard chrom processing)
  if (is.na(chrom_name)) {
    sample_g_file <- file.path(dirname(output_file), "sample_g.txt")
    sample_g_data <- data.frame(
      ID_1 = c(0, "INDIVI1"),
      ID_2 = c(0, "INDIVI1"),
      missing = c(0, 0),
      sex = c("D", 2)
    )
    data.table::fwrite(sample_g_data, file = sample_g_file, sep = " ", row.names = FALSE, col_names = TRUE, quote = FALSE)
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
  recalc_corr_afterwards = FALSE
) {
  if (is.null(gc_content_file_prefix)) {
    stop("GC content reference files must be supplied to WGS GC content correction")
  }

  Tumor_LogR <- read_logr(Tumour_LogR_file)

  # Processing GC data
  print("Processing GC content data")
  gc_files <- paste0(gc_content_file_prefix, chrom_names, ".txt.gz")
  GC_data <- do.call(rbind, lapply(gc_files, read_gccontent))
  colnames(GC_data) <- c(
    "chr", "Position", paste0(c(25, 50, 100, 200, 500), "bp"),
    paste0(c(1, 2, 5, 10, 20, 50, 100), "kb")
  )

  # Processing replication data
  has_replic <- !is.null(replic_timing_file_prefix)
  if (has_replic) {
    print("Processing replication timing data")
    replic_files <- paste0(replic_timing_file_prefix, chrom_names, ".txt.gz")
    replic_data <- do.call(rbind, lapply(replic_files, read_replication))
  }

  # Matching loci - using a more efficient matching key
  # Pixel-perfect match to: paste0(Tumor_LogR$Chromosome, "_", Tumor_LogR$Position)
  logr_key <- paste0(Tumor_LogR$Chromosome, "_", Tumor_LogR$Position)
  gc_key <- paste0(GC_data$chr, "_", GC_data$Position)
  locimatches <- match(logr_key, gc_key)

  valid_idx <- which(!is.na(locimatches))
  matched_gc <- locimatches[valid_idx]

  Tumor_LogR <- Tumor_LogR[valid_idx, ]
  GC_data <- GC_data[matched_gc, ]

  if (has_replic) {
    replic_data <- replic_data[matched_gc, ]
  }
  rm(logr_key, gc_key, locimatches, valid_idx, matched_gc)

  # Initial Correlations
  corr <- abs(collapse::fcor(GC_data[, 3:ncol(GC_data)], Tumor_LogR[, 3], use = "complete.obs")[, 1])
  if (has_replic) {
    corr_rep <- abs(collapse::fcor(replic_data[, 3:ncol(replic_data)], Tumor_LogR[, 3], use = "complete.obs")[, 1])
  }

  # Identify best windows
  index_1kb <- which(names(corr) == "1kb")
  maxGCcol_insert <- names(which.max(corr[1:index_1kb]))
  index_100kb <- which(names(corr) == "100kb")
  maxGCcol_amplic <- names(which.max(corr[(index_1kb + 2):index_100kb]))

  if (has_replic) {
    maxreplic <- names(which.max(corr_rep))
    cat("Replication timing correlation: ", paste(names(corr_rep), format(corr_rep, digits = 2), ";"), "\n")
    cat("Replication dataset: ", maxreplic, "\n")
  }
  cat("GC correlation: ", paste(names(corr), format(corr, digits = 2), ";"), "\n")
  cat("Short window size: ", maxGCcol_insert, "\n")
  cat("Long window size: ", maxGCcol_amplic, "\n")

  # Write 'before' correlations
  corr_df_save <- if (has_replic) {
    data.frame(windowsize = c(names(corr), names(corr_rep)), correlation = c(corr, corr_rep))
  } else {
    data.frame(windowsize = names(corr), correlation = corr)
  }
  data.table::fwrite(corr_df_save, file = gsub(".txt", "_beforeCorrection.txt", correlations_outfile), sep = "\t")

  # Setup Design Matrix (X) for Linear Model
  # This replaces the lm() formula interface
  if (has_replic) {
    X <- stats::model.matrix(~ splines::ns(GC_data[[maxGCcol_insert]], df = 5, intercept = TRUE) +
      splines::ns(GC_data[[maxGCcol_amplic]], df = 5, intercept = TRUE) +
      splines::ns(replic_data[[maxreplic]], df = 5, intercept = TRUE))
  } else {
    X <- stats::model.matrix(~ splines::ns(GC_data[[maxGCcol_insert]], df = 5, intercept = TRUE) +
      splines::ns(GC_data[[maxGCcol_amplic]], df = 5, intercept = TRUE))
  }

  # Pixel-perfect NA handling (na.exclude behavior)
  y <- Tumor_LogR[, 3, drop = TRUE]
  keep_idx <- stats::complete.cases(X) & !is.na(y)

  # Solve OLS using fast C++ backend
  y_clean <- y[keep_idx]
  X_clean <- X[keep_idx, , drop = FALSE]
  betas <- collapse::flm(y_clean, X_clean)

  # Reconstruct residuals (Observed - Predicted)
  # Pre-filling with NA matches 'na.exclude' padding
  resids <- rep(NA, length(y))
  resids[keep_idx] <- y_clean - as.vector(X_clean %*% betas)

  # Update LogR and clean up predictors if requested
  Tumor_LogR[, 3] <- resids

  if (!recalc_corr_afterwards) {
    rm(GC_data)
    if (has_replic) rm(replic_data)
  }
  rm(X, X_clean, y_clean, betas, resids)

  # Write corrected LogR
  readr::write_tsv(x = Tumor_LogR[!is.na(Tumor_LogR[, 3]), ], file = outfile)

  # Post-correction processing
  if (recalc_corr_afterwards) {
    corr_post <- abs(collapse::fcor(GC_data[, 3:ncol(GC_data)], Tumor_LogR[, 3], use = "complete.obs")[, 1])
    if (has_replic) {
      corr_rep_post <- abs(collapse::fcor(replic_data[, 3:ncol(replic_data)], Tumor_LogR[, 3], use = "complete.obs")[, 1])
      cat("Replication timing correlation post correction: ", paste(names(corr_rep_post), format(corr_rep_post, digits = 2), ";"), "\n")

      corr_final <- data.frame(
        windowsize = c(names(corr_post), names(corr_rep_post)),
        correlation = c(corr_post, corr_rep_post)
      )
    } else {
      cat("GC correlation post correction: ", paste(names(corr_post), format(corr_post, digits = 2), ";"), "\n")
      corr_final <- data.frame(windowsize = names(corr_post), correlation = corr_post)
    }
    data.table::fwrite(corr_final, file = gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep = "\t")
  } else {
    # If not recalculating, set correlation to NA as per original code
    corr_df_save$correlation <- NA
    data.table::fwrite(corr_df_save, file = gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep = "\t")
  }
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
  skip_allele_counting_normal = FALSE
) {
  `%dopar%` <- foreach::`%dopar%`
  if (!skip_allele_counting) {
    # Obtain allele counts for 1000 Genomes locations for both tumour and normal
    foreach::foreach(i = seq_along(chrom_names)) %dopar% {
      getAlleleCounts(
        bam.file = tumourbam,
        output_file = paste(tumourname, "_alleleFrequencies_chr", chrom_names[i], ".txt", sep = ""),
        g1000.loci = paste(g1000prefix, chrom_names[i], ".txt", sep = ""),
        min.base.qual = min_base_qual,
        min.map.qual = min_map_qual,
        allelecounter.exe = allelecounter_exe
      )

      if (!skip_allele_counting_normal) {
        getAlleleCounts(
          bam.file = normalbam,
          output_file = paste(normalname, "_alleleFrequencies_chr", chrom_names[i], ".txt", sep = ""),
          g1000.loci = paste(g1000prefix, chrom_names[i], ".txt", sep = ""),
          min.base.qual = min_base_qual,
          min.map.qual = min_map_qual,
          allelecounter.exe = allelecounter_exe
        )
      }
    }
  }

  # Obtain BAF and LogR from the raw allele counts
  getBAFsAndLogRs(
    tumourAlleleCountsFile.prefix = paste(tumourname, "_alleleFrequencies_chr", sep = ""),
    normalAlleleCountsFile.prefix = paste(normalname, "_alleleFrequencies_chr", sep = ""),
    figuresFile.prefix = paste(tumourname, "_", sep = ""),
    BAFnormalFile = paste(tumourname, "_normalBAF.tab", sep = ""),
    BAFmutantFile = paste(tumourname, "_mutantBAF.tab", sep = ""),
    logRnormalFile = paste(tumourname, "_normalLogR.tab", sep = ""),
    logRmutantFile = paste(tumourname, "_mutantLogR.tab", sep = ""),
    combinedAlleleCountsFile = paste(tumourname, "_alleleCounts.tab", sep = ""),
    chr_names = chrom_names,
    g1000file.prefix = g1000allelesprefix,
    minCounts = min_normal_depth,
    samplename = tumourname
  )
  # Perform GC correction
  gc_correct_wgs(
    Tumour_LogR_file = paste(tumourname, "_mutantLogR.tab", sep = ""),
    outfile = paste(tumourname, "_mutantLogR_gcCorrected.tab", sep = ""),
    correlations_outfile = paste(tumourname, "_GCwindowCorrelations.txt", sep = ""),
    gc_content_file_prefix = gccorrectprefix,
    replic_timing_file_prefix = repliccorrectprefix,
    chrom_names = chrom_names
  )
}
