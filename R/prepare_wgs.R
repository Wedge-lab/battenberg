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


#' Obtain BAF and LogR from the allele counts
#'
#' @param tumourAlleleCountsFile.prefix Prefix of the allele counts files for the tumour.
#' @param normalAlleleCountsFile.prefix Prefix of the allele counts files for the normal.
#' @param figuresFile.prefix Prefix for output figures file names.
#' @param BAFnormalFile File where BAF from the normal will be written.
#' @param BAFmutantFile File where BAF from the tumour will be written.
#' @param logRnormalFile File where LogR from the normal will be written.
#' @param logRmutantFile File where LogR from the tumour will be written.
#' @param combinedAlleleCountsFile File where combined allele counts for tumour and normal will be written.
#' @param chr_names A vector with allowed chromosome names.
#' @param g1000file.prefix Prefix to where 1000 Genomes reference files can be found.
#' @param minCounts Integer, minimum depth required for a SNP to be included (optional, default=NA).
#' @param samplename String, name of the sample (optional, default=sample1).
#' @param seed A seed to be set for when randomising the alleles.
#' @author dw9, sd11
#' @export
getBAFsAndLogRs <- function(tumourAlleleCountsFile.prefix, normalAlleleCountsFile.prefix, figuresFile.prefix, BAFnormalFile, BAFmutantFile, logRnormalFile, logRmutantFile, combinedAlleleCountsFile, chr_names, g1000file.prefix, minCounts = NA, samplename = "sample1", seed = as.integer(Sys.time())) {
  set.seed(seed)

  input_data <- concatenateAlleleCountFiles(tumourAlleleCountsFile.prefix, ".txt", chr_names)
  normal_input_data <- concatenateAlleleCountFiles(normalAlleleCountsFile.prefix, ".txt", chr_names)
  allele_data <- concatenateG1000SnpFiles(g1000file.prefix, ".txt", chr_names)

  # We're no longer stripping out the "chr", which is causing problems
  allele_data[, 1] <- gsub("chr", "", allele_data[, 1])
  normal_input_data[, 1] <- gsub("chr", "", normal_input_data[, 1])
  input_data[, 1] <- gsub("chr", "", input_data[, 1])

  # Synchronise all the data frames
  chrpos_allele <- paste(allele_data[, 1], "_", allele_data[, 2], sep = "")
  chrpos_normal <- paste(normal_input_data[, 1], "_", normal_input_data[, 2], sep = "")
  chrpos_tumour <- paste(input_data[, 1], "_", input_data[, 2], sep = "")
  matched_data <- Reduce(intersect, list(chrpos_allele, chrpos_normal, chrpos_tumour))

  allele_data <- allele_data[chrpos_allele %in% matched_data, ]
  normal_input_data <- normal_input_data[chrpos_normal %in% matched_data, ]
  input_data <- input_data[chrpos_tumour %in% matched_data, ]

  # Clean up and reduce amount of unneeded data
  names(input_data)[1] <- "CHR"
  names(normal_input_data)[1] <- "CHR"

  normal_data <- normal_input_data[, 3:6]
  mutant_data <- input_data[, 3:6]

  # Obtain depth for both alleles for tumour and normal
  len <- nrow(normal_data)
  normCount1 <- normal_data[cbind(1:len, allele_data[, 3])]
  normCount2 <- normal_data[cbind(1:len, allele_data[, 4])]
  totalNormal <- normCount1 + normCount2
  mutCount1 <- mutant_data[cbind(1:len, allele_data[, 3])]
  mutCount2 <- mutant_data[cbind(1:len, allele_data[, 4])]
  totalMutant <- mutCount1 + mutCount2

  # Clean up a few unused variables to save some memory
  rm(normal_data, mutant_data, allele_data, normal_input_data)

  # Clear SNPs where there is not enough coverage
  indices <- seq_len(nrow(input_data))
  if (!is.na(minCounts)) {
    print(paste("minCount=", minCounts, sep = ""))
    # Only normal has to have min coverage, mutant must have at least 1 read to prevent division by zero
    indices <- which(totalNormal >= minCounts & totalMutant >= 1)

    totalNormal <- totalNormal[indices]
    totalMutant <- totalMutant[indices]
    normCount1 <- normCount1[indices]
    normCount2 <- normCount2[indices]
    mutCount1 <- mutCount1[indices]
    mutCount2 <- mutCount2[indices]
  }
  n <- length(indices)

  normalBAF <- vector(length = n, mode = "numeric")
  mutantBAF <- vector(length = n, mode = "numeric")
  normalLogR <- vector(length = n, mode = "numeric")
  mutantLogR <- vector(length = n, mode = "numeric")

  # randomise A and B alleles
  selector <- round(runif(n))
  normalBAF[which(selector == 0)] <- normCount1[which(selector == 0)] / totalNormal[which(selector == 0)]
  normalBAF[which(selector == 1)] <- normCount2[which(selector == 1)] / totalNormal[which(selector == 1)]
  mutantBAF[which(selector == 0)] <- mutCount1[which(selector == 0)] / totalMutant[which(selector == 0)]
  mutantBAF[which(selector == 1)] <- mutCount2[which(selector == 1)] / totalMutant[which(selector == 1)]

  normalLogR <- vector(length = n, mode = "integer") # assume that normallogR is 0, and normalise mutantLogR to normalLogR
  mutantLogR <- totalMutant / totalNormal
  rm(selector)

  # Create the output data.frames
  germline.BAF <- data.frame(Chromosome = input_data$CHR[indices], Position = input_data$POS[indices], baf = normalBAF)
  germline.LogR <- data.frame(Chromosome = input_data$CHR[indices], Position = input_data$POS[indices], samplename = normalLogR)
  tumor.BAF <- data.frame(Chromosome = input_data$CHR[indices], Position = input_data$POS[indices], baf = mutantBAF)
  tumor.LogR <- data.frame(Chromosome = input_data$CHR[indices], Position = input_data$POS[indices], samplename = log2(mutantLogR / mean(mutantLogR, na.rm = TRUE)))
  alleleCounts <- data.frame(Chromosome = input_data$CHR[indices], Position = input_data$POS[indices], mutCountT1 = mutCount1, mutCountT2 = mutCount2, mutCountN1 = normCount1, mutCountN2 = normCount2)

  # Save data.frames to disk
  data.table::fwrite(germline.BAF, file = BAFnormalFile, row.names = FALSE, quote = FALSE, sep = "\t", col_names = c("Chromosome", "Position", samplename))
  data.table::fwrite(tumor.BAF, file = BAFmutantFile, row.names = FALSE, quote = FALSE, sep = "\t", col_names = c("Chromosome", "Position", samplename))
  data.table::fwrite(germline.LogR, file = logRnormalFile, row.names = FALSE, quote = FALSE, sep = "\t", col_names = c("Chromosome", "Position", samplename))
  data.table::fwrite(tumor.LogR, file = logRmutantFile, row.names = FALSE, quote = FALSE, sep = "\t", col_names = c("Chromosome", "Position", samplename))
  data.table::fwrite(alleleCounts, file = combinedAlleleCountsFile, row.names = FALSE, quote = FALSE, sep = "\t")

  # Plot the raw data using ASCAT
  # Manually create an ASCAT object, which saves reading in the above files again
  SNPpos <- germline.BAF[, c("Chromosome", "Position")]
  ch <- list()
  for (i in seq_along(chr_names)) {
    temp <- which(SNPpos$Chromosome == chr_names[i])
    if (length(temp) == 0) {
      ch[[i]] <- 0
    } else {
      ch[[i]] <- temp[1]:temp[length(temp)]
    }
  }

  ascat.bc <- list(
    Tumor_LogR = as.data.frame(tumor.LogR[, 3]), Tumor_BAF = as.data.frame(tumor.BAF[, 3]),
    Germline_LogR = as.data.frame(germline.LogR[, 3]), Germline_BAF = as.data.frame(germline.BAF[, 3]),
    Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL, Tumor_counts = NULL, Germline_counts = NULL,
    SNPpos = tumor.LogR[, 1:2], chrs = chr_names, samples = c(samplename), chrom = split_genome(tumor.LogR[, 1:2]),
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

  print("Processing GC content data")
  gc_files <- paste0(gc_content_file_prefix, chrom_names, ".txt.gz")
  GC_data <- do.call(rbind, lapply(gc_files, read_gccontent))
  colnames(GC_data) <- c(
    "chr", "Position", paste0(c(25, 50, 100, 200, 500), "bp"),
    paste0(c(1, 2, 5, 10, 20, 50, 100), "kb")
  )

  if (!is.null(replic_timing_file_prefix)) {
    print("Processing replication timing data")
    replic_files <- paste0(replic_timing_file_prefix, chrom_names, ".txt.gz")
    replic_data <- do.call(rbind, lapply(replic_files, read_replication))
  }

  # omit non-matching loci, replication data generated at exactly same GC loci
  locimatches <- match(
    x = paste0(Tumor_LogR$Chromosome, "_", Tumor_LogR$Position),
    table = paste0(GC_data$chr, "_", GC_data$Position)
  )
  Tumor_LogR <- Tumor_LogR[which(!is.na(locimatches)), ]
  GC_data <- GC_data[na.omit(locimatches), ]
  if (!is.null(replic_timing_file_prefix)) {
    replic_data <- replic_data[na.omit(locimatches), ]
  }
  rm(locimatches)

  corr <- abs(cor(GC_data[, 3:ncol(GC_data)], Tumor_LogR[, 3], use = "complete.obs")[, 1])
  if (!is.null(replic_timing_file_prefix)) {
    corr_rep <- abs(cor(replic_data[, 3:ncol(replic_data)], Tumor_LogR[, 3], use = "complete.obs")[, 1])
  }

  index_1kb <- which(names(corr) == "1kb")
  maxGCcol_insert <- names(which.max(corr[1:index_1kb]))
  index_100kb <- which(names(corr) == "100kb")
  # start large window sizes at 5kb rather than 2kb to avoid overly correlated expl variables
  maxGCcol_amplic <- names(which.max(corr[(index_1kb + 2):index_100kb]))
  if (!is.null(replic_timing_file_prefix)) {
    maxreplic <- names(which.max(corr_rep))
  }

  if (!is.null(replic_timing_file_prefix)) {
    cat("Replication timing correlation: ", paste(names(corr_rep), format(corr_rep, digits = 2), ";"), "\n")
    cat("Replication dataset: ", maxreplic, "\n")
  }
  cat("GC correlation: ", paste(names(corr), format(corr, digits = 2), ";"), "\n")
  cat("Short window size: ", maxGCcol_insert, "\n")
  cat("Long window size: ", maxGCcol_amplic, "\n")

  if (!is.null(replic_timing_file_prefix)) {
    # Multiple regression - with replication timing
    corrdata <- data.frame(
      logr = Tumor_LogR[, 3, drop = TRUE],
      GC_insert = GC_data[, maxGCcol_insert, drop = TRUE],
      GC_amplic = GC_data[, maxGCcol_amplic, drop = TRUE],
      replic = replic_data[, maxreplic, drop = TRUE]
    )
    colnames(corrdata) <- c("logr", "GC_insert", "GC_amplic", "replic")
    if (!recalc_corr_afterwards) {
      rm(GC_data, replic_data)
    }

    model <- lm(logr ~ splines::ns(x = GC_insert, df = 5, intercept = TRUE) + splines::ns(x = GC_amplic, df = 5, intercept = TRUE) + splines::ns(x = replic, df = 5, intercept = TRUE), y = FALSE, model = FALSE, data = corrdata, na.action = "na.exclude")

    corr <- data.frame(windowsize = c(names(corr), names(corr_rep)), correlation = c(corr, corr_rep))
    data.table::fwrite(corr, file = gsub(".txt", "_beforeCorrection.txt", correlations_outfile), sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    # Multiple regression  - without replication timing
    corrdata <- data.frame(
      logr = Tumor_LogR[, 3, drop = TRUE],
      GC_insert = GC_data[, maxGCcol_insert, drop = TRUE],
      GC_amplic = GC_data[, maxGCcol_amplic, drop = TRUE]
    )
    colnames(corrdata) <- c("logr", "GC_insert", "GC_amplic")
    if (!recalc_corr_afterwards) {
      rm(GC_data)
    }

    model <- lm(logr ~ splines::ns(x = GC_insert, df = 5, intercept = TRUE) + splines::ns(x = GC_amplic, df = 5, intercept = TRUE), y = FALSE, model = FALSE, data = corrdata, na.action = "na.exclude")

    corr <- data.frame(windowsize = names(corr), correlation = corr)
    data.table::fwrite(corr, file = gsub(".txt", "_beforeCorrection.txt", correlations_outfile), sep = "\t", quote = FALSE, row.names = FALSE)
  }

  Tumor_LogR[, 3] <- residuals(model)
  rm(model, corrdata)

  readr::write_tsv(x = Tumor_LogR[which(!is.na(Tumor_LogR[, 3])), ], file = outfile)

  if (recalc_corr_afterwards) {
    # Recalculate the correlations to see how much there is left
    corr <- abs(cor(GC_data[, 3:ncol(GC_data)], Tumor_LogR[, 3], use = "complete.obs")[, 1])
    if (!is.null(replic_timing_file_prefix)) {
      corr_rep <- abs(cor(replic_data[, 3:ncol(replic_data)], Tumor_LogR[, 3], use = "complete.obs")[, 1])
      cat("Replication timing correlation post correction: ", paste(names(corr_rep), format(corr_rep, digits = 2), ";"), "\n")
    }
    cat("GC correlation post correction: ", paste(names(corr), format(corr, digits = 2), ";"), "\n")

    if (!is.null(replic_timing_file_prefix)) {
      corr <- data.frame(windowsize = c(names(corr), names(corr_rep)), correlation = c(corr, corr_rep))
      data.table::fwrite(corr, file = gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      corr <- data.frame(windowsize = c(names(corr)), correlation = corr)
      data.table::fwrite(corr, file = gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  } else {
    corr$correlation <- NA
    data.table::fwrite(corr, file = gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep = "\t", quote = FALSE, row.names = FALSE)
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
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")

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
