#' Morphs phased SNPs from SNP6 input into haplotype blocks
#'
#' This function matches allele frequencies and halplotype info, reverses frequencies by haplotype, combines the output and saves it to disk.
#' @param chrom The chromosome number for which this function should run.
#' @param alleleFreqFile File containing allele frequency information.
#' @param haplotypeFile File containing haplotype information.
#' @param samplename Identifier for the sample (used in header of output file).
#' @param outputfile Full path pointing to where output should be written.
#' @param chr_names Vector of chromosome names
#' @author dw9
#' @export
GetChromosomeBAFs_SNP6 <- function(chrom, alleleFreqFile, haplotypeFile, samplename, outputfile, chr_names) {
  # Read in the allele frequencies and variant info
  alleleFreqData <- data.table::fread(alleleFreqFile, header = TRUE, data.table = FALSE)
  variant_data <- data.table::fread(haplotypeFile, header = FALSE, data.table = FALSE)

  # Match the two
  alleleFreqData <- alleleFreqData[alleleFreqData[, 1] %in% variant_data[, 3], ]
  select <- match(alleleFreqData[, 1], variant_data[, 3])
  variant_data <- variant_data[select, ]

  chr_name <- chrom
  log_info("Processing: {chr_name}")

  # Switch the haplotypes where required
  alleleFreqs <- alleleFreqData$allele.frequency
  reversedHaplotypes <- variant_data[, 6] == 1
  alleleFreqs[reversedHaplotypes] <- 1.0 - alleleFreqs[reversedHaplotypes]

  log_info("{nrow(variant_data)},{length(alleleFreqs)}")
  # Combine the allele frequencies and variant info and save output
  knownMutBAFs <- cbind(chr_name, variant_data[, 3], alleleFreqs)
  data.table::fwrite(knownMutBAFs, outputfile, sep = "\t", col.names = c("Chromosome", "Position", samplename), quote = FALSE)
}

#' Morphs phased SNPs from WGS input into haplotype blocks
#'
#' @param chrom The chromosome number for which this function is called.
#' @param SNP_file File containing allele counts for each SNP location.
#' @param haplotypeFile File containing impute phasing output.
#' @param samplename Name of the sample (used in header of output file).
#' @param outfile Full path to where the output will be written.
#' @param chr_names Names of all allowed chromosomes as a Vector.
#' @param minCounts An integer describing the minimum number of reads covering this position to be included in the output.
#' @author dw9
#' @importFrom data.table :=
#' @export
GetChromosomeBAFs <- function(
  chrom,
  SNP_file,
  haplotypeFile,
  samplename,
  outfile,
  chr_names,
  minCounts = 1L
) {
  # Input validation
  if (!chrom %in% chr_names) {
    log_failure("chrom must be one of the allowed chromosomes specified in chr_names")
  }
  if (!file.exists(SNP_file)) log_failure("SNP_file not found: {SNP_file}")
  if (!file.exists(haplotypeFile)) log_failure("haplotypeFile not found: {haplotypeFile}")
  minCounts <- as.integer(minCounts)

  log_info("Reading SNP file: {SNP_file}")
  log_info("Reading haplotype file: {haplotypeFile}")
  log_info("Minimum counts: {minCounts} {class(minCounts)}")
  # Load data with explicit column classes to prevent join type mismatches
  # SNP_file (allele frequencies) columns: CHR, POS, A, C, G, T, DEPTH
  snp_dt <- data.table::fread(
    SNP_file,
    sep = "\t",
    header = FALSE,
    skip = "#",
    colClasses = list(character = 1, integer = 2:7)
  )
  # haplotypeFile (phasing) columns: V1..V5 are meta, V6..V7+ are haplotypes. V3 is position.
  phase_dt <- data.table::fread(
    haplotypeFile,
    header = FALSE,
    colClasses = list(integer = 3)
  )

  # If header = FALSE was used but file had a header, the first row might contain NAs
  # due to colClasses. We remove those rows.
  snp_dt <- snp_dt[!is.na(snp_dt[[2]])]
  phase_dt <- phase_dt[!is.na(phase_dt[[3]])]

  # FORCE conversion using character midway to break any factor/weird metadata bonds
  # We use set() to be more robust than := in some parallel environments
  data.table::set(snp_dt, j = "V2", value = as.integer(as.character(snp_dt[["V2"]])))
  data.table::set(phase_dt, j = "V3", value = as.integer(as.character(phase_dt[["V3"]])))

  # Also force count columns to integer to avoid "non-numeric argument" errors later
  for (col in paste0("V", 3:6)) {
    if (col %in% names(snp_dt)) {
      data.table::set(snp_dt, j = col, value = as.integer(as.character(snp_dt[[col]])))
    }
  }

  # Remove any rows that failed conversion
  snp_dt <- snp_dt[!is.na(snp_dt[["V2"]])]
  phase_dt <- phase_dt[!is.na(phase_dt[["V3"]])]

  log_info("VERIFIED types - SNP V2: {class(snp_dt$V2)}, Phase V3: {class(phase_dt$V3)}, SNP V3: {class(snp_dt$V3)}")

  if (nrow(snp_dt) == 0) {
    log_failure("SNP file is empty after filtering/type conversion: {SNP_file}")
  }
  if (nrow(phase_dt) == 0) {
    log_failure("Haplotype file is empty after filtering/type conversion: {haplotypeFile}")
  }

  # Use [[ indexing to explicitly reference columns by name (strings)
  # This avoids "no visible binding" warnings
  het_phase <- phase_dt[phase_dt[["V6"]] != phase_dt[["V7"]]]

  if (nrow(het_phase) == 0) {
    write_empty_output(chrom, samplename, outfile)
    return(invisible(NULL))
  }

  # Match positions using setkeyv (the string-based version of setkey)
  data.table::setkeyv(snp_dt, "V2")

  # Use list() instead of .() to avoid global function warnings
  matched <- snp_dt[list(het_phase[["V3"]]), nomatch = NULL]

  if (nrow(matched) == 0) {
    write_empty_output(chrom, samplename, outfile)
    return(invisible(NULL))
  }

  # Ensure count columns are numeric before matrix conversion
  for (col in names(matched)[3:6]) {
    matched[[col]] <- as.numeric(as.character(matched[[col]]))
  }

  if (nrow(matched) == 0) {
    write_empty_output(chrom, samplename, outfile)
    return(invisible(NULL))
  }

  # Filter het_phase based on matched positions
  het_phase <- het_phase[het_phase[["V3"]] %in% matched[["V2"]]]

  # Map nucleotide characters to column offsets (A=3, C=4, G=5, T=6)
  nuc_to_col <- c(A = 3L, C = 4L, G = 5L, "T" = 6L)

  # Extract phased alleles as characters
  ref_allele <- ifelse(het_phase[["V6"]] == 0, het_phase[["V4"]], het_phase[["V5"]])
  alt_allele <- ifelse(het_phase[["V6"]] == 1, het_phase[["V4"]], het_phase[["V5"]])

  # Use matrix indexing to get counts safely without dynamic column warnings
  # We select only the count columns (3 through 6)
  count_matrix <- as.matrix(matched[, 3:6, with = FALSE])

  # ref_allele and alt_allele map to 1:4 relative to the count_matrix
  ref_idx <- nuc_to_col[ref_allele] - 2L
  alt_idx <- nuc_to_col[alt_allele] - 2L

  row_indices <- seq_len(nrow(count_matrix))
  ref_count <- count_matrix[cbind(row_indices, ref_idx)]
  alt_count <- count_matrix[cbind(row_indices, alt_idx)]

  total_depth <- ref_count + alt_count
  valid <- total_depth >= minCounts

  if (!any(valid)) {
    write_empty_output(chrom, samplename, outfile)
    return(invisible(NULL))
  }

  # Construct output data.table
  output_dt <- data.table::data.table(
    Chromosome = chrom,
    Position   = matched[["V2"]][valid],
    BAF        = alt_count[valid] / total_depth[valid]
  )
  data.table::setnames(output_dt, "BAF", samplename)
  data.table::fwrite(output_dt, file = outfile, sep = "\t", quote = FALSE)
}

# Helper function to avoid code duplication
write_empty_output <- function(chrom, samplename, outfile) {
  empty_dt <- data.table::data.table(
    Chromosome = character(),
    Position   = integer(),
    dummy      = numeric()
  )
  data.table::setnames(empty_dt, "dummy", samplename)
  data.table::fwrite(empty_dt, file = outfile, sep = "\t", quote = FALSE)
}

#' Plot haplotyped BAF values for a single chromosome
#'
#' Reads a tab-separated file produced by GetChromosomeBAFs() (columns: Chromosome, Position, <samplename>)
#' and creates a high-resolution PNG showing the B Allele Frequency (BAF) mirrored around 0.5
#' (standard haplotype/ASCAT-style plot).
#'
#' @param haplotyped_baf_file Path to the input TSV file with haplotyped BAF data.
#' @param image_file_name Path to the output PNG file.
#' @param samplename Name of the sample (used in plot title).
#' @param chrom Chromosome identifier (used only for validation and title if data is empty).
#'
#' @return Invisibly returns NULL; side effect is writing the PNG file.
#' @author Original: dw9; Modernized version
#' @export
plot_haplotype_data <- function(haplotyped_baf_file,
                                image_file_name,
                                samplename,
                                chrom) {
  # Input validation
  if (!file.exists(haplotyped_baf_file)) {
    log_failure("Input file not found: ", haplotyped_baf_file)
  }

  # Read data (expecting columns: Chromosome, Position, <samplename>)
  baf_dt <- data.table::fread(haplotyped_baf_file, header = TRUE)

  # Determine x-axis limits
  if (nrow(baf_dt) == 0) {
    log_info("No data in '{haplotyped_baf_file}' - creating empty plot")
    x_min <- 1
    x_max <- 2
    positions <- numeric()
    baf_values <- numeric()
    plot_chrom <- chrom
  } else {
    x_min <- min(baf_dt$Position, na.rm = TRUE)
    x_max <- max(baf_dt$Position, na.rm = TRUE)
    positions <- baf_dt$Position
    # third column is the sample BAF
    baf_values <- baf_dt[[3]]
    plot_chrom <- baf_dt$Chromosome[1]
  }

  # Open PNG device with reasonable size and resolution
  grDevices::png(
    filename = image_file_name,
    width = 1200, height = 600, res = 150, type = "cairo"
  )

  # Assuming create_haplotype_plot is a custom function available in your package/environment
  create_haplotype_plot(
    chrom_position = positions,
    points.blue    = baf_values,
    points.red     = 1 - baf_values,
    x_min          = x_min,
    x_max          = x_max,
    title          = paste(samplename, ", chromosome", plot_chrom),
    xlab           = "Position",
    ylab           = "BAF"
  )

  grDevices::dev.off()
  invisible(NULL)
}
#' Combine per-chromosome BAF files into a single table
#'
#' @param prefix   File path prefix before chromosome name
#' @param suffix   File path suffix after chromosome name
#' @param chroms   Character vector of chromosome names
#' @param output   Path to output TSV file
#'
#' @return Invisibly returns the combined data.frame
#' @export
concatenate_baf_files <- function(
  input_start,
  input_end,
  output_file,
  chr_names
) {
  files <- paste0(input_start, chr_names, input_end)

  log_info("Starting concatenation for {length(chr_names)} expected BAF files {files}")


  # Filter for existing and non-empty files
  valid_files <- files[
    fs::file_exists(files) &
      fs::file_size(files) > 0
  ]

  exists_mask <- fs::file_exists(files)
  size_mask <- fs::file_size(files) > 0
  missing_chrs <- chr_names[!exists_mask]
  if (base::length(missing_chrs) > 0) {
    log_info("Chromosomes missing files: {base::paste(missing_chrs, collapse = ', ')}")
  }

  empty_chrs <- chr_names[exists_mask & !size_mask]
  if (base::length(empty_chrs) > 0) {
    log_info("DATA ISSUE: Chromosomes with 0-byte files: {base::paste(empty_chrs, collapse = ', ')}")
  }

  valid_files <- files[exists_mask & size_mask]

  if (base::length(valid_files) == 0) {
    log_info("CRITICAL: Zero valid BAF files found across all chromosomes.")
  }

  log_info("Proceeding to combine {base::length(valid_files)} valid files")

  # Force first column to character
  # Use column index 1 to avoid needing names(vroom(...)) twice
  first_file_cols <- names(vroom::vroom(
    valid_files[1],
    n_max = 0,
    progress = FALSE,
    show_col_types = FALSE
  ))
  col_spec <- vroom::cols(
    .default = vroom::col_guess(),
    !!!stats::setNames(list(vroom::col_character()), first_file_cols[1])
  )

  log_info("Reading data using column spec based on {fs::path_file(valid_files[1])}")

  combined <- vroom::vroom(
    valid_files,
    id = "file_path",
    delim = "\t",
    col_types = col_spec,
    progress = FALSE,
    show_col_types = FALSE,
    .name_repair = "minimal"
  ) |>
    dplyr::select(-dplyr::any_of("file_path"))

  total_rows <- base::nrow(combined)
  if (total_rows == 0) {
    log_failure("DATA ISSUE: Files were read but the resulting table is empty.")
  }

  log_info("Total combined rows: {base::format(total_rows, big.mark = ',')}")

  # Ensure output directory exists
  fs::dir_create(fs::path_dir(output_file), recurse = TRUE)

  # Write output
  vroom::vroom_write(
    combined,
    path = output_file,
    delim = "\t",
    na = "NA",
    quote = "none"
  )

  log_info("BAF concatenation complete. Final file size: {fs::file_size(output_file)}")
}
