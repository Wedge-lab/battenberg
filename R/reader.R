########################################################################################
# Generic table reader
########################################################################################
#' Generic reading function using the readr R package, tailored for reading in genomic data
#' @param file Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @param row.names Whether the file contains row names (Default: FALSE)
#' @param stringsAsFactor Legacy parameter that is no longer used (Default: FALSE)
#' @param sep Column separator (Default: \\t)
#' @param chrom_col The column number that contains chromosome denominations. This column will automatically be cast as a character. Should be counted including the row.names (Default: 1)
#' @param skip The number of rows to skip before reading (Default: 0)
#' @return A data frame with contents of the file
#' @export
read_table_generic <- function(file, header = TRUE, stringsAsFactor = FALSE, sep = "\t", chrom_col = 1, skip = 0) {
  # We use a named character vector to force the chromosome column(s) to character
  # This prevents loss of leading zeros or scientific notation issues
  col_classes <- "character"
  names(col_classes) <- as.character(chrom_col)
  log_info("Reading read_table_generic from: {normalizePath(file, mustWork = FALSE)}")

  # fread is the fastest modern parser for large genomic tables
  d <- data.table::fread(
    file = file,
    sep = sep,
    header = header,
    skip = skip,
    colClasses = col_classes,
    check.names = TRUE,
    data.table = TRUE,
    nThread = 4
  )
  log_info("Verified headers generic {paste(colnames(d), collapse = ', ')}")
  return(d)
}


#' Parser for logR data
#' @param filename Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @return A data frame with logR content
read_logr <- function(filename, header = TRUE) {
  log_info("Reading LogR data from: {normalizePath(filename, mustWork = FALSE)}")
  dt <- data.table::fread(
    file = filename,
    header = header,
    colClasses = c("character", "integer", "numeric")
  )
  log_info("Verified headers read_logr {paste(colnames(dt), collapse = ', ')}")
  return(dt)
}

#' Parser for BAF data
#' @param filename Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @return A data frame with BAF content
read_baf_as_data_frame <- function(filename, header = TRUE) {
  log_info("Reading BAF data from: {normalizePath(filename, mustWork = FALSE)}")
  output <- data.table::fread(
    file = filename,
    header = header,
    colClasses = c("character", "integer", "numeric")
  )
  data.table::setDF(output)
  log_info("Verified headers read_baf_as_data_frame {paste(colnames(output), collapse = ', ')}")
  return(output)
}

#' Parser for GC content reference data
#' @param filename Filename of the file to read in
#' @return A data frame with GC content
read_gccontent <- function(filename) {
  log_info("Reading gccontent from: {normalizePath(filename, mustWork = FALSE)}")
  dt <- data.table::fread(
    file = filename,
    header = TRUE,
    sep = "\t",
    skip = "chr",
    check.names = FALSE,
    fill = TRUE,
    select = 1:20
  )
  log_info("Verified headers gccontent {paste(colnames(dt), collapse = ', ')}")
  return(dt)
}

#' Parser for replication timing reference data
#' @param filename Filename of the file to read in
#' @return A data frame with replication timing
read_replication <- function(filename) {
  log_info("Reading replication timing data from: {normalizePath(filename, mustWork = FALSE)}")
  dt <- data.table::fread(
    file = filename,
    header = TRUE,
    sep = "\t",
    skip = "chr"
  )
  log_info("Verified headers replication {paste(colnames(dt), collapse = ', ')}")
  return(dt)
}

#' Parser for BAFsegmented data
#' @param filename Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @return A data frame with BAFsegmented content
read_bafsegmented <- function(filename, header = TRUE) {
  log_info("Reading BAFsegmented data from: {normalizePath(filename, mustWork = FALSE)}")

  dt <- data.table::fread(
    file = filename,
    header = header,
    sep = "\t",
    # Force column types to prevent the coercion warnings
    colClasses = c(Chromosome = "character", Position = "integer")
  )
  # If the file uses 'chr', 'chrom', or 'CHR', we standardize it to 'Chromosome'
  if ("CHR" %in% colnames(dt)) {
    data.table::setnames(dt, "CHR", "Chromosome")
  } else if ("chr" %in% colnames(dt)) {
    data.table::setnames(dt, "chr", "Chromosome")
  }

  log_info("Verified headers bafsegmented: {paste(colnames(dt), collapse = ', ')}")
  return(dt)
}
#' Parser for imputed genotype data
#' @param filename Filename of the file to read in
#' @return A data frame with the imputed genotype output
read_imputed_output <- function(filename) {
  log_info("Reading imputed genotype data from: {normalizePath(filename, mustWork = FALSE)}")
  dt <- data.table::fread(
    file = filename,
    col.names = c("snpidx", "rsidx", "pos", "ref", "alt", "hap1", "hap2"),
    colClasses = c("character", "character", "integer", "character", "character", "integer", "integer"),
    header = FALSE
  )
  log_info("Verified headers read_imputed_output {paste(colnames(dt), collapse = ', ')}")
  return(dt)
}

#' Parser for allele frequencies data
#' @param filename Filename of the file to read in
#' @return A data frame with the alleleCounter output
read_alleleFrequencies <- function(filename) {
  log_info("Reading allele frequencies data from: {normalizePath(filename, mustWork = FALSE)}")
  # skip = "#" handles the comment lines typically found in alleleCounter output
  dt <- data.table::fread(
    file = filename,
    col.names = c("CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth"),
    colClasses = c("character", "integer", "integer", "integer", "integer", "integer", "integer"),
    skip = "#"
  )
  log_info("Verified headers read_alleleFrequencies {paste(colnames(dt), collapse = ', ')}")
  return(dt)
}

#' Parser for impute input data
#' @param filename Filename of the file to read in
#' @return A data frame with the input for impute
read_impute_input <- function(filename) {
  # :: syntax used for log_info or other package calls
  log_info("Reading impute input data from: {normalizePath(filename, mustWork = FALSE)}")

  # Read with data.table for speed
  dt <- data.table::fread(
    file = filename,
    header = FALSE,
    sep = "auto"
  )
  # Convert to data.frame to ensure compatibility with legacy indexing
  dt_df <- as.data.frame(dt)

  # Force column names to start with 'X' instead of 'V'
  # This fixes the 'inp$X6' NULL issue in the Beagle converter
  colnames(dt_df) <- paste0("X", seq_len(ncol(dt_df)))
  log_info("Verified headers read_impute_input: {paste(colnames(dt_df), collapse = ', ')}")
  return(dt_df)
}

#' Parser for beagle5 output data
#' @param filename Filename of the file to read in
#' @return A data frame with the beagle5 output
read_beagle_output <- function(filename) {
  # :: syntax and pure comments
  log_info("Reading beagle5 output data from: {normalizePath(filename, mustWork = FALSE)}")

  # Check if file exists and has content before trying to read
  if (!file.exists(filename) || file.info(filename)$size < 100) {
    log_info("Beagle output file is missing or too small (likely no SNPs phased).")
    # Return an empty data table with the expected structure to prevent dimnames errors
    empty_dt <- data.table::data.table(
      "#CHROM" = character(), POS = integer(), ID = character(),
      REF = character(), ALT = character(), QUAL = character(),
      FILTER = character(), INFO = character(), FORMAT = character(),
      SAMP001 = character()
    )
    return(empty_dt)
  }
  dt <- data.table::fread(
    file = filename,
    skip = "#CHROM",
    header = FALSE
  )

  colnames(dt) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMP001")
  log_info("Successfully read {nrow(dt)} phased SNPs from Beagle output.")
  return(dt)
}

#' Load the rho and psi estimates from a file.
#' @noRd
load_rho_psi_file <- function(rho_psi_file) {
  log_info("Reading rho and psi estimates from: {normalizePath(rho_psi_file, mustWork = FALSE)}")
  rho_psi_info <- data.table::fread(rho_psi_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Always use best solution from grid search - reference segment sometimes gives strange results
  rho <- rho_psi_info$rho[rownames(rho_psi_info) == "FRAC_GENOME"] # rho = tumour percentage (called tp in previous versions)
  psit <- rho_psi_info$psi[rownames(rho_psi_info) == "FRAC_GENOME"] # psi of tumour cells
  goodness <- rho_psi_info$distance[rownames(rho_psi_info) == "FRAC_GENOME"] # goodness of fit
  return(list(rho = rho, psit = psit, goodness = goodness))
}

#' Parse the reference info file
#' @param snp6_reference_info_file A SNP6 reference info master file
#' @noRd
parse_snp6_ref_file <- function(snp6_reference_info_file) {
  log_info("Reading SNP6 reference info from: {normalizePath(snp6_reference_info_file, mustWork = FALSE)}")
  return(data.table::fread(snp6_reference_info_file, header = TRUE))
}

#' Infer the gender using the birdseed report file
#' @param birdseed_report_file The birdseed report file
#' @export
infer_gender_birdseed <- function(birdseed_report_file) {
  log_info("Reading birdseed report from: {normalizePath(birdseed_report_file, mustWork = FALSE)}")
  z <- data.table::fread(birdseed_report_file)
  return(as.character(z$em.cluster.chrX.het.contrast_gender))
}
