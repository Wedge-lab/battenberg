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

  # fread is the fastest modern parser for large genomic tables
  d <- data.table::fread(
    file = file,
    sep = sep,
    header = header,
    skip = skip,
    colClasses = col_classes,
    check.names = TRUE,
    data.table = FALSE,
    nThread = 4
  )

  return(d)
}


#' Parser for logR data
#' @param filename Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @return A data frame with logR content
read_logr <- function(filename, header = TRUE) {
  data.table::fread(
    file = filename,
    header = header,
    colClasses = c("character", "integer", "numeric")
  )
}

#' Parser for BAF data
#' @param filename Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @return A data frame with BAF content
read_baf_as_data_frame <- function(filename, header = TRUE) {
  output <- data.table::fread(
    file = filename,
    header = header,
    colClasses = c("character", "integer", "numeric")
  )
  data.table::setDF(output)
  return(output)
}

#' Parser for GC content reference data
#' @param filename Filename of the file to read in
#' @return A data frame with GC content
read_gccontent <- function(filename) {
  data.table::fread(
    file = filename,
    skip = 1,
    header = FALSE,
    select = 2:14,
    colClasses = list(character = 2, integer = 3, numeric = 4:14)
  )
}

#' Parser for replication timing reference data
#' @param filename Filename of the file to read in
#' @return A data frame with replication timing
read_replication <- function(filename) {
  data.table::fread(
    file = filename,
    header = FALSE,
    colClasses = list(character = 1, integer = 2, numeric = 3:17)
  )
}

#' Parser for BAFsegmented data
#' @param filename Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @return A data frame with BAFsegmented content
read_bafsegmented <- function(filename, header = TRUE) {
  data.table::fread(
    file = filename,
    header = header,
    sep = "\t",
    colClasses = c("character", "integer", "numeric", "numeric", "numeric")
  )
}

#' Parser for imputed genotype data
#' @param filename Filename of the file to read in
#' @return A data frame with the imputed genotype output
read_imputed_output <- function(filename) {
  data.table::fread(
    file = filename,
    col_names = c("snpidx", "rsidx", "pos", "ref", "alt", "hap1", "hap2"),
    colClasses = c("character", "character", "integer", "character", "character", "integer", "integer"),
    header = FALSE
  )
}

#' Parser for allele frequencies data
#' @param filename Filename of the file to read in
#' @return A data frame with the alleleCounter output
read_alleleFrequencies <- function(filename) {
  # skip = "#" handles the comment lines typically found in alleleCounter output
  data.table::fread(
    file = filename,
    col_names = c("CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth"),
    colClasses = c("character", "integer", "integer", "integer", "integer", "integer", "integer"),
    skip = "#"
  )
}

#' Parser for impute input data
#' @param filename Filename of the file to read in
#' @return A data frame with the input for impute
read_impute_input <- function(filename) {
  # Automatically detects delimiters (like space or tab) while forcing column types
  data.table::fread(
    file = filename,
    col_names = NULL, # Uses default or looks for header
    colClasses = c("character", "character", "integer", "character", "character", "integer", "integer", "integer"),
    header = FALSE
  )
}

#' Parser for beagle5 output data
#' @param filename Filename of the file to read in
#' @return A data frame with the beagle5 output
read_beagle_output <- function(filename) {
  # Efficiently skips VCF-style headers using the '#' skip pattern
  data.table::fread(
    file = filename,
    col_names = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMP001"),
    colClasses = c("character", "integer", "character", "character", "character", "character", "character", "character", "character", "character"),
    skip = "#"
  )
}

#' Load the rho and psi estimates from a file.
#' @noRd
load_rho_psi_file <- function(rho_psi_file) {
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
  return(data.table::fread(snp6_reference_info_file, header = TRUE))
}

#' Infer the gender using the birdseed report file
#' @param birdseed_report_file The birdseed report file
#' @export
infer_gender_birdseed <- function(birdseed_report_file) {
  z <- data.table::fread(birdseed_report_file)
  return(as.character(z$em.cluster.chrX.het.contrast_gender))
}
