########################################################################################
# Concatenate files
########################################################################################
#' Function to concatenate Impute output
#' @noRd
concatenateImputeFiles <- function(inputStart, boundaries) {
  # Generate the list of potential filenames
  # Using paste0 and vectorized division for a bit more speed
  infiles <- paste0(inputStart, "_", boundaries[, 1] / 1000, "K_", boundaries[, 2] / 1000, "K.txt_haps")

  # Filter for existing files with data
  # This uses vectorized checks instead of a for-loop
  existing_files <- infiles[file.exists(infiles) & file.info(infiles)$size > 0]

  # Check if we actually have files to read
  if (length(existing_files) == 0) {
    return(NULL)
  }
  result <- vroom::vroom(existing_files, delim = " ")
  return(data.table::as.data.table(result))
}

#' Function to concatenate allele counter output
#' @noRd
concatenateAlleleCountFiles <- function(inputStart, inputEnd, chr_names) {
  # Vectorized filename generation
  all_files <- paste0(inputStart, chr_names, inputEnd)

  # Vectorized file checking (much faster than a for-loop)
  # This filters the list to only existing, non-empty files
  infiles <- all_files[file.exists(all_files) & file.info(all_files)$size > 0]
  if (length(infiles) == 0) {
    return(data.frame())
  }
  log_info("Using infiles in concatenateAlleleCountFiles: {infiles}")

  # Use rbindlist for the merge
  # We read them as data.tables first (internal to rbindlist)
  # then convert to data.frame at the very end.
  combined <- data.table::rbindlist(
    lapply(infiles, function(f) {
      dt <- read_table_generic(f)
      if (nrow(dt) == 0) {
        log_failure("Allele count file is empty: {f}")
      }
      if (ncol(dt) < 6) {
        log_failure("Allele count file has fewer than 6 columns: {f}")
      }
      return(dt)
    })
  )
  data.table::setDF(combined)
  return(combined)
}

#' Function to concatenate 1000 Genomes SNP reference files
#' @noRd
concatenateG1000SnpFiles <- function(inputStart, inputEnd, chr_names) {
  # Vectorized filename generation
  filenames <- paste0(inputStart, chr_names, inputEnd)
  names(filenames) <- chr_names

  # Filter for valid files
  existing_files <- filenames[file.exists(filenames) & file.info(filenames)$size > 0]

  if (length(existing_files) == 0) {
    return(data.frame())
  }

  # Read files into a named list
  data_list <- lapply(existing_files, read_table_generic)

  # idcol = "chromosome" prepends the list names (chr_names) as the first column
  # This matches the original: cbind(chromosome=chrom, read_table_generic(filename))
  combined <- data.table::rbindlist(data_list, idcol = "chromosome")

  # Convert back to data.frame for index compatibility [[4]]
  data.table::setDF(combined)

  return(combined)
}
