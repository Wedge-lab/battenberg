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

  # Use rbindlist for the merge
  # We read them as data.tables first (internal to rbindlist)
  # then convert to data.frame at the very end.
  combined <- data.table::rbindlist(
    lapply(infiles, read_table_generic())
  )
  return(as.data.frame(combined))
}

#' Function to concatenate 1000 Genomes SNP reference files
#' @noRd
concatenateG1000SnpFiles <- function(inputStart, inputEnd, chr_names) {
  # Generate all potential filenames at once
  filenames <- paste0(inputStart, chr_names, inputEnd)
  names(filenames) <- chr_names # Keep names so rbindlist knows the ID
  # Filter for files that exist and are not empty
  existing_files <- filenames[file.exists(filenames) & file.info(filenames)$size > 0]

  if (length(existing_files) == 0) {
    return(data.frame())
  }

  # read_table_generic should ideally return a data.table for this to be fastest
  # We use lapply to read them into a list
  data_list <- lapply(existing_files, read_table_generic)
  # rbindlist with 'idcol' automatically creates the 'chromosome' column
  # based on the names of our list (which are the chr_names)
  combined <- data.table::rbindlist(data_list, idcol = "chromosome")

  return(as.data.frame(combined))
}
