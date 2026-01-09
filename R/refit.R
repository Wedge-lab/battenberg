########################################################################################
# Refitting functions
########################################################################################
#' Calculate rho and psi values from a refit suggestion
#'
#' Use this function to calculate the refit values from a refit suggestion.
#' @param refBAF BAF of the segment
#' @param refLogR logR of the segment
#' @param refMajor Major allele copy number
#' @param refMinor Minor allele copy number
#' @param rho Sample rho parameter
#' @param gamma_param Platform gamma parameter
#' @return A list with a field for rho and psi_t
#' @author sd11
#' @export
calc_rho_psi_refit <- function(refBAF, refLogR, refMajor, refMinor, rho, gamma_param) {
  rho <- (2 * refBAF - 1) / (2 * refBAF - refBAF * (refMajor + refMinor) - 1 + refMajor)
  psi <- (rho * (refMajor + refMinor) + 2 - 2 * rho) / (2^(refLogR / gamma_param))
  psi_t <- psi2psit(rho, psi)
  return(list(rho = rho, psi_t = psi_t))
}

#' Calculate refit values from a refit suggestion
#'
#' Use this function to calculate the refit values from a refit suggestion.
#' @param subclones_file A Battenberg subclones.txt file
#' @param segment_chrom Chromsome of the segment to use for refitting
#' @param segment_pos Position within the start/end coordinates of the segment to use for refitting
#' @param new_nMaj Major allele copy number
#' @param new_nMin Minor allele copy number
#' @param rho Sample rho parameter
#' @param gamma_param Platform gamma parameter
#' @return A list with a field for rho and psi_t
#' @author sd11
#' @export
suggest_refit <- function(subclones_file, segment_chrom, segment_pos, new_nMaj, new_nMin, rho, gamma_param) {
  subclones <- data.table::fread(subclones_file, header = TRUE, stringsAsFactors = FALSE)
  segment <- subclones[subclones$chr == segment_chrom & subclones$startpos <= segment_pos & subclones$endpos >= segment_pos, ]
  segment_BAF <- segment$BAF
  segment_LogR <- segment$LogR
  return(calc_rho_psi_refit(segment_BAF, segment_LogR, new_nMaj, new_nMin, rho, gamma_param))
}

#' Create refit suggestions for a fit copy number profile
#'
#' This function takes a fit copy number profile and generates refit suggestions for a future rerun.
#' If there are clonal alterations above a specified size, then those written out as supplied as suggestions,
#' otherwise a refit suggestion of an external purity value will be saved.
#' @param samplename Samplename for the output file
#' @param subclones_file File containing a fit copy number profile
#' @param rho_psi_file File with rho and psi values
#' @param gamma_param Platform gamma parameter
#' @param min_segment_size_mb Minimum size of a segment in Mb to be considered for a refit suggestion (Default: 2)
#' @author sd11
#' @export
cnfit_to_refit_suggestions <- function(samplename, subclones_file, rho_psi_file, gamma_param, min_segment_size_mb = 2) {
  subclones <- read_table_generic(subclones_file)
  subclones$len <- subclones$endpos / 1000000 - subclones$startpos / 1000000
  subclones$is_cna <- subclones$nMaj1_A != subclones$nMin1_A

  print(min_segment_size_mb)
  print(subclones$is_cna)
  if (any(subclones$len > min_segment_size_mb & subclones$is_cna)) {
    # There are large scale alterations, save the top couple as suggestions
    rho_psi <- utils::read.table(rho_psi_file, header = TRUE, stringsAsFactors = FALSE)
    rho <- rho_psi["FRAC_GENOME", "rho"]
    psi_t <- rho_psi["FRAC_GENOME", "psi"]

    # Take only segments that are clonal and are an alteration
    is_subclonal <- subclones$frac1_A < 1
    subclones_clonal_cna <- subset(subclones, !is_subclonal & subclones$is_cna)
    subclones_clonal_cna <- subclones_clonal_cna[with(subclones_clonal_cna, order(len, decreasing = TRUE)), ]

    if (nrow(subclones_clonal_cna) == 0) {
      output <- data.table::data.table(
        project = NA, samplename = samplename,
        qc = NA, cellularity_refit = TRUE,
        chrom = NA, pos = NA, maj = NA,
        min = NA, baf = NA, logr = NA,
        rho_estimate = NA, psi_t_estimate = NA,
        rho_diff = NA, psi_t_diff = NA
      )
      data.table::setDF(output)
    } else {
      # Generate a couple of solutions, but not more than are possibly available
      max_solutions <- ifelse(nrow(subclones_clonal_cna) >= 5, 5, nrow(subclones_clonal_cna))
      subclones_clonal_cna <- subclones_clonal_cna[1:max_solutions, , drop = FALSE]

      # Determine position in Mb within the segment
      position <- subclones_clonal_cna$startpos + (subclones_clonal_cna$endpos - subclones_clonal_cna$startpos) / 2
      position <- position / 1000000
      position_round_up <- ceiling(position)
      position_round_down <- floor(position)
      position <- ifelse(position_round_up < subclones_clonal_cna$endpos, position_round_up, position_round_down)

      output <- data.frame(
        project = rep(NA, max_solutions),
        samplename = rep(samplename, max_solutions),
        qc = rep(NA, max_solutions),
        cellularity_refit = rep(FALSE, max_solutions),
        chrom = subclones_clonal_cna$chr[1:max_solutions],
        pos = paste(position, "M", sep = ""),
        maj = subclones_clonal_cna$nMaj1_A[1:max_solutions],
        min = subclones_clonal_cna$nMin1_A[1:max_solutions],
        baf = subclones_clonal_cna$BAF[1:max_solutions],
        logr = subclones_clonal_cna$LogR[1:max_solutions]
      )

      # refBAF, refLogR, refMajor, refMinor, rho, gamma_param
      res <- calc_rho_psi_refit(output$baf, output$logr, output$maj, output$min, rho, gamma_param)
      output$rho_estimate <- res$rho
      output$psi_t_estimate <- res$psi_t
      output$rho_diff <- abs(rho - output$rho_estimate)
      output$psi_t_diff <- abs(psi_t - output$psi_t_estimate)
    }
  } else {
    # No large clonal alteration, save a suggestion that should use an external purity value
    output <- data.frame(project = NA, samplename = samplename, qc = NA, cellularity_refit = TRUE, chrom = NA, pos = NA, maj = NA, min = NA, baf = NA, logr = NA, rho_estimate = NA, psi_t_estimate = NA, rho_diff = NA, psi_t_diff = NA)
  }
  data.table::fwrite(output, file = paste0(samplename, "_refit_suggestion.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
}
