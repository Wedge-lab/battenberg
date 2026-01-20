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
      if (!is.na(m_baf) && m_baf != 0 && m_baf != 1) {
        NEW_BAF <- m_baf
      } else {
        NEW_BAF <- collapse::fmean(baf_vals, na.rm = TRUE)
      }
    }

    # LogR update with safety for empty segments
    logr_idx <- S4Vectors::subjectHits(
      GenomicRanges::findOverlaps(subclones[INDEX], logR)
    )

    if (length(logr_idx) == 0) {
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
    subclones$ID <- seq_along(subclones)

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
            logr_p <- fast_p(logr_curr, logr_other)
            baf_p <- fast_p(baf_curr, baf_other)
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
  to_mask_idx <- which(subclones$nMaj1_A > max_allowed_state | subclones$nMin1_A > max_allowed_state)

  if (length(to_mask_idx) == 0) {
    return(list(
      subclones = subclones,
      bafsegmented = bafsegmented,
      masked_count = 0,
      masked_size = 0
    ))
  }

  count <- length(to_mask_idx)
  masked_size <- sum(subclones$endpos[to_mask_idx] - subclones$startpos[to_mask_idx])

  # Identify segments to mask in the BAFsegmented file
  # Use GenomicRanges for O(N+M) overlap detection instead of the O(N*M) loop
  segs_to_mask <- subclones[to_mask_idx, ]
  gr_segs <- GenomicRanges::GRanges(
    seqnames = segs_to_mask$chr,
    # Original logic: startpos < Position <= endpos
    ranges = IRanges::IRanges(start = segs_to_mask$startpos + 1, end = segs_to_mask$endpos)
  )

  gr_snps <- GenomicRanges::GRanges(
    seqnames = bafsegmented$Chromosome,
    ranges = IRanges::IRanges(start = bafsegmented$Position, end = bafsegmented$Position)
  )

  # Find SNPs that fall within any masked segment
  overlaps <- GenomicRanges::findOverlaps(gr_snps, gr_segs)
  if (length(overlaps) > 0) {
    bafsegmented$BAFseg[unique(S4Vectors::queryHits(overlaps))] <- NA
  }

  # Now mask the subclones table
  subclones[to_mask_idx, c("nMaj1_A", "nMin1_A", "nMaj2_A", "nMin2_A")] <- NA

  return(list(
    subclones = subclones,
    bafsegmented = bafsegmented,
    masked_count = count,
    masked_size = masked_size
  ))
}
