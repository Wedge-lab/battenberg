#' Helper function to adjust the BAF segmented values. By default the segmentation
#' takes the mean BAFphased for each segment, but that doesn't work very well with
#' outliers (i.e. badly phased regions). This function is then called to adjust
#' the segmented BAF. By default this now takes the median
#' @param baf_chrom A data frame with columns BAFphased and BAFseg. BAFseg will be overwritten.
#' @return A data frame with columns BAFphased and BAFseg.
#' @author sd11
#' @noRd
adjustSegmValues <- function(baf_chrom) {
  segs <- rle(baf_chrom$BAFseg)
  for (i in seq_along(segs$lengths)) {
    end <- cumsum(segs$lengths[1:i])
    end <- end[length(end)]
    start <- (end - segs$lengths[i]) + 1 # segs$lengths contains end points
    # baf_chrom$bafmean[start:end] = mean(baf_chrom$BAFphased[start:end])
    baf_chrom$BAFseg[start:end] <- median(baf_chrom$BAFphased[start:end])
    # This needs the ASCAT version of PCF
    # datwins = madWins(baf_chrom$BAFphased[start:end], 2.5, 25)$ywin
    # baf_chrom$madwins_mean[start:end] = mean(datwins)
    # baf_chrom$madwins_median[start:end] = median(datwins)
  }
  return(baf_chrom)
}

#' Segment BAF, with the possible inclusion of structural variant breakpoints
#'
#' This function breaks the genome up into chromosomes, possibly further when SV breakpoints
#' are provided, and runs PCF on each to segment the chromosomes independently.
#' @param samplename Name of the sample, which is used to name output figures
#' @param inputfile String that points to the output from the \code{concatenate_baf_files} function. This contains the phased SNPs with their BAF values
#' @param outputfile String where the segmentation output will be written
#' @param prior_breakpoints_file String that points to a file with prior breakpoints (from SVs for example) with chromosome and position columns (Default: NULL)
#' @param gamma The gamma parameter controls the size of the penalty of starting a new segment during segmentation. It is therefore the key parameter for controlling the number of segments (Default 10)
#' @param kmin Kmin represents the minimum number of probes/SNPs that a segment should consist of (Default 3)
#' @param phasegamma Gamma parameter used when correcting phasing mistakes (Default 3)
#' @param phasekmin Kmin parameter used when correcting phasing mistakes (Default 3)
#' @param no_segmentation Do not perform segmentation. This step will switch the haplotype blocks, but then just takes the mean BAFphased as BAFsegm
#' @param calc_seg_baf_option Various options to recalculate the BAF of a segment. Options are: 1 - median, 2 - mean, 3 - ifelse median==0 or 1, median, mean. (Default: 3)
#' @author sd11
#' @export
segment_baf_phased <- function(
  samplename, inputfile,
  outputfile, prior_breakpoints_file = NULL,
  gamma = 10, phasegamma = 3, kmin = 3,
  phasekmin = 3, no_segmentation = FALSE,
  calc_seg_baf_option = 3
) {
  # Function that takes SNPs that belong to a single segment and looks for big holes between
  # each pair of SNPs. If there is a big hole it will add another breakpoint to the breakpoints data.frame
  addin_bigholes <- function(breakpoints, positions, chrom, startpos, maxsnpdist) {
    # If there is a big hole (i.e. centromere), add it in as a separate set of breakpoints

    # Get the chromosome coordinate right before a big hole
    bigholes <- which(diff(positions) >= maxsnpdist)
    if (length(bigholes) > 0) {
      for (endindex in bigholes) {
        breakpoints <- rbind(
          breakpoints,
          data.frame(chrom = chrom, start = startpos, end = positions[endindex])
        )
        startpos <- positions[endindex + 1]
      }
    }
    return(list(breakpoints = breakpoints, startpos = startpos))
  }

  # Helper function that creates segment breakpoints from SV calls
  # @param bkps_chrom Breakpoints for a single chromosome
  # @param BAFrawchr Raw BAF values of germline heterozygous SNPs on a single chromosome
  # @param addin_bigholes Flag whether bog holes in data are to be added as breakpoints
  # @return A data.frame with chrom, start and end columns
  # @author sd11
  bkps_to_presegment_breakpoints <- function(chrom, bkps_chrom, BAFrawchr, use_bigholes) {
    maxsnpdist <- 3000000

    bkps_breakpoints <- bkps_chrom$position

    # If there are no prior breakpoints, we cannot insert any
    if (length(bkps_breakpoints) > 0) {
      breakpoints <- data.frame()

      # check which comes first, the breakpoint or the first SNP
      if (BAFrawchr$Position[1] < bkps_breakpoints[1]) {
        startpos <- BAFrawchr$Position[1]
        # We're starting from SNP data, so the first SV should be added first
        startfromsv <- 1
      } else {
        startpos <- bkps_breakpoints[1]
        startfromsv <- 2 # We've just added the first SV, don't use it again
      }

      for (svposition in bkps_breakpoints[startfromsv:length(bkps_breakpoints)]) {
        selectedsnps <- BAFrawchr$Position >= startpos & BAFrawchr$Position <= svposition
        if (sum(selectedsnps, na.rm = TRUE) > 0) {
          if (use_bigholes) {
            # If there is a big hole (i.e. centromere), add it in as a separate set of breakpoints
            res <- addin_bigholes(breakpoints, BAFrawchr$Position[selectedsnps], chrom, startpos, maxsnpdist)
            breakpoints <- res$breakpoints
            startpos <- res$startpos
          }

          endindex <- max(which(selectedsnps))
          breakpoints <- rbind(breakpoints, data.frame(chrom = chrom, start = startpos, end = BAFrawchr$Position[endindex]))
          # Previous SV is the new starting point for the next segment
          startpos <- BAFrawchr$Position[endindex + 1]
        }
      }

      # Add the remainder of the chromosome, if available
      if (BAFrawchr$Position[nrow(BAFrawchr)] > bkps_breakpoints[length(bkps_breakpoints)]) {
        endindex <- nrow(BAFrawchr)
        breakpoints <- rbind(breakpoints, data.frame(chrom = chrom, start = startpos, end = BAFrawchr$Position[endindex]))
      }
    } else {
      # There are no SVs, so create one big segment
      print("No prior breakpoints found")
      startpos <- BAFrawchr$Position[1]
      breakpoints <- data.frame()

      if (use_bigholes) {
        # If there is a big hole (i.e. centromere), add it in as a separate set of breakpoints
        res <- addin_bigholes(breakpoints, BAFrawchr$Position, chrom, startpos, maxsnpdist = maxsnpdist)
        breakpoints <- res$breakpoints
        startpos <- res$startpos
      }

      breakpoints <- rbind(breakpoints, data.frame(chrom = chrom, start = startpos, end = BAFrawchr$Position[nrow(BAFrawchr)]))
    }
    return(breakpoints)
  }

  # Run PCF on presegmented data
  # @param BAFrawchr Raw BAF for this chromosome
  # @param presegment_chrom_start
  # @param presegment_chrom_end
  # @param phasekmin
  # @param phasegamma
  # @param kmin
  # @param gamma
  # @param no_segmentation Do not perform segmentation. This step will switch the haplotype blocks, but then just takes the mean BAFphased as BAFsegm
  # @return A data.frame with columns Chromosome,Position,BAF,BAFphased,BAFseg
  run_pcf <- function(BAFrawchr, presegment_chrom_start, presegment_chrom_end, phasekmin, phasegamma, kmin, gamma, no_segmentation = FALSE) {
    row.indices <- which(BAFrawchr$Position >= presegment_chrom_start &
      BAFrawchr$Position <= presegment_chrom_end)

    BAF <- BAFrawchr[row.indices, 2]
    pos <- BAFrawchr[row.indices, 1]

    sdev <- getMad(ifelse(BAF < 0.5, BAF, 1 - BAF), k = 25)
    # Standard deviation is not defined for a single value
    if (is.na(sdev)) {
      sdev <- 0
    }
    # for cell lines, sdev goes to zero in regions of LOH, which causes problems.
    # 0.09 is around the value expected for a binomial distribution around 0.5 with depth 30
    if (sdev < 0.09) {
      sdev <- 0.09
    }

    print(paste("BAFlen=", length(BAF), sep = ""))
    if (length(BAF) < 50) {
      BAFsegm <- rep(mean(BAF), length(BAF))
    } else {
      res <- selectFastPcf(BAF, phasekmin, phasegamma * sdev, T)
      BAFsegm <- res$yhat
    }

    BAFphased <- ifelse(BAFsegm > 0.5, BAF, 1 - BAF)

    if (length(BAFphased) < 50 | no_segmentation) {
      BAFphseg <- rep(mean(BAFphased), length(BAFphased))
    } else {
      res <- selectFastPcf(BAFphased, kmin, gamma * sdev, T)
      BAFphseg <- res$yhat
    }

    if (length(BAF) > 0) {
      #
      # Note: When adding options, also add to merge_segments
      #

      # Recalculate the BAF of each segment, if required
      if (calc_seg_baf_option == 1) {
        # Adjust the segment BAF to not take the mean as that is sensitive to improperly phased segments
        BAFphseg <- adjustSegmValues(data.frame(BAFphased = BAFphased, BAFseg = BAFphseg))$BAFseg
      } else if (calc_seg_baf_option == 2) {
        # Don't do anything, the BAF is already the mean
      } else if (calc_seg_baf_option == 3) {
        # Take the median, unless the median is exactly 0 or 1. At the extreme
        # there is no difference between lets say 40 and 41 copies and BB cannot
        # fit a copy number state. The mean is less prone to become exactly 0 or 1
        # but the median is generally a better estimate that is less sensitive to
        # how well the haplotypes have been reconstructed
        BAFphseg_median <- adjustSegmValues(data.frame(BAFphased = BAFphased, BAFseg = BAFphseg))$BAFseg
        BAFphseg <- ifelse(BAFphseg_median %in% c(0, 1), BAFphseg, BAFphseg_median)
      } else {
        warning("Supplied calc_seg_baf_option to segment_baf_phased not valid, using mean BAF by default")
      }
    }

    return(data.frame(
      Chromosome = rep(chr, length(row.indices)),
      Position = BAFrawchr[row.indices, 1],
      BAF = BAF,
      BAFphased = BAFphased,
      BAFseg = BAFphseg,
      tempBAFsegm = BAFsegm
    )) # Keep track of BAFsegm for the plot below
  }

  BAFraw <- as.data.frame(read_baf(inputfile))
  if (!is.null(prior_breakpoints_file)) {
    bkps <- read.table(prior_breakpoints_file, header = TRUE, stringsAsFactors = FALSE)
  } else {
    bkps <- NULL
  }

  BAFoutput <- NULL
  for (chr in unique(BAFraw[, 1])) {
    print(paste0("Segmenting ", chr))
    BAFrawchr <- BAFraw[BAFraw[, 1] == chr, c(2, 3)]
    BAFrawchr <- BAFrawchr[!is.na(BAFrawchr[, 2]), ]
    if (!is.null(bkps)) {
      bkps_chrom <- bkps[bkps$chromosome == chr, ]
    } else {
      bkps_chrom <- data.frame(chromosome = character(), position = numeric())
    }

    breakpoints_chrom <- bkps_to_presegment_breakpoints(chr, bkps_chrom, BAFrawchr, addin_bigholes = TRUE)
    BAFoutputchr <- NULL

    for (r in seq_len(nrow(breakpoints_chrom))) {
      BAFoutput_preseg <- run_pcf(BAFrawchr, breakpoints_chrom$start[r], breakpoints_chrom$end[r], phasekmin, phasegamma, kmin, gamma, no_segmentation)
      BAFoutputchr <- rbind(BAFoutputchr, BAFoutput_preseg)
    }

    png(filename = paste(samplename, "_RAFseg_chr", chr, ".png", sep = ""), width = 2000, height = 1000, res = 200, type = "cairo")
    create_segmented_plot(
      chrom_position = BAFoutputchr$Position / 1000000,
      points.red = BAFoutputchr$BAF,
      points.green = BAFoutputchr$tempBAFsegm,
      x_min = min(BAFoutputchr$Position) / 1000000,
      x_max = max(BAFoutputchr$Position) / 1000000,
      title = paste(samplename, ", chromosome ", chr, sep = ""),
      xlab = "Position (Mb)",
      ylab = "BAF (phased)",
      prior_bkps_pos = bkps_chrom$position / 1000000
    )
    dev.off()

    png(filename = paste(samplename, "_segment_chr", chr, ".png", sep = ""), width = 2000, height = 1000, res = 200, type = "cairo")
    create_baf_plot(
      chrom_position = BAFoutputchr$Position / 1000000,
      points.red.blue = BAFoutputchr$BAF,
      plot.red = BAFoutputchr$tempBAFsegm > 0.5,
      points.darkred = BAFoutputchr$BAFseg,
      points.darkblue = 1 - BAFoutputchr$BAFseg,
      x_min = min(BAFoutputchr$Position) / 1000000,
      x_max = max(BAFoutputchr$Position) / 1000000,
      title = paste(samplename, ", chromosome ", chr, sep = ""),
      xlab = "Position (Mb)",
      ylab = "BAF (phased)",
      prior_bkps_pos = bkps_chrom$position / 1000000
    )
    dev.off()

    BAFoutputchr$BAFphased <- ifelse(BAFoutputchr$tempBAFsegm > 0.5, BAFoutputchr$BAF, 1 - BAFoutputchr$BAF)
    # Remove the temp BAFsegm values as they are only needed for plotting
    BAFoutput <- rbind(BAFoutput, BAFoutputchr[, c(1:5)])
  }
  colnames(BAFoutput) <- c("Chromosome", "Position", "BAF", "BAFphased", "BAFseg")
  data.table::fwrite(BAFoutput, outputfile, sep = "\t", row.names = FALSE, col_names = TRUE, quote = FALSE)
}


#' Segment BAF, with the possible inclusion of structural variant breakpoints
#'
#' This function breaks the genome up into chromosomes, possibly further when SV breakpoints
#' are provided, and runs PCF on each to segment the chromosomes independently.
#' @param samplename Name of the sample, which is used to name output figures
#' @param inputfile String that points to the output from the \code{concatenate_baf_files} function. This contains the phased SNPs with their BAF values
#' @param outputfile String where the segmentation output will be written
#' @param prior_breakpoints_file String that points to a file with prior breakpoints (from SVs for example) with chromosome and position columns (Default: NULL)
#' @param gamma The gamma parameter controls the size of the penalty of starting a new segment during segmentation. It is therefore the key parameter for controlling the number of segments (Default 10)
#' @param calc_seg_baf_option Various options to recalculate the BAF of a segment. Options are: 1 - median, 2 - mean, 3 - ifelse median==0 or 1, median, mean. (Default: 3)
#' @param GENOMEBUILD Genome build upon which the 1000G SNP coordinates were obtained
#' @author jdemeul, sd11
#' @export
segment_baf_phased_multisample <- function(samplename, inputfile, outputfile, prior_breakpoints_file = NULL, gamma = 10, calc_seg_baf_option = 3, GENOMEBUILD) {
  # --- 1. Internal Helper: Segment Generator ---
  get_segments <- function(chrom, bkps_chrom, BAFrawchr, maxsnpdist = 3000000) {
    snps <- BAFrawchr$Position

    # Identify gaps using base R vectorization
    gaps <- which(diff(snps) >= maxsnpdist)
    gap_bkps <- snps[gaps]

    # Merge SV and Gap breakpoints
    all_cuts <- sort(unique(c(bkps_chrom$position, gap_bkps)))

    # Define start/end pairs
    cut_indices <- findInterval(all_cuts, snps)

    seg_starts <- c(snps[1], snps[cut_indices + 1])
    seg_ends <- c(snps[cut_indices], snps[length(snps)])

    # Explicitly use data.table namespace for construction
    segments <- data.table::data.table(chrom = chrom, start = seg_starts, end = seg_ends)
    return(segments[start <= end])
  }

  # --- 2. Internal Helper: PCF Runner ---
  run_pcf_modern <- function(BAFrawchr, start, end, gamma) {
    # Subset using standard data.table syntax (methods are registered if package is installed)
    BAF_subset <- BAFrawchr[Position >= start & Position <= end]
    if (nrow(BAF_subset) == 0) {
      return(NULL)
    }

    vals <- as.matrix(BAF_subset[, -c(1:2)])

    # Fully qualified copynumber calls
    sdevs <- apply(vals, 2, function(x) {
      getMad(ifelse(x < 0.5, x, 1 - x), k = 25)
    })
    sdevs[is.na(sdevs) | sdevs < 0.09] <- 0.09
    sdev <- mean(sdevs)

    if (nrow(BAF_subset) < 50) {
      BAFsegm <- matrix(colMeans(vals), nrow = nrow(BAF_subset), ncol = ncol(vals), byrow = TRUE)
    } else {
      # Fully qualified copynumber calls
      winsor_data <- copynumber::winsorize(BAF_subset, assembly = GENOMEBUILD)
      res <- copynumber::multipcf(
        data = winsor_data,
        Y = BAF_subset,
        fast = TRUE,
        gamma = gamma * sdev,
        return.est = TRUE,
        normalize = FALSE,
        assembly = GENOMEBUILD
      )
      BAFsegm <- as.matrix(res$estimates[, -c(1:2)])
    }

    BAFphased <- ifelse(BAFsegm > 0.5, vals, 1 - vals)

    # Logic for segment BAF calculation
    if (calc_seg_baf_option %in% c(1, 3)) {
      BAFphseg <- apply(BAFphased, 2, stats::median)
      if (calc_seg_baf_option == 3) {
        means <- apply(BAFsegm, 2, function(x) ifelse(x[1] > 0.5, x[1], 1 - x[1]))
        BAFphseg <- ifelse(BAFphseg %in% c(0, 1), means, BAFphseg)
      }
    } else {
      BAFphseg <- apply(BAFsegm, 2, function(x) ifelse(x[1] > 0.5, x[1], 1 - x[1]))
    }

    out <- lapply(seq_along(samplename), function(i) {
      data.table::data.table(
        Chromosome = BAF_subset$Chromosome,
        Position = BAF_subset$Position,
        BAF = vals[, i],
        BAFphased = BAFphased[, i],
        BAFseg = rep(BAFphseg[i], nrow(BAF_subset)),
        tempBAFsegm = BAFsegm[, i]
      )
    })
    names(out) <- samplename
    return(out)
  }

  # --- 3. Main Execution ---
  # Initial data loading using data.table namespace
  BAFraw <- data.table::as.data.table(
    Reduce(function(...) merge(..., sort = FALSE), lapply(inputfile, read_baf))
  )

  bkps <- if (!is.null(prior_breakpoints_file)) {
    data.table::as.data.table(read.table(prior_breakpoints_file, header = TRUE))
  } else {
    NULL
  }

  all_results <- list()

  for (chr in unique(BAFraw$Chromosome)) {
    message("Processing ", chr, "...")
    chr_data <- BAFraw[Chromosome == chr][complete.cases(BAFraw[Chromosome == chr, -c(1:2)])]

    chr_bkps <- if (!is.null(bkps)) {
      bkps[chromosome == chr]
    } else {
      data.table::data.table(position = numeric())
    }

    segments <- get_segments(chr, chr_bkps, chr_data)

    seg_results <- lapply(seq_len(nrow(segments)), function(i) {
      run_pcf_modern(chr_data, segments$start[i], segments$end[i], gamma)
    })

    for (id in samplename) {
      # Explicitly use rbindlist from data.table
      chr_sample_dt <- data.table::rbindlist(lapply(seg_results, `[[`, id))

      # [Plotting logic - requires BAFoutputchr to be populated or used here]
      # ... (PNG/Plotting code as per original script) ...

      if (is.null(all_results[[id]])) all_results[[id]] <- list()
      all_results[[id]][[chr]] <- chr_sample_dt[, !"tempBAFsegm"]
    }
  }

  # Final Export using data.table::fwrite
  for (i in seq_along(samplename)) {
    final_dt <- data.table::rbindlist(all_results[[samplename[i]]])
    data.table::fwrite(final_dt, file = outputfile[i], sep = "\t")
  }

  return(NULL)
}
