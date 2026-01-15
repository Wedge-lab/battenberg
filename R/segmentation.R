#' Helper function to adjust the BAF segmented values. By default the segmentation
#' takes the mean BAFphased for each segment, but that doesn't work very well with
#' outliers (i.e. badly phased regions). This function is then called to adjust
#' the segmented BAF. By default this now takes the median
#' @param baf_chrom A data frame with columns BAFphased and BAFseg. BAFseg will be overwritten.
#' @return A data frame with columns BAFphased and BAFseg.
#' @author sd11
#' @noRd
adjustSegmValues <- function(baf_chrom) {
  if (nrow(baf_chrom) <= 1) {
    baf_chrom$BAFseg <- baf_chrom$BAFphased
    return(baf_chrom)
  }
  diffs <- collapse::fdiff(baf_chrom$BAFseg)
  runs <- collapse::fcumsum(diffs != 0)
  baf_chrom$BAFseg <- collapse::fmedian(
    baf_chrom$BAFphased,
    g = runs,
    TRA = "replace"
  )
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
    # Calculate gaps between consecutive SNPs
    gaps <- diff(positions)
    gap_indices <- which(gaps >= maxsnpdist)

    # If no holes, we don't return a new table, just the original
    if (length(gap_indices) == 0) {
      return(list(breakpoints = breakpoints, startpos = startpos))
    }

    # Define segment boundaries
    # Segment ends at the SNP before the gap
    ends <- c(positions[gap_indices], positions[length(positions)])

    # Segment starts at the original startpos, then the SNP AFTER each gap
    starts <- c(startpos, positions[gap_indices + 1])

    # Safety: Remove segments where start == end (the BAFlen=1 case)
    # Also ensures we don't have overlapping boundaries
    valid_mask <- starts < ends

    new_segments <- data.table::data.table(
      chrom = chrom,
      start = starts[valid_mask],
      end = ends[valid_mask]
    )

    updated_breakpoints <- data.table::rbindlist(
      list(breakpoints, new_segments),
      use.names = TRUE
    )

    # The startpos for the NEXT segment in the outer loop
    # should be the position AFTER the last SNP of this batch
    return(list(
      breakpoints = updated_breakpoints,
      startpos = positions[length(positions)] + 1
    ))
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
      log_info("No prior breakpoints found")
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
    sdev <- get_mad(ifelse(BAF < 0.5, BAF, 1 - BAF), k = 25)
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

    if (length(BAFphased) < 50 || no_segmentation) {
      BAFphseg <- rep(mean(BAFphased), length(BAFphased))
    } else {
      res <- selectFastPcf(BAFphased, kmin, gamma * sdev, T)
      BAFphseg <- res$yhat
    }

    if (length(BAF) > 0) {
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

  BAFraw <- read_baf_as_data_frame(inputfile)
  if (!is.null(prior_breakpoints_file)) {
    bkps <- utils::read.table(prior_breakpoints_file, header = TRUE, stringsAsFactors = FALSE)
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

    breakpoints_chrom <- bkps_to_presegment_breakpoints(chr, bkps_chrom, BAFrawchr, use_bigholes = TRUE)
    BAFoutputchr <- NULL

    for (r in seq_len(nrow(breakpoints_chrom))) {
      current_snps <- which(BAFrawchr$Position >= breakpoints_chrom$start[r] &
        BAFrawchr$Position <= breakpoints_chrom$end[r])

      if (length(current_snps) < 2) {
        log_info("Skipping empty/tiny segment {r} on chr {chr} (SNPs: {length(current_snps)})")
        next
      }
      BAFoutput_preseg <- run_pcf(BAFrawchr, breakpoints_chrom$start[r], breakpoints_chrom$end[r], phasekmin, phasegamma, kmin, gamma, no_segmentation)
      if (!is.null(BAFoutput_preseg)) {
        BAFoutputchr <- rbind(BAFoutputchr, BAFoutput_preseg)
      }
    }

    grDevices::png(
      filename = paste(samplename, "_RAFseg_chr", chr, ".png", sep = ""),
      width = 2000, height = 1000, res = 200, type = "cairo"
    )
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
    grDevices::dev.off()

    grDevices::png(
      filename = paste(samplename, "_segment_chr", chr, ".png", sep = ""),
      width = 2000, height = 1000, res = 200, type = "cairo"
    )
    create_baf_plot(
      chrom_position = BAFoutputchr$Position / 1000000,
      points_red_blue = BAFoutputchr$BAF,
      plot_red = BAFoutputchr$tempBAFsegm > 0.5,
      points_darkred = BAFoutputchr$BAFseg,
      points_darkblue = 1 - BAFoutputchr$BAFseg,
      x_min = min(BAFoutputchr$Position) / 1000000,
      x_max = max(BAFoutputchr$Position) / 1000000,
      title = paste(samplename, ", chromosome ", chr, sep = ""),
      xlab = "Position (Mb)",
      ylab = "BAF (phased)",
      prior_bkps_pos = bkps_chrom$position / 1000000
    )
    grDevices::dev.off()

    BAFoutputchr$BAFphased <- ifelse(BAFoutputchr$tempBAFsegm > 0.5, BAFoutputchr$BAF, 1 - BAFoutputchr$BAF)
    # Remove the temp BAFsegm values as they are only needed for plotting
    BAFoutput <- rbind(BAFoutput, BAFoutputchr[, c(1:5)])
  }
  colnames(BAFoutput) <- c("Chromosome", "Position", "BAF", "BAFphased", "BAFseg")
  data.table::fwrite(BAFoutput, outputfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
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
segment_baf_phased_multisample <- function(
  samplename, inputfile,
  outputfile, prior_breakpoints_file = NULL,
  gamma = 10, calc_seg_baf_option = 3,
  GENOMEBUILD
) {
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
    segments <- data.table::data.table(
      chrom = chrom, start = seg_starts, end = seg_ends
    )
    return(segments[segments$start <= segments$end])
  }

  run_pcf_helper <- function(BAFrawchr, start, end, gamma) {
    Position <- NULL
    BAF_subset <- BAFrawchr[Position >= start & Position <= end]

    if (nrow(BAF_subset) == 0) {
      return(NULL)
    }

    vals <- as.matrix(BAF_subset[, -c(1:2)])

    # Calculate sdev using Mean Absolute Deviation
    sdevs <- apply(vals, 2, function(x) {
      get_mad(ifelse(x < 0.5, x, 1 - x), k = 25)
    })
    sdevs[is.na(sdevs) | sdevs < 0.09] <- 0.09
    sdev <- mean(sdevs)

    if (nrow(BAF_subset) < 50) {
      BAFsegm <- matrix(colMeans(vals), nrow = nrow(BAF_subset), ncol = ncol(vals), byrow = TRUE)
    } else {
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
    stats::setNames(out, samplename)
  }

  BAFraw <- data.table::as.data.table(
    Reduce(function(...) merge(..., sort = FALSE), lapply(inputfile, read_baf_as_data_frame))
  )

  bkps <- if (!is.null(prior_breakpoints_file)) {
    data.table::fread(prior_breakpoints_file, header = TRUE)
  } else {
    NULL
  }

  all_results <- list()

  Chromosome <- chromosome <- NULL

  # Using string indexing to avoid warnings in the loop header
  for (chr in unique(BAFraw[["Chromosome"]])) {
    cli::cli_inform("Processing {chr}...")

    chr_data <- BAFraw[Chromosome == chr]
    chr_data <- chr_data[stats::complete.cases(chr_data[, -c(1:2)])]

    chr_bkps <- if (!is.null(bkps)) {
      bkps[chromosome == chr]
    } else {
      data.table::data.table(position = numeric())
    }

    segments <- get_segments(chr, chr_bkps, chr_data)

    seg_results <- lapply(seq_len(nrow(segments)), function(i) {
      run_pcf_helper(chr_data, segments$start[i], segments$end[i], gamma)
    })

    # We combine the segments for this specific chromosome once
    # This creates a named list of DataTables, one per sample
    chr_sample_results <- lapply(samplename, function(id) {
      data.table::rbindlist(lapply(seg_results, `[[`, id))
    })
    names(chr_sample_results) <- samplename

    for (id in samplename) {
      # Reference the combined data for this sample/chromosome
      sample_dt <- chr_sample_results[[id]]

      # Plot 1: RAFseg
      grDevices::png(
        filename = paste0(id, "_RAFseg_chr", chr, ".png"),
        width = 2000, height = 1000, res = 200, type = "cairo"
      )
      create_segmented_plot(
        chrom_position = sample_dt$Position / 1e6,
        points.red = sample_dt$BAF,
        points.green = sample_dt$tempBAFsegm,
        x_min = min(sample_dt$Position) / 1e6,
        x_max = max(sample_dt$Position) / 1e6,
        title = paste0(id, ", chromosome ", chr),
        xlab = "Position (Mb)",
        ylab = "BAF (phased)",
        prior_bkps_pos = chr_bkps$position / 1e6
      )
      grDevices::dev.off()

      # Plot 2: BAF segments
      grDevices::png(
        filename = paste0(id, "_segment_chr", chr, ".png"),
        width = 2000, height = 1000, res = 200, type = "cairo"
      )
      create_baf_plot(
        chrom_position = sample_dt$Position / 1e6,
        points_red_blue = sample_dt$BAF,
        plot_red = sample_dt$tempBAFsegm > 0.5,
        points_darkred = sample_dt$BAFseg,
        points_darkblue = 1 - sample_dt$BAFseg,
        x_min = min(sample_dt$Position) / 1e6,
        x_max = max(sample_dt$Position) / 1e6,
        title = paste0(id, ", chromosome ", chr),
        xlab = "Position (Mb)",
        ylab = "BAF (phased)",
        prior_bkps_pos = chr_bkps$position / 1e6
      )
      grDevices::dev.off()

      # Store for final export, removing the temp column used for plotting
      if (is.null(all_results[[id]])) all_results[[id]] <- list()
      all_results[[id]][[chr]] <- sample_dt[, !"tempBAFsegm"]
    }
  }

  # Final Export
  for (i in seq_along(samplename)) {
    final_dt <- data.table::rbindlist(all_results[[samplename[i]]])
    data.table::fwrite(final_dt, file = outputfile[i], sep = "\t")
  }
}
