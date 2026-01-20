#' Obtain BAF and LogR from the Cell line (tumour only) allele counts
#'
#' Function to generate BAF and LogR files based on allele counts of the Cell line.
#' It also generates the input data required by the following 'cell_line_reconstruct_normal' function.
#' @param TUMOURNAME The tumour name used for Battenberg (i.e. the cell line BAM file name without the '.bam' extension).
#' @param g1000alleles_prefix Prefix to where 1000 Genomes allele files can be found.
#' @param chrom_names A vector with allowed chromosome names.
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export

cell_line_baf_logR <- function(TUMOURNAME, g1000alleles_prefix, chrom_names) {
  # read heterozygous SNPs per chromosome for alleleCounter files & 1000G allele files####
  AC <- list() # alleleCounts
  AL <- list() # 1000G alleles
  MaC <- list() # matched alleleCounts
  OHET <- list() # HET SNP data

  for (chr in chrom_names) {
    # read in alleleCounter output for each chromosome (FAST)
    ac_file <- paste0(TUMOURNAME, "_alleleFrequencies_chr", chr, ".txt")
    if (!file.exists(ac_file) || file.size(ac_file) == 0) {
      log_failure("Allele count file '{ac_file}' is missing or empty. Preprocessing cannot continue.")
    }
    ac <- data.table::fread(ac_file, header = FALSE, stringsAsFactors = FALSE)
    if (nrow(ac) == 0) {
      log_failure("Allele count file '{ac_file}' contains no data.")
    }
    data.table::setorder(ac, V2)
    AC[[chr]] <- ac
    log_info("length(AC): '{length(AC)}'")

    # match allele counts with respective SNP alleles
    al_file <- paste0(g1000alleles_prefix, chr, ".txt")
    if (!file.exists(al_file) || file.size(al_file) == 0) {
      log_failure("1000G alleles file '{al_file}' is missing or empty.")
    }
    al <- data.table::fread(al_file, header = TRUE, stringsAsFactors = FALSE)
    if (nrow(al) == 0) {
      log_failure("1000G alleles file '{al_file}' contains no data.")
    }
    AL[[chr]] <- al
    log_info("length(AL): '{length(AL)}'")

    ref <- al$a0
    alt <- al$a1

    # Matrix indexing for lightning-fast extraction
    m_ac <- as.matrix(ac)
    REF <- m_ac[cbind(seq_len(nrow(al)), ref + 2)]
    ALT <- m_ac[cbind(seq_len(nrow(al)), alt + 2)]

    mac <- data.frame(ref = REF, alt = ALT)
    mac$depth <- as.numeric(mac$ref) + as.numeric(mac$alt)
    mac$baf <- as.numeric(mac$alt) / as.numeric(mac$depth)

    if (nrow(mac) == 0) {
      log_failure("No matching SNPs found between allele counts and 1000G alleles for chromosome {chr}.")
    }

    o <- cbind(al, mac)
    names(o) <- c("Position", "a0", "a1", "ref", "alt", "depth", "baf")
    MaC[[chr]] <- o

    # Extract HET SNPs
    ohet <- o[which(o$baf >= 0.10 & o$baf <= 0.90 & o$depth > 10), ]
    if (nrow(ohet) < 50) {
      log_warning("Extremely low heterozygosity detected on chromosome {chr} (n={nrow(ohet)}). Results may be unreliable.")
    }
    if (nrow(ohet) > 0) {
      ohet$Position2 <- c(
        ohet$Position[2:nrow(ohet)],
        2 * ohet$Position[nrow(ohet)] - ohet$Position[nrow(ohet) - 1]
      )
      ohet$Position_dist <- ohet$Position2 - ohet$Position
      ohet$Position_dist_percent <- ohet$Position_dist / max(ohet$Position_dist)
    }
    OHET[[chr]] <- ohet
    log_info("chromosome {chr} file read")
  }

  # CREATE mutantBAF and mutantLogR *.tab files #
  cellline <- TUMOURNAME

  # Assemble MAC efficiently (O(N))
  MAC_list <- lapply(chrom_names, function(chr) {
    data.frame(chr = chr, MaC[[chr]], stringsAsFactors = FALSE)
  })
  MAC <- collapse::rowbind(MAC_list)
  names(MAC) <- c("chr", "position", "a0", "a1", "ref", "alt", "coverage", "baf")

  log_info("Sync complete. dim(MAC): {paste(dim(MAC), collapse = ' ')}")

  # LogR calculation
  MAC$logr <- log2(MAC$coverage / mean(MAC$coverage, na.rm = TRUE))
  MACC <- MAC[which(!is.na(MAC$baf)), ]

  # Prepare and save BAF
  BAF <- data.frame(
    Chromosome = MACC$chr,
    Position = MACC$position,
    cellline = MACC$baf
  )
  names(BAF)[3] <- cellline
  # Standardization
  BAF$Chromosome[BAF$Chromosome %in% c("23", 23)] <- "X"
  data.table::setorder(BAF, Chromosome, Position)
  data.table::fwrite(BAF, paste0(cellline, "_mutantBAF.tab"), sep = "\t")
  rm(BAF)

  # Prepare and save LogR
  LogR_out <- data.frame(
    Chromosome = MACC$chr,
    Position = MACC$position,
    cellline = MACC$logr
  )
  names(LogR_out)[3] <- cellline
  LogR_out$Chromosome[LogR_out$Chromosome %in% c("23", 23)] <- "X"
  data.table::setorder(LogR_out, Chromosome, Position)
  data.table::fwrite(LogR_out, paste0(cellline, "_mutantLogR.tab"), sep = "\t")

  return(list(
    OHET = OHET,
    AL   = AL,
    AC   = AC,
    LogR = LogR_out
  ))
}

#' Reconstruct normal-pair allele count files for cell lines
#'
#' Function to generate normal-pair allele count files based on IVD-PCF and inter-hetSNP logR-based LOH detection (IVD: Inter-Variant Distance, het: heterozygote)
#' This method reconstructs the normal-pair counts by using the allele counts of the Cell line as template.
#' It fills the detected LOH regions with evenly-distributed hetSNPs with the density estimated based on each chromosome in each tumour sample.
#' It essentially informs Battenberg of the location of hetSNPs across the genome in the tumour sample.
#' @param TUMOURNAME The tumour name used for Battenberg (i.e. the cell line BAM file name without the '.bam' extension).
#' @param NORMALNAME The normal name used for naming the generated normal-pair allele counts files.
#' @param chrom_coord Full path to the file with chromosome coordinates including start, end and left/right centromere positions
#' @param chrom Chromosome number for which normal-pair will be reconstructed (1,2, etc.)
#' @param CL_OHET List of observed heterozygous SNPs across all chromosomes generated within the cell_line_baf_logR function
#' @param CL_AL List of alleles at SNPs across all chromosomes generated within the cell_line_baf_logR function
#' @param CL_AC List of allele counts at SNPs across all chromosomes generated within the cell_line_baf_logR function
#' @param CL_LogR Dataframe of genomewide LogR values for SNPs across all chromosomes generated within the cell_line_baf_logR function
#' @param GAMMA_IVD The PCF gamma value for segmentation of 1000G hetSNP IVD values (Default 1e5).
#' @param KMIN_IVD The min number of SNPs to support a segment in PCF of 1000G hetSNP IVD values (Default 50)
#' @param CENTROMERE_DIST The minimum distance from the centromere to ignore in analysis due to the noisy nature of data in the vicinity of centromeres (Default 5e5)
#' @param CENTROMERE_NOISE_SEG_SIZE The maximum size of PCF segment to be removed as noise when it overlaps with the centromere due to the noisy nature of data (Default 1e6)
#' @param MIN_HET_DIST The minimum distance for detecting higher resolution inter-hetSNP regions with potential LOH while accounting for inherent homozygote stretches (Default 1e5)
#' @param GAMMA_LOGR The PCF gamma value for confirming LOH within each inter-hetSNP candidate segment (Default 100)
#' @param LENGTH_ADJACENT The length of adjacent regions either side of a candidate inter-hetSNP LOH region to be plotted (Default 5e4)
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export

cell_line_reconstruct_normal <- function(
  TUMOURNAME, NORMALNAME,
  chrom_coord, chrom,
  CL_OHET, CL_AL,
  CL_AC, CL_LogR,
  GAMMA_IVD, KMIN_IVD,
  CENTROMERE_NOISE_SEG_SIZE,
  CENTROMERE_DIST, MIN_HET_DIST,
  GAMMA_LOGR, LENGTH_ADJACENT
) {
  # IDENTIFY REGIONS OF LOH ####
  colClasses <- c(chr = "numeric", start = "numeric", cen.left.base = "numeric", cen.right.base = "numeric", end = "numeric")
  # Use fast I/O
  chr_loc <- data.table::fread(chrom_coord, colClasses = colClasses, header = TRUE, stringsAsFactors = FALSE)
  data.table::setDF(chr_loc)
  chr_loc$length <- (chr_loc$cen.left.base - chr_loc$start) + (chr_loc$end - chr_loc$cen.right.base)

  # identify LOH by IVD-PCF
  LOH <- list()
  PCF_folder <- "PCF_plots"
  if (!dir.exists(PCF_folder)) {
    dir.create(PCF_folder)
  }
  i <- chrom
  log_info("chrom={i}")
  pcf_input <- data.frame(chr = i, position = CL_OHET[[i]]$Position, IVD = (CL_OHET[[i]]$Position_dist_percent))
  pcf_input <- pcf_input[which(pcf_input$position < chr_loc[i, "cen.left.base"] - CENTROMERE_DIST | pcf_input$position > chr_loc[i, "cen.right.base"] + CENTROMERE_DIST), ]
  pcf_input <- pcf_input[which(pcf_input$position >= chr_loc[i, "start"] & pcf_input$position <= chr_loc[i, "end"]), ] # use only regions covered with gcCorrect LogR range
  PCF <- copynumber::pcf(pcf_input, gamma = GAMMA_IVD, kmin = KMIN_IVD)
  grDevices::pdf(paste0(PCF_folder, "/", TUMOURNAME, "_chr", i, "_PCF_plot.pdf"))
  copynumber::plotChrom(pcf_input, PCF)
  grDevices::dev.off()
  PCF$diff <- PCF$end.pos - PCF$start.pos

  # Decide if there is any LOH based on PCF and chr_snp_density
  # density of HET SNPs across the region covered by HET SNPs
  chr_snp_density <- nrow(pcf_input) / (pcf_input$position[nrow(pcf_input)] - pcf_input$position[1])
  min_normal_snp_density <- 0.0001
  # LOH regions
  loh_regions <- PCF[which(round(PCF$mean, 3) > 0.001), ]
  # only keep segments with minimum of 2 probes (SNPs) in PCF jump
  loh_regions <- loh_regions[which(loh_regions$n.probes > 1), ]
  if (nrow(loh_regions) > 0) {
    # can change chr_snp_density from 0.00005 to 0.0001 as conservative measure - done
    if (mean(pcf_input$IVD) > 0.01 && chr_snp_density < min_normal_snp_density) {
      # mean(pcf_input$IVD) or mean(PCF$mean) indicates presence of jumps in IVD
      loh_regions <- loh_regions # LOH regions
      log_info("full-length chromosomal loss at chr {i}")
    } else if (sum(loh_regions$diff) >= ((pcf_input$position[nrow(pcf_input)] - pcf_input$position[1])) * 0.9 && chr_snp_density > min_normal_snp_density) {
      # do PCF regions cover >=90% of the chromosome & is the chromosome snp density above the minimum
      loh_regions <- 0 # LOH regions
      log_info("no PCF jumps at chr {i}")
    } else {
      loh_regions <- loh_regions # LOH regions
      log_info("likely partial LOH(s) at chr {i}")
    }
  } else {
    loh_regions <- 0
  }

  # filter regions for those next to the centromere and 'short'
  noise <- NULL
  if (!is.null(nrow(loh_regions))) {
    for (j in seq_len(nrow(loh_regions))) {
      if (loh_regions$arm[j] == "p") {
        if (loh_regions$end.pos[j] > chr_loc$cen.left.base[i] && loh_regions$diff[j] < CENTROMERE_NOISE_SEG_SIZE) {
          # FOR EXCLUSION: segment is short IVD region (default<1Mb) and endpos is over the p-arm limit (ending point)
          noise <- c(noise, j)
        }
      }
      if (loh_regions$arm[j] == "q") {
        if (loh_regions$start.pos[j] < chr_loc$cen.right.base[i] && loh_regions$diff[j] < CENTROMERE_NOISE_SEG_SIZE) {
          # FOR EXCLUSION: segment is short IVD region (default<1Mb) and startpos is below the q-arm limit (starting point)
          noise <- c(noise, j)
        }
        if (loh_regions$start.pos[j] < (chr_loc$cen.right.base[i] + 1e5) && loh_regions$diff[j] > CENTROMERE_NOISE_SEG_SIZE && !is.na(match(chrom, c(1, 9, 16)))) {
          # qARM of Chr 1,9,16 have large heterochromatin region next to centromere + 100kb tolerance for start of heterochromatin region
          noise <- c(noise, j)
        }
      }
    }
  } else {
    log_info("no 'centromere noise' calculation")
  }
  if (!is.null(noise)) {
    LOH_regions <- loh_regions[-noise, ]
  } else {
    LOH_regions <- loh_regions
  }
  # remove LOH regions in the p arm of acrocentric chromosomes 13,14,15,21 and 22
  if (!is.na(match(i, c(13:15, 21:22))) && !is.null(nrow(LOH_regions))) {
    LOH_regions <- LOH_regions[which(LOH_regions$arm != "p"), ]
  }
  # remove LOH regions which do not have negative LogR and are essentially stretches of homozygosity
  if (!is.null(nrow(LOH_regions)) && nrow(LOH_regions) > 0) {
    logr <- CL_LogR[which(CL_LogR$Chromosome == i), ]
    colnames(logr)[3] <- "LogR"
    logr$Position <- as.numeric(logr$Position)

    # Use findInterval for O(M) mapping to segments
    snp_to_loh <- findInterval(logr$Position, LOH_regions$start.pos)
    valid_mask <- snp_to_loh > 0 & logr$Position <= LOH_regions$end.pos[pmax(1, snp_to_loh)]

    if (any(valid_mask)) {
      stats <- collapse::fgroup_by(logr[valid_mask, ], snp_to_loh[valid_mask]) |>
        collapse::fsummarise(medcov = fmedian(LogR), cov = fmean(LogR), n = fnobs(LogR))

      # Only keep regions that meet the LOH criteria
      keep_regions <- stats$g[stats$cov < -0.8 & stats$medcov < -0.8 & stats$n >= 10]
      if (length(keep_regions) > 0) {
        LOH_regions <- LOH_regions[keep_regions, ]
      } else {
        LOH_regions <- data.frame()
      }
    } else {
      LOH_regions <- data.frame()
    }
  }
  if (is.null(dim(LOH_regions))) {
    log_info("no LOH detected in chr {i}")
    LOH[[i]] <- 0
  } else if (dim(LOH_regions)[1] != 0 && dim(LOH_regions)[2] != 0) {
    log_info("we have LOH for {sum(LOH_regions$diff)} bp in chr {i}")
    LOH[[i]] <- data.frame(chr = i, LOH_regions)
  } else if (dim(LOH_regions)[1] == 0) {
    log_info("no LOH regions remained after noise correction for chr {i}")
    LOH[[i]] <- 0
  } else {
    log_info("unkown issue!")
  }
  log_info("chrom={i} IVD-PCF finished")

  # STEP 2 - get higher resolution LOH regions
  log_info("chrom={i}")
  # Use list for efficient non_LOH construction
  ac <- CL_AC[[i]]
  al <- CL_AL[[i]]
  names(ac) <- c("chr", "position", "A", "C", "G", "T", "depth")

  chr_interval <- c(chr_loc[i, "start"], chr_loc[i, "end"])
  if (!is.null(nrow(LOH[[i]])) && nrow(LOH[[i]]) > 0) {
    non_LOH_list <- list()
    for (j in 1:(nrow(LOH[[i]]) + 1)) {
      if (j == 1 && chr_interval[1] >= LOH[[i]]$start.pos[j]) {} else if (j == 1) {
        non_LOH_list[[length(non_LOH_list) + 1]] <- data.frame(start = chr_interval[1], end = LOH[[i]]$start.pos[j] - 1)
      } else if (j <= nrow(LOH[[i]]) && LOH[[i]]$arm[j] == LOH[[i]]$arm[j - 1]) {
        non_LOH_list[[length(non_LOH_list) + 1]] <- data.frame(start = LOH[[i]]$end.pos[j - 1] + 1, end = LOH[[i]]$start.pos[j] - 1)
      } else if (j <= nrow(LOH[[i]])) {
        non_LOH_list[[length(non_LOH_list) + 1]] <- data.frame(
          start = c(LOH[[i]]$end.pos[j - 1] + 1, chr_loc[i, ]$cen.right.base),
          end = c(chr_loc[i, ]$cen.left.base, LOH[[i]]$start.pos[j] - 1)
        )
      } else if ((LOH[[i]]$end.pos[j - 1] + 1) < chr_interval[2]) {
        non_LOH_list[[length(non_LOH_list) + 1]] <- data.frame(start = LOH[[i]]$end.pos[j - 1] + 1, end = chr_interval[2])
      }
    }
    non_LOH <- collapse::rowbind(non_LOH_list)
  } else {
    non_LOH <- data.frame(start = chr_interval[1], end = chr_interval[2])
  }

  if (nrow(non_LOH) > 0) {
    # Check for centromere crossing and split if necessary
    cross_idx <- which(non_LOH$start < chr_loc[i, ]$cen.left.base & non_LOH$end > chr_loc[i, ]$cen.right.base)
    if (length(cross_idx) > 0) {
      to_split <- non_LOH[cross_idx, ]
      non_LOH <- non_LOH[-cross_idx, ]
      split_list <- list(
        non_LOH,
        data.frame(start = to_split$start, end = chr_loc[i, ]$cen.left.base),
        data.frame(start = chr_loc[i, ]$cen.right.base, end = to_split$end)
      )
      non_LOH <- collapse::rowbind(split_list)
    }
    non_LOH$diff <- non_LOH$end - non_LOH$start
    non_LOH <- non_LOH[non_LOH$diff > 0, ]
  }

  non_LOH <- non_LOH[order(non_LOH$start), ]

  # identify LOH by inter-het regions
  ohet <- CL_OHET[[i]]
  nSNPs <- as.numeric(nrow(CL_LogR))
  logr <- CL_LogR[which(CL_LogR$Chromosome == i), ]
  colnames(logr)[3] <- "LogR"
  logr$Position <- as.numeric(logr$Position)

  pLOH_collector_list <- list() # to collect results of p-arm analysis
  if (!is.null(non_LOH)) {
    if (is.na(match(i, c(13, 14, 15, 21, 22)))) {
      log_info("START {i}, p ARM")
      PARM <- non_LOH[which(non_LOH$end <= chr_loc[i, ]$cen.left.base), ]
      if (nrow(PARM) > 0) {
        parm <- PARM
      } else if (nrow(PARM) == 0 && sum(non_LOH$diff) != 0) {
        parm <- data.frame(start = chr_interval[1], end = chr_loc[i, ]$cen.left.base - CENTROMERE_DIST)
      } else {
        log_info("unknown issue")
      }

      if (parm[nrow(parm), 1] < (parm[nrow(parm), 2] - CENTROMERE_DIST)) {
        # to exclude the last CENTROMERE_DIST segment next to the centromere (left side) - too noisy
        parm[nrow(parm), 2] <- parm[nrow(parm), 2] - CENTROMERE_DIST
      } else {
        parm <- parm[-nrow(parm), ]
      }
      #
      for (seg in seq_len(nrow(parm))) {
        LoH_list <- list()
        # IVD-based breakpoints for small regions#
        seg_ivd <- ohet[which(ohet$Position_dist >= MIN_HET_DIST & ohet$Position >= parm$start[seg] & ohet$Position <= parm$end[seg]), ]
        if (nrow(seg_ivd) > 0) {
          logr_in_seg_idx <- which(logr$Position >= parm$start[seg] & logr$Position <= parm$end[seg])
          if (length(logr_in_seg_idx) > 0) {
            logr_seg <- logr[logr_in_seg_idx, ]
            starts_idx <- findInterval(seg_ivd$Position, logr_seg$Position) + 1
            ends_idx <- findInterval(seg_ivd$Position + seg_ivd$Position_dist, logr_seg$Position)

            for (j in seq_len(nrow(seg_ivd))) {
              if (starts_idx[j] > ends_idx[j]) next
              COV <- logr_seg[starts_idx[j]:ends_idx[j], ]
              medcov <- collapse::fmedian(COV$LogR)
              cov <- mean(COV$LogR)
              denSNP <- nrow(COV) / (nSNPs / sum(chr_loc$length) * seg_ivd$Position_dist[j])

              if (!is.na(cov) && cov < -0.8 && medcov < -0.8 && denSNP > 0.5) {
                jpcf <- copynumber::pcf(COV, gamma = GAMMA_LOGR, verbose = FALSE)
                jpcf_loh <- jpcf[which(jpcf$mean < -0.8), ]
                if (nrow(jpcf_loh) > 0) {
                  loh <- data.frame(
                    start = jpcf_loh$start.pos[1],
                    end = jpcf_loh$end.pos[nrow(jpcf_loh)],
                    LogR = mean(jpcf_loh$mean),
                    denSNP = denSNP,
                    stringsAsFactors = FALSE
                  )
                  loh$N <- sum(COV$Position >= loh$start & COV$Position <= loh$end)
                  if (loh$N >= 10) LoH_list[[length(LoH_list) + 1]] <- loh
                }
              }
            }
          }
        }

        if (length(LoH_list) > 0) {
          LoH <- collapse::rowbind(LoH_list)
          LoH_regions_list <- list()
          start <- LoH$start[1]
          end <- LoH$end[1]
          if (nrow(LoH) > 1) {
            for (j in 2:nrow(LoH)) {
              if (LoH$start[j] <= end) {
                end <- max(end, LoH$end[j])
              } else {
                LoH_regions_list[[length(LoH_regions_list) + 1]] <- data.frame(chrom = i, arm = "p", start.pos = start, end.pos = end)
                start <- LoH$start[j]
                end <- LoH$end[j]
              }
            }
          }
          LoH_regions_list[[length(LoH_regions_list) + 1]] <- data.frame(chrom = i, arm = "p", start.pos = start, end.pos = end)
          pLOH_collector_list[[length(pLOH_collector_list) + 1]] <- collapse::rowbind(LoH_regions_list)
        }
      }
      pLOH_regions <- collapse::rowbind(pLOH_collector_list)

      if (nrow(pLOH_regions) > 0) {
        grDevices::pdf(paste0(TUMOURNAME, "_chr", i, "_", MIN_HET_DIST / 1e3, "k_based_pLOH_events.pdf"))
        suppressWarnings(
          for (s in seq_len(nrow(pLOH_regions))) {
            sBAF <- ggplot2::ggplot(ohet, ggplot2::aes(rlang::.data$Position, rlang::.data$baf)) +
              ggplot2::geom_jitter() +
              ggplot2::ylim(0, 1) +
              ggplot2::geom_vline(xintercept = c(pLOH_regions$start.pos[s], pLOH_regions$end.pos[s]), col = "red", linetype = "longdash") +
              ggplot2::xlim(pLOH_regions$start.pos[s] - LENGTH_ADJACENT, pLOH_regions$end.pos[s] + LENGTH_ADJACENT) +
              ggplot2::ggtitle(paste("pARM LOH region", s)) +
              ggplot2::labs(y = "BAF")
            sLogR <- ggplot2::ggplot(logr, ggplot2::aes(rlang::.data$Position, rlang::.data$LogR)) +
              ggplot2::geom_jitter() +
              ggplot2::ylim(-5.2, 1.2) +
              ggplot2::geom_vline(xintercept = c(pLOH_regions$start.pos[s], pLOH_regions$end.pos[s]), col = "red", linetype = "longdash") +
              ggplot2::xlim(pLOH_regions$start.pos[s] - LENGTH_ADJACENT, pLOH_regions$end.pos[s] + LENGTH_ADJACENT)
            grid::grid.newpage()
            grid::grid.draw(rbind(ggplot2::ggplotGrob(sBAF), ggplot2::ggplotGrob(sLogR), size = "last"))
          }
        )
        grDevices::dev.off()
        #
        log_info("Candidate LOH regions plotted for pARM")
      }
    } else {
      pLOH_regions <- data.frame() # ensure it exists
      log_info("chr {i} is acrocentric - no p arm analysis")
    }
    # Q ARM RUN:
    log_info("START {i} q ARM")
    qLOH_collector_list <- list()
    QARM <- non_LOH[which(non_LOH$start >= chr_loc[i, ]$cen.right.base), ]
    if (nrow(QARM) > 0) {
      qarm <- QARM
    } else if (nrow(QARM) == 0 && sum(non_LOH$diff) != 0) {
      qarm <- data.frame(start = chr_loc[i, ]$cen.right.base, end = chr_interval[2])
    } else {
      log_info("unknown issue")
    }
    # to exclude the first CENTROMERE_DIST segment next to the centromere (right side) - noisy
    qarm[1, 1] <- qarm[1, 1] + CENTROMERE_DIST
    qarm$diff <- qarm$end - qarm$start
    #
    # search per non_LOH segment
    for (seg in seq_len(nrow(qarm))) {
      LoH_list <- list()
      # IVD-based breakpoints for small regions#
      seg_ivd <- ohet[which(ohet$Position_dist >= MIN_HET_DIST & ohet$Position >= qarm$start[seg] & ohet$Position <= qarm$end[seg]), ]
      if (nrow(seg_ivd) > 0) {
        logr_in_seg_idx <- which(logr$Position >= qarm$start[seg] & logr$Position <= qarm$end[seg])
        if (length(logr_in_seg_idx) > 0) {
          logr_seg <- logr[logr_in_seg_idx, ]
          starts_idx <- findInterval(seg_ivd$Position, logr_seg$Position) + 1
          ends_idx <- findInterval(seg_ivd$Position + seg_ivd$Position_dist, logr_seg$Position)

          for (j in seq_len(nrow(seg_ivd))) {
            if (starts_idx[j] > ends_idx[j]) next
            COV <- logr_seg[starts_idx[j]:ends_idx[j], ]
            cov <- mean(COV$LogR)
            medcov <- collapse::fmedian(COV$LogR)
            denSNP <- nrow(COV) / (nSNPs / sum(chr_loc$length) * seg_ivd$Position_dist[j])
            if (!is.na(cov) && cov < -0.8 && medcov < -0.8 && denSNP > 0.5) {
              jpcf <- copynumber::pcf(COV, gamma = GAMMA_LOGR, verbose = FALSE)
              jpcf_loh <- jpcf[which(jpcf$mean < -0.8), ]
              if (nrow(jpcf_loh) > 0) {
                loh <- data.frame(
                  start = jpcf_loh$start.pos[1],
                  end = jpcf_loh$end.pos[nrow(jpcf_loh)],
                  LogR = mean(jpcf_loh$mean),
                  denSNP = denSNP,
                  stringsAsFactors = FALSE
                )
                loh$N <- sum(COV$Position >= loh$start & COV$Position <= loh$end)
                if (loh$N >= 10) LoH_list[[length(LoH_list) + 1]] <- loh
              }
            }
          }
        }
      }

      if (length(LoH_list) > 0) {
        LoH <- collapse::rowbind(LoH_list)
        LoH_regions_list <- list()
        start <- LoH$start[1]
        end <- LoH$end[1]
        if (nrow(LoH) > 1) {
          for (j in 2:nrow(LoH)) {
            if (LoH$start[j] <= end) {
              end <- max(end, LoH$end[j])
            } else {
              LoH_regions_list[[length(LoH_regions_list) + 1]] <- data.frame(chrom = i, arm = "q", start.pos = start, end.pos = end)
              start <- LoH$start[j]
              end <- LoH$end[j]
            }
          }
        }
        LoH_regions_list[[length(LoH_regions_list) + 1]] <- data.frame(chrom = i, arm = "q", start.pos = start, end.pos = end)
        qLOH_collector_list[[length(qLOH_collector_list) + 1]] <- collapse::rowbind(LoH_regions_list)
      }
    }
    qLOH_regions <- collapse::rowbind(qLOH_collector_list)

    if (nrow(qLOH_regions) > 0) {
      grDevices::pdf(paste0(TUMOURNAME, "_chr", i, "_", MIN_HET_DIST / 1e3, "k_based_qLOH_events.pdf"))
      suppressWarnings(
        for (s in seq_len(nrow(qLOH_regions))) {
          sBAF <- ggplot2::ggplot(ohet, ggplot2::aes(rlang::.data$Position, rlang::.data$baf)) +
            ggplot2::geom_jitter() +
            ggplot2::ylim(0, 1) +
            ggplot2::geom_vline(xintercept = c(qLOH_regions$start.pos[s], qLOH_regions$end.pos[s]), col = "red", linetype = "longdash") +
            ggplot2::xlim(qLOH_regions$start.pos[s] - LENGTH_ADJACENT, qLOH_regions$end.pos[s] + LENGTH_ADJACENT) +
            ggplot2::ggtitle(paste("qARM LOH region", s))
          sLogR <- ggplot2::ggplot(logr, ggplot2::aes(rlang::.data$Position, rlang::.data$LogR)) +
            ggplot2::geom_jitter() +
            ggplot2::ylim(-5.2, 1.2) +
            ggplot2::geom_vline(xintercept = c(qLOH_regions$start.pos[s], qLOH_regions$end.pos[s]), col = "red", linetype = "longdash") +
            ggplot2::xlim(qLOH_regions$start.pos[s] - LENGTH_ADJACENT, qLOH_regions$end.pos[s] + LENGTH_ADJACENT)
          grid::grid.newpage()
          grid::grid.draw(rbind(ggplot2::ggplotGrob(sBAF), ggplot2::ggplotGrob(sLogR), size = "last"))
        }
      )
      grDevices::dev.off()
      #
      log_info("Candidate LOH regions plotted for qARM")
    }
    # merge LOH regions of both methods
    LOH_merge_list <- list()
    if (!is.null(pLOH_regions) && nrow(pLOH_regions) > 0) {
      LOH_merge_list[[length(LOH_merge_list) + 1]] <- pLOH_regions
    }
    if (!is.null(qLOH_regions) && nrow(qLOH_regions) > 0) {
      LOH_merge_list[[length(LOH_merge_list) + 1]] <- qLOH_regions
    }

    if (length(LOH_merge_list) > 0) {
      LOH_regions_final <- collapse::rowbind(LOH_merge_list)
      if (!is.null(LOH[[i]]) && !is.null(nrow(LOH[[i]])) && nrow(LOH[[i]]) > 0) {
        LOH[[i]] <- collapse::rowbind(LOH[[i]][, c("chrom", "arm", "start.pos", "end.pos")], LOH_regions_final)
        LOH[[i]] <- LOH[[i]][order(LOH[[i]]$start.pos), ]
      } else {
        LOH[[i]] <- LOH_regions_final
      }
    }

    # combine adjacent regions into larger regions of LOH
    if (!is.null(LOH[[i]]) && !is.null(nrow(LOH[[i]])) && nrow(LOH[[i]]) > 0) {
      LOH[[i]] <- LOH[[i]][!duplicated(LOH[[i]]), ]
      LOHall_list <- list()
      ChrArms <- unique(LOH[[i]]$arm)
      for (arm in ChrArms) {
        LOHarm <- LOH[[i]][LOH[[i]]$arm == arm, ]
        if (nrow(LOHarm) > 1) {
          start <- LOHarm$start.pos[1]
          end <- LOHarm$end.pos[1]
          for (j in 2:nrow(LOHarm)) {
            if (LOHarm$start.pos[j] <= end) {
              end <- max(end, LOHarm$end.pos[j])
            } else {
              LOHall_list[[length(LOHall_list) + 1]] <- data.frame(chrom = i, arm = arm, start.pos = start, end.pos = end)
              start <- LOHarm$start.pos[j]
              end <- LOHarm$end.pos[j]
            }
          }
          LOHall_list[[length(LOHall_list) + 1]] <- data.frame(chrom = i, arm = arm, start.pos = start, end.pos = end)
        } else {
          LOHall_list[[length(LOHall_list) + 1]] <- LOHarm[, c("chrom", "arm", "start.pos", "end.pos")]
        }
      }
      LOHall <- collapse::rowbind(LOHall_list)
    } else {
      LOHall <- LOH[[i]]
    }
  } else {
    # no non_LOH region was found - all chromosome is called as LOH
    if (!is.null(LOH[[i]]) && !is.null(nrow(LOH[[i]])) && nrow(LOH[[i]]) > 0) {
      LOHall <- LOH[[i]][, c("chrom", "arm", "start.pos", "end.pos")]
    } else {
      LOHall <- NULL
    }
  }

  if (!is.null(LOHall) && !is.null(nrow(LOHall)) && nrow(LOHall) > 0) {
    LOHall <- LOHall[!duplicated(LOHall), ]
    LOHall$diff <- LOHall$end.pos - LOHall$start.pos
  }

  # RECONSTRUCT alleleCounter files for the pseudo-NORMAL sample
  if (!is.null(LOHall) && !is.null(nrow(LOHall)) && nrow(LOHall) > 0) {
    names(ac) <- c("chr", "position", "A", "C", "G", "T", "depth")
    chr_interval <- c(ac$position[1], ac$position[nrow(ac)])

    # Get non_LOH regions based on LOHall
    non_LOH_list <- list()
    for (j in 1:(nrow(LOHall) + 1)) {
      if (j == 1 && chr_interval[1] >= LOHall$start.pos[j]) {} else if (j == 1) {
        non_LOH_list[[length(non_LOH_list) + 1]] <- data.frame(start = chr_interval[1], end = LOHall$start.pos[j] - 1)
      } else if (j <= nrow(LOHall) && LOHall$arm[j] == LOHall$arm[j - 1]) {
        non_LOH_list[[length(non_LOH_list) + 1]] <- data.frame(start = LOHall$end.pos[j - 1] + 1, end = LOHall$start.pos[j] - 1)
      } else if (j <= nrow(LOHall)) {
        non_LOH_list[[length(non_LOH_list) + 1]] <- data.frame(
          start = c(min(LOHall$end.pos[j - 1] + 1, chr_loc[i, ]$cen.left.base), chr_loc[i, ]$cen.right.base),
          end = c(chr_loc[i, ]$cen.left.base, LOHall$start.pos[j] - 1)
        )
      } else if ((LOHall$end.pos[j - 1] + 1) < chr_interval[2]) {
        non_LOH_list[[length(non_LOH_list) + 1]] <- data.frame(start = LOHall$end.pos[j - 1] + 1, end = chr_interval[2])
      }
    }
    non_LOH <- collapse::rowbind(non_LOH_list)
    non_LOH <- non_LOH[non_LOH$end >= non_LOH$start, ]

    if (nrow(non_LOH) > 0) {
      non_LOH$length <- non_LOH$end - non_LOH$start
      non_LOH_length <- sum(non_LOH$length)
      if (non_LOH_length > 1e6) {
        SNP_interval <- non_LOH_length / max(1, nrow(CL_OHET[[i]]))
      } else {
        SNP_interval <- 2000
      }
    } else {
      SNP_interval <- 2000
    }

    # Spike in heterozygotes in LOH regions
    lohs_list <- list()
    for (j in seq_len(nrow(LOHall))) {
      loh_idx <- which(ac$position >= LOHall$start.pos[j] & ac$position <= LOHall$end.pos[j])
      if (length(loh_idx) == 0) next
      loh <- ac[loh_idx, ]

      # Merge with alleles
      m <- merge(loh, al, by = "position")

      hetSNP_number <- max(floor(LOHall$diff[j] / SNP_interval), 10)
      if (nrow(m) >= hetSNP_number) {
        spike <- unique(c(1, floor(seq(1, nrow(m), length.out = hetSNP_number)), nrow(m)))
        for (k in spike) {
          m$depth[k] <- max(m$depth[k], 10)
          a0_col <- match(as.character(m$a0[k]), c("1", "2", "3", "4")) + 2
          a1_col <- match(as.character(m$a1[k]), c("1", "2", "3", "4")) + 2
          if (!is.na(a0_col)) m[k, a0_col] <- ceiling(m$depth[k] / 2)
          if (!is.na(a1_col)) m[k, a1_col] <- floor(m$depth[k] / 2)
        }
      } else {
        for (k in seq_len(nrow(m))) {
          m$depth[k] <- max(m$depth[k], 10)
          a0_col <- match(as.character(m$a0[k]), c("1", "2", "3", "4")) + 2
          a1_col <- match(as.character(m$a1[k]), c("1", "2", "3", "4")) + 2
          if (!is.na(a0_col)) m[k, a0_col] <- ceiling(m$depth[k] / 2)
          if (!is.na(a1_col)) m[k, a1_col] <- floor(m$depth[k] / 2)
        }
      }
      # Reorder columns to match ac
      lohs_list[[j]] <- m[, c("chr", "position", "A", "C", "G", "T", "depth")]
    }
    lohs <- collapse::rowbind(lohs_list)

    # Combine non_LOH regions
    non_lohs_list <- list()
    if (nrow(non_LOH) > 0) {
      for (j in seq_len(nrow(non_LOH))) {
        non_lohs_list[[j]] <- ac[ac$position >= non_LOH$start[j] & ac$position <= non_LOH$end[j], ]
      }
    }
    non_lohs <- collapse::rowbind(non_lohs_list)

    # Final assembly
    ac_out_list <- list(non_lohs, lohs)
    covered_pos <- c(lohs$position, non_lohs$position)
    missing_ac <- ac[!(ac$position %in% covered_pos), ]
    if (nrow(missing_ac) > 0) {
      ac_out_list[[3]] <- missing_ac
    }

    ac_out <- collapse::rowbind(ac_out_list)
    ac_out <- ac_out[order(ac_out$position), ]
    ac_out <- ac_out[!duplicated(ac_out$position), ]

    data.table::fwrite(ac_out, paste0(NORMALNAME, "_alleleFrequencies_chr", i, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    log_info("reconstruction OK - new alleleCounts file generated for chr {i}")
  } else {
    # No LOH identified
    data.table::fwrite(ac, paste0(NORMALNAME, "_alleleFrequencies_chr", i, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    log_info("No change to allele frequencies for chr {i}")
  }
}

#' Prepare WGS data of cell line for haplotype construction
#'
#' This function performs part of the Battenberg WGS pipeline: Counting alleles, generating BAF and logR,
#' reconstructing normal-pair allele counts for the cell line and performing GC content correction.
#'
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param tumourbam Full path to the tumour BAM file
#' @param tumourname Identifier to be used for tumour output files (i.e. the cell line BAM file name without the '.bam' extension).
#' @param chrom_coord Path to the chromosome coordinates file
#' @param g1000lociprefix Prefix path to the 1000 Genomes loci reference files
#' @param g1000allelesprefix Prefix path to the 1000 Genomes SNP allele reference files
#' @param gamma_ivd The PCF gamma value for segmentation of 1000G hetSNP IVD values (Default 1e5).
#' @param kmin_ivd The min number of SNPs to support a segment in PCF of 1000G hetSNP IVD values (Default 50)
#' @param centromere_noise_seg_size The maximum size of PCF segment to be removed as noise when it overlaps with the centromere due to the noisy nature of data (Default 1e6)
#' @param centromere_dist The minimum distance from the centromere to ignore in analysis due to the noisy nature of data in the vicinity of centromeres (Default 5e5)
#' @param min_het_dist The minimum distance for detecting higher resolution inter-hetSNP regions with potential LOH while accounting for inherent homozygote stretches (Default 1e5)
#' @param gamma_logr The PCF gamma value for confirming LOH within each inter-hetSNP candidate segment (Default 100)
#' @param length_adjacent The length of adjacent regions either side of a candidate inter-hetSNP LOH region to be plotted (Default 5e4)
#' @param gccorrectprefix Prefix path to GC content reference data
#' @param repliccorrectprefix Prefix path to replication timing reference data (supply NULL if no replication timing correction is to be applied)
#' @param min_base_qual Minimum base quality required for a read to be counted
#' @param min_map_qual Minimum mapping quality required for a read to be counted
#' @param allele_counts_dir Directory containing the allele counts files
#' @param min_normal_depth Minimum depth required in the normal for a SNP to be included
#' @param libs Path to the R libraries to be used by parallel workers
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
prepare_wgs_cell_line <- function(
  chrom_names, chrom_coord, tumourbam, tumourname,
  g1000lociprefix, g1000allelesprefix, gamma_ivd = 1e5,
  kmin_ivd = 50, centromere_noise_seg_size = 1e6,
  centromere_dist = 5e5, min_het_dist = 1e5, gamma_logr = 100,
  length_adjacent = 5e4, gccorrectprefix, repliccorrectprefix,
  min_base_qual, min_map_qual, allele_counts_dir, min_normal_depth,
  libs
) {
  # Standardise Chr notation (removes 'chr' string if present; essential for cell_line_baf_logR)
  # Skipping modification of external files. Assuming files are correct or handled in R reading.

  tumour_prefix <- file.path(allele_counts_dir, tumourname)

  # Check existence of at least one file
  first_file <- paste0(tumour_prefix, "_alleleFrequencies_chr", chrom_names[1], ".txt")
  if (!file.exists(first_file)) {
    log_failure("Expected allele counts file not found: {first_file}")
  }

  # Obtain BAF and LogR from the raw allele counts of the cell line
  cl_data <- cell_line_baf_logR(
    TUMOURNAME = tumour_prefix,
    g1000alleles_prefix = g1000allelesprefix,
    chrom_names = chrom_names
  )
  # Reconstruct normal-pair allele count files for the cell line

  run_parallel_or_serial(seq_along(chrom_names), function(i) {
    # If we are in parallel mode, ensure the packages are loaded on the worker
    if (FALSE) {
      # The least shit way to load dependencies inside a worker
      # This replaces the .packages argument from foreach
      requireNamespace("copynumber", quietly = TRUE)
      requireNamespace("ggplot2", quietly = TRUE)
      requireNamespace("grid", quietly = TRUE)
    }

    # Execute the reconstruction
    cell_line_reconstruct_normal(
      TUMOURNAME = tumourname,
      NORMALNAME = paste(tumourname, "_normal", sep = ""),
      chrom_coord = chrom_coord,
      chrom = i,
      CL_OHET = cl_data$OHET,
      CL_AL = cl_data$AL,
      CL_AC = cl_data$AC,
      CL_LogR = cl_data$LogR,
      GAMMA_IVD = gamma_ivd,
      KMIN_IVD = kmin_ivd,
      CENTROMERE_NOISE_SEG_SIZE = centromere_noise_seg_size,
      CENTROMERE_DIST = centromere_dist,
      MIN_HET_DIST = min_het_dist,
      GAMMA_LOGR = gamma_logr,
      LENGTH_ADJACENT = length_adjacent
    )
  }, libs)

  if (length(list.files(pattern = "normal_alleleFrequencies")) == length(chrom_names)) {
    log_info("STEP 2 - Normal allelecounts reconstruction - completed")
  } else {
    log_failure("Missing 'normal' allelecount files - all chromosomes NOT reconstructed")
  }

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
