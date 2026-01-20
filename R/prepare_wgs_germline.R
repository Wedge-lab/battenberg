#' Obtain BAF and LogR from the Germline allele counts
#'
#' Function to generate BAF and LogR files based on allele counts of the Germline.
#' It also generates the input data required by the following 'germline_reconstruct_normal' function.
#' @param GERMLINENAME The germline name used for Battenberg (i.e. the Germline BAM file name without the '.bam' extension).
#' @param g1000alleles_prefix Prefix to where 1000 Genomes allele files can be found.
#' @param chrom_names A vector with allowed chromosome names.
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export

germline_baf_logR <- function(GERMLINENAME, g1000alleles_prefix, chrom_names) {
  # read heterozygous SNPs per chromosome for alleleCounter files & 1000G allele files####
  AC <- list() # alleleCounts
  AL <- list() # 1000G alleles
  MaC <- list() # matched alleleCounts
  OHET <- list() # HET SNP data

  for (chr in chrom_names) {
    # read in alleleCounter output for each chromosome (FAST)
    ac_file <- paste0(GERMLINENAME, "_alleleFrequencies_chr", chr, ".txt")
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
  germline <- GERMLINENAME

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
    germline = MACC$baf
  )
  names(BAF)[3] <- germline
  # Standardization
  BAF$Chromosome[BAF$Chromosome %in% c("23", 23)] <- "X"
  data.table::setorder(BAF, Chromosome, Position)
  data.table::fwrite(BAF, paste0(germline, "_mutantBAF.tab"), sep = "\t")
  rm(BAF)

  # Prepare and save LogR
  LogR_out <- data.frame(
    Chromosome = MACC$chr,
    Position = MACC$position,
    germline = MACC$logr
  )
  names(LogR_out)[3] <- germline
  LogR_out$Chromosome[LogR_out$Chromosome %in% c("23", 23)] <- "X"
  data.table::setorder(LogR_out, Chromosome, Position)
  data.table::fwrite(LogR_out, paste0(germline, "_mutantLogR.tab"), sep = "\t")

  return(list(
    OHET = OHET,
    AL   = AL,
    AC   = AC,
    LogR = LogR_out
  ))
}

#' Reconstruct normal-pair allele count files for Germlines
#'
#' Function to generate normal-pair allele count files based on IVD-PCF and inter-hetSNP logR-based LOH detection (IVD: Inter-Variant Distance, het: heterozygote)
#' This method reconstructs the normal-pair counts by using the allele counts of the Germline as template
#' It fills the detected LOH regions with evenly-distributed hetSNPs with the density estimated based on each chromosome in each germline sample
#' It essentially informs Battenberg of the location of hetSNPs across the genome in the germline sample
#' @param GERMLINENAME The germline name used for Battenberg (i.e. the germline BAM file name without the '.bam' extension)
#' @param NORMALNAME The normal name used for naming the generated normal-pair allele counts files
#' @param chrom_coord Full path to the file with chromosome coordinates including start, end and left/right centromere positions
#' @param chrom Chromosome number for which normal-pair will be reconstructed
#' @param GL_OHET List of observed heterozygous SNPs across all chromosomes generated within the germline_baf_logR function
#' @param GL_AL List of alleles at SNPs across all chromosomes generated within the germline_baf_logR function
#' @param GL_AC List of allele counts at SNPs across all chromosomes generated within the germline_baf_logR function
#' @param GL_LogR Dataframe of genomewide LogR values for SNPs across all chromosomes generated within the germline_baf_logR function
#' @param GAMMA_IVD The PCF gamma value for segmentation of 1000G hetSNP IVD values (Default 1e5)
#' @param KMIN_IVD The min number of SNPs to support a segment in PCF of 1000G hetSNP IVD values (Default 50)
#' @param CENTROMERE_DIST The minimum distance from the centromere to ignore in analysis due to the noisy nature of data in the vicinity of centromeres (Default 5e5)
#' @param MIN_HET_DIST The minimum distance for detecting higher resolution inter-hetSNP regions with potential LOH while accounting for inherent homozygote stretches (Default 1e5)
#' @param GAMMA_LOGR The PCF gamma value for confirming LOH within each inter-hetSNP candidate segment (Default 100)
#' @param LENGTH_ADJACENT The length of adjacent regions either side of a candidate inter-hetSNP LOH region to be plotted (Default 5e4)
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export

germline_reconstruct_normal <- function(
  GERMLINENAME, NORMALNAME,
  chrom_coord, chrom,
  GL_OHET, GL_AL, GL_AC,
  GL_LogR, GAMMA_IVD, KMIN_IVD,
  CENTROMERE_NOISE_SEG_SIZE,
  CENTROMERE_DIST, MIN_HET_DIST,
  GAMMA_LOGR, LENGTH_ADJACENT
) {
  # IDENTIFY REGIONS OF LOH #
  colClasses <- c(chr = "numeric", start = "numeric", cen.left.base = "numeric", cen.right.base = "numeric", end = "numeric")
  # Use fast I/O
  chr_loc <- data.table::fread(chrom_coord, colClasses = colClasses, header = TRUE, stringsAsFactors = FALSE)
  data.table::setDF(chr_loc)
  chr_loc$length <- (chr_loc$cen.left.base - chr_loc$start) + (chr_loc$end - chr_loc$cen.right.base)

  # STEP 2.0: identify LOH by IVD-PCF
  LOH <- list()
  PCF_folder <- "PCF_plots"
  if (!dir.exists(PCF_folder)) {
    dir.create(PCF_folder)
  }
  i <- chrom
  log_info("chrom={i}")
  pcf_input <- data.frame(chr = i, position = GL_OHET[[i]]$Position, IVD = (GL_OHET[[i]]$Position_dist_percent))
  pcf_input <- pcf_input[which(pcf_input$position < chr_loc[i, "cen.left.base"] - CENTROMERE_DIST | pcf_input$position > chr_loc[i, "cen.right.base"] + CENTROMERE_DIST), ]
  # use only regions covered with gcCorrect LogR range
  pcf_input <- pcf_input[which(pcf_input$position >= chr_loc[i, "start"] & pcf_input$position <= chr_loc[i, "end"]), ]
  PCF <- copynumber::pcf(pcf_input, gamma = GAMMA_IVD, kmin = KMIN_IVD)
  grDevices::pdf(paste0(
    PCF_folder, "/", GERMLINENAME, "_chr", i, "_PCF_plot.pdf"
  ))
  copynumber::plotChrom(pcf_input, PCF)
  grDevices::dev.off()
  PCF$diff <- PCF$end.pos - PCF$start.pos

  # Decide if there is any LOH based on PCF and chr_snp_density
  chr_snp_density <- nrow(pcf_input) / (pcf_input$position[nrow(pcf_input)] - pcf_input$position[1]) # density of HET SNPs across the region covered by HET SNPs
  # CALCULATE min_normal_snp_density#
  # minimum normal density for SNPs (in bps) is 3 x 10^-4 with median of 7 x 10^-4
  ####
  min_normal_snp_density <- 0.0001
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
        # FOR EXCLUSION: segment is short IVD region (default<1Mb) and endpos is over the p-arm limit (ending point)
        if (loh_regions$end.pos[j] > chr_loc$cen.left.base[i] && loh_regions$diff[j] < CENTROMERE_NOISE_SEG_SIZE) {
          noise <- append(noise, j)
        }
      }
      if (loh_regions$arm[j] == "q") {
        # FOR EXCLUSION: segment is short IVD region (default<1Mb) and startpos is below the q-arm limit (starting point)
        if (loh_regions$start.pos[j] < chr_loc$cen.right.base[i] && loh_regions$diff[j] < CENTROMERE_NOISE_SEG_SIZE) {
          noise <- append(noise, j)
        }
        # qARM of Chr 1,9,16 have large heterochromatin region next to centromere + 100kb tolerance for start of heterochromatin region
        if (loh_regions$start.pos[j] < (chr_loc$cen.right.base[i] + 1e5) && loh_regions$diff[j] > CENTROMERE_NOISE_SEG_SIZE && !is.na(match(chrom, c(1, 9, 16)))) {
          noise <- append(noise, j)
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
  #
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
  #
  ##
  # STEP 2 - get higher resolution LOH regions
  log_info("chrom={i}")
  # Use list for efficient non_LOH construction
  ac <- GL_AC[[i]]
  al <- GL_AL[[i]]
  names(ac) <- c("chr", "position", "1", "2", "3", "4", "depth")

  chr_interval <- c(chr_loc[i, "start"], chr_loc[i, "end"])
  if (!is.null(nrow(LOH[[i]])) && nrow(LOH[[i]]) > 0) {
    non_LOH_list <- list()
    for (j in 1:(nrow(LOH[[i]]) + 1)) {
      if (j == 1 && chr_interval[1] >= LOH[[i]]$start.pos[j]) {
        # LOH starts at or before interval start
      } else if (j == 1) {
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
    # Filter out invalid segments
    non_LOH <- non_LOH[non_LOH$diff > 0, ]
  }

  non_LOH <- non_LOH[order(non_LOH$start), ]

  # STEP 2.1: identify LOH by inter-HET SNP regions
  ohet <- GL_OHET[[i]]
  nSNPs <- as.numeric(nrow(GL_LogR))
  logr <- GL_LogR[which(GL_LogR$Chromosome == i), ]
  colnames(logr)[3] <- "LogR"
  logr$Position <- as.numeric(logr$Position)
  # if regions of non_LOH exist after IVD-PCF, run window-based search
  if (!is.null(non_LOH)) {
    pLOH_regions <- data.frame()
    if (is.na(match(i, c(13, 14, 15, 21, 22)))) {
      log_info("START {i} p ARM")
      PARM <- non_LOH[which(non_LOH$end <= chr_loc[i, ]$cen.left.base), ]
      if (nrow(PARM) > 0) {
        parm <- PARM
      } else if (nrow(PARM) == 0 && sum(non_LOH$diff) != 0) {
        parm <- data.frame(start = chr_interval[1], end = chr_loc[i, ]$cen.left.base - CENTROMERE_DIST)
      } else {
        log_info("unknown issue")
      }

      if (parm[nrow(parm), 1] < (parm[nrow(parm), 2] - CENTROMERE_DIST)) {
        # exclude the last CENTROMERE_DIST segment next to the centromere (left side) - too noisy
        parm[nrow(parm), 2] <- parm[nrow(parm), 2] - CENTROMERE_DIST
      } else {
        parm <- parm[-nrow(parm), ]
      }
      parm$diff <- parm$end - parm$start

      # search per non_LOH segment
      for (seg in seq_len(nrow(parm))) {
        LoH_list <- list()
        # IVD-based breakpoints for small regions#
        seg_ivd <- ohet[which(ohet$Position_dist >= MIN_HET_DIST & ohet$Position >= parm$start[seg] & ohet$Position <= parm$end[seg]), ]
        if (nrow(seg_ivd) > 0) {
          # Pre-calculate indices for logr to avoid repeated subsetting
          logr_in_seg_idx <- which(logr$Position >= parm$start[seg] & logr$Position <= parm$end[seg])
          if (length(logr_in_seg_idx) > 0) {
            logr_seg <- logr[logr_in_seg_idx, ]
            # Using findInterval to quickly get boundaries for each window
            starts_idx <- findInterval(seg_ivd$Position, logr_seg$Position) + 1
            ends_idx <- findInterval(seg_ivd$Position + seg_ivd$Position_dist, logr_seg$Position)

            for (j in seq_len(nrow(seg_ivd))) {
              if (starts_idx[j] > ends_idx[j]) next

              COV <- logr_seg[starts_idx[j]:ends_idx[j], ]
              cov <- mean(COV$LogR)
              denSNP <- nrow(COV) / (nSNPs / sum(chr_loc$length) * seg_ivd$Position_dist[j])

              if (!is.na(cov) && denSNP > 0.5) {
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
                  # Count SNPs in the actual LOH region
                  loh$N <- sum(COV$Position >= loh$start & COV$Position <= loh$end)
                  if (loh$N >= 10) {
                    LoH_list[[length(LoH_list) + 1]] <- loh
                  }
                }
              }
            }
          }
        }

        if (length(LoH_list) > 0) {
          LoH <- collapse::rowbind(LoH_list)
          # Combine smaller regions into larger regions of LOH
          LoH_regions_list <- list()
          if (nrow(LoH) > 0) {
            start <- LoH$start[1]
            end <- LoH$end[1]
            if (nrow(LoH) > 1) {
              for (j in 2:nrow(LoH)) {
                if (LoH$start[j] == end) {
                  end <- LoH$end[j]
                } else {
                  LoH_regions_list[[length(LoH_regions_list) + 1]] <- data.frame(chrom = i, arm = "p", start.pos = start, end.pos = end)
                  start <- LoH$start[j]
                  end <- LoH$end[j]
                }
              }
            }
            LoH_regions_list[[length(LoH_regions_list) + 1]] <- data.frame(chrom = i, arm = "p", start.pos = start, end.pos = end)
          }
          pLOH_regions <- rbind(pLOH_regions, collapse::rowbind(LoH_regions_list))
        }
      }
      if (nrow(pLOH_regions) > 0) {
        # pARM BAF/LogR plot(s)
        grDevices::pdf(paste0(GERMLINENAME, "_chr", i, "_", MIN_HET_DIST / 1e3, "k_based_pLOH_events.pdf"))
        suppressWarnings(
          for (s in seq_len(nrow(pLOH_regions))) {
            sBAF <- ggplot2::ggplot(
              ohet, ggplot2::aes(rlang::.data$Position, rlang::.data$baf)
            ) +
              ggplot2::geom_jitter() +
              ggplot2::ylim(0, 1) +
              ggplot2::geom_vline(
                xintercept = c(pLOH_regions$start.pos[s], pLOH_regions$end.pos[s]),
                col = "red", linetype = "longdash"
              ) +
              ggplot2::xlim(
                pLOH_regions$start.pos[s] - LENGTH_ADJACENT,
                pLOH_regions$end.pos[s] + LENGTH_ADJACENT
              ) +
              ggplot2::ggtitle(paste("pARM LOH region", s)) +
              ggplot2::labs(y = "BAF")
            sLogR <- ggplot2::ggplot(
              logr,
              ggplot2::aes(rlang::.data$Position, rlang::.data$LogR)
            ) +
              ggplot2::geom_jitter() +
              ggplot2::ylim(-5.2, 1.2) +
              ggplot2::geom_vline(
                xintercept = c(pLOH_regions$start.pos[s], pLOH_regions$end.pos[s]),
                col = "red", linetype = "longdash"
              ) +
              ggplot2::xlim(
                pLOH_regions$start.pos[s] - LENGTH_ADJACENT,
                pLOH_regions$end.pos[s] + LENGTH_ADJACENT
              )
            grid::grid.newpage()
            grid::grid.draw(
              rbind(ggplot2::ggplotGrob(sBAF),
                ggplot2::ggplotGrob(sLogR),
                size = "last"
              )
            )
          }
        )
        grDevices::dev.off()
        #
        log_info("Candidate LOH regions plotted for pARM")
      }
    } else {
      log_info("chr {i} is acrocentric - no p arm analysis")
    }
    # Q ARM RUN:
    log_info("START {i} q ARM")
    qLOH_regions <- data.frame()
    QARM <- non_LOH[which(non_LOH$start >= chr_loc[i, ]$cen.right.base), ]
    if (nrow(QARM) > 0) {
      qarm <- QARM
    } else if (nrow(QARM) == 0 && sum(non_LOH$diff) != 0) {
      qarm <- data.frame(start = chr_loc[i, ]$cen.right.base, end = chr_interval[2])
    } else {
      log_info("unknown issue")
    }
    # to exclude the first CENTROMERE_DIST next to the centromere (right side) - noisy
    qarm[1, 1] <- qarm[1, 1] + CENTROMERE_DIST
    qarm$diff <- qarm$end - qarm$start
    #
    # search per non_LOH segment
    for (seg in seq_len(nrow(qarm))) {
      LoH_list <- list()
      # IVD-based breakpoints for small regions#
      seg_ivd <- ohet[which(ohet$Position_dist >= MIN_HET_DIST & ohet$Position >= qarm$start[seg] & ohet$Position <= qarm$end[seg]), ]
      if (nrow(seg_ivd) > 0) {
        # Pre-calculate indices for logr
        logr_in_seg_idx <- which(logr$Position >= qarm$start[seg] & logr$Position <= qarm$end[seg])
        if (length(logr_in_seg_idx) > 0) {
          logr_seg <- logr[logr_in_seg_idx, ]
          starts_idx <- findInterval(seg_ivd$Position, logr_seg$Position) + 1
          ends_idx <- findInterval(seg_ivd$Position + seg_ivd$Position_dist, logr_seg$Position)

          for (j in seq_len(nrow(seg_ivd))) {
            if (starts_idx[j] > ends_idx[j]) next

            COV <- logr_seg[starts_idx[j]:ends_idx[j], ]
            cov <- mean(COV$LogR)
            denSNP <- nrow(COV) / (nSNPs / sum(chr_loc$length) * seg_ivd$Position_dist[j])

            if (!is.na(cov) && denSNP > 0.5) {
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
                if (loh$N >= 10) {
                  LoH_list[[length(LoH_list) + 1]] <- loh
                }
              }
            }
          }
        }
      }

      if (length(LoH_list) > 0) {
        LoH <- collapse::rowbind(LoH_list)
        LoH_regions_list <- list()
        if (nrow(LoH) > 0) {
          start <- LoH$start[1]
          end <- LoH$end[1]
          if (nrow(LoH) > 1) {
            for (j in 2:nrow(LoH)) {
              if (LoH$start[j] == end) {
                end <- LoH$end[j]
              } else {
                LoH_regions_list[[length(LoH_regions_list) + 1]] <- data.frame(chrom = i, arm = "q", start.pos = start, end.pos = end)
                start <- LoH$start[j]
                end <- LoH$end[j]
              }
            }
          }
          LoH_regions_list[[length(LoH_regions_list) + 1]] <- data.frame(chrom = i, arm = "q", start.pos = start, end.pos = end)
        }
        qLOH_regions <- rbind(qLOH_regions, collapse::rowbind(LoH_regions_list))
      }
    }
    if (nrow(qLOH_regions) > 0) {
      # qARM BAF/LogR plot(s)
      grDevices::pdf(paste0(GERMLINENAME, "_chr", i, "_", MIN_HET_DIST / 1e3, "k_based_qLOH_events.pdf"))
      suppressWarnings(
        for (s in seq_len(nrow(qLOH_regions))) {
          sBAF <- ggplot2::ggplot(ohet, ggplot2::aes(rlang::.data$Position, rlang::.data$baf)) +
            ggplot2::geom_jitter() +
            ggplot2::ylim(0, 1) +
            ggplot2::geom_vline(
              xintercept = c(
                qLOH_regions$start.pos[s],
                qLOH_regions$end.pos[s]
              ),
              col = "red", linetype = "longdash"
            ) +
            ggplot2::xlim(
              qLOH_regions$start.pos[s] - LENGTH_ADJACENT,
              qLOH_regions$end.pos[s] + LENGTH_ADJACENT
            ) +
            ggplot2::ggtitle(paste("qARM LOH region", s))
          sLogR <- ggplot2::ggplot(
            logr, ggplot2::aes(rlang::.data$Position, rlang::.data$LogR)
          ) +
            ggplot2::geom_jitter() +
            ggplot2::ylim(-5.2, 1.2) +
            ggplot2::geom_vline(
              xintercept = c(
                qLOH_regions$start.pos[s],
                qLOH_regions$end.pos[s]
              ),
              col = "red", linetype = "longdash"
            ) +
            ggplot2::xlim(
              qLOH_regions$start.pos[s] - LENGTH_ADJACENT,
              qLOH_regions$end.pos[s] + LENGTH_ADJACENT
            )
          grid::grid.newpage()
          grid::grid.draw(rbind(
            ggplot2::ggplotGrob(sBAF),
            ggplot2::ggplotGrob(sLogR),
            size = "last"
          ))
        }
      )
      grDevices::dev.off()
      #
      log_info("Candidate LOH regions plotted for qARM")
    }
    # STEP 2.2: clean-up LOH[[i]] and merge LOH regions of both methods

    noLOH <- NULL
    if (!is.null(nrow(LOH[[i]]))) {
      for (j in seq_len(nrow(LOH[[i]]))) {
        LOH[[i]]$logR[j] <- mean(logr[which(logr$Position >= LOH[[i]]$start.pos[j] & logr$Position <= LOH[[i]]$end.pos[j]), ][, 3])
        LOH[[i]]$nSNP[j] <- nrow(logr[which(logr$Position >= LOH[[i]]$start.pos[j] & logr$Position <= LOH[[i]]$end.pos[j]), ])
        LOH[[i]]$denSNP[j] <- LOH[[i]]$nSNP[j] / ((LOH[[i]]$end.pos[j] - LOH[[i]]$start.pos[j]) * nSNPs / sum(chr_loc$length))
        if (LOH[[i]]$logR[j] > -0.8 || LOH[[i]]$denSNP[j] < 0.5) {
          noLOH <- append(noLOH, j)
          log_info("j: '{j}'")
        }
      }
      LOH[[i]] <- LOH[[i]][-noLOH, ]
      LOH[[i]] <- ifelse(nrow(LOH[[i]]) == 0, 0, LOH[[i]])
    }

    LOH_regions <- data.frame()
    if (nrow(pLOH_regions) > 0) {
      LOH_regions <- rbind(LOH_regions, pLOH_regions)
    } else {
      log_info("no window-based LOH regions identified in p arm of non_LOH of IVD-PCF")
    }
    if (nrow(qLOH_regions) > 0) {
      if (nrow(LOH_regions) > 0) {
        LOH_regions <- collapse::rowbind(LOH_regions, qLOH_regions)
      } else {
        LOH_regions <- qLOH_regions
      }
    }

    if (nrow(LOH_regions) > 0) {
      if (!is.null(LOH[[i]]) && !is.null(nrow(LOH[[i]])) && nrow(LOH[[i]]) > 0) {
        LOH[[i]] <- collapse::rowbind(LOH[[i]][, c("chrom", "arm", "start.pos", "end.pos")], LOH_regions)
        LOH[[i]] <- LOH[[i]][order(LOH[[i]]$start.pos), ]
      } else {
        LOH[[i]] <- LOH_regions
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

  # STEP 3: RECONSTRUCT alleleCounter files for the pseudo-NORMAL sample
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
        SNP_interval <- non_LOH_length / max(1, nrow(GL_OHET[[i]]))
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
    # Check for centromeric SNPs not covered by LOH/non-LOH
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
    log_info("No changes made to the alleleCounter file - no LOH in chr {i}")
  }
}

#' Prepare data for impute
#'
#' @param chrom The chromosome for which impute input should be generated.
#' @param germline_allele_counts_file Output from the allele counter on the matched germline for this chromosome.
#' @param normal_allele_counts_file Output from the allele counter on the matched normal for this chromosome.
#' @param output_file File where the impute input for this chromosome will be written.
#' @param imputeinfofile Info file with impute reference information.
#' @param is_male Boolean denoting whether this sample is male (TRUE), or female (FALSE).
#' @param problem_loci_file A file containing genomic locations that must be discarded (optional).
#' @param use_loci_file A file containing genomic locations that must be included (optional).
#' @param heterozygous_filter The cutoff where a SNP will be considered as heterozygous (default 0.01).
#' @author dw9, sd11, Naser Ansari-Pour (BDI, Oxford)
#' @export
generate_impute_input_wgs_germline <- function(
  chrom,
  germline_allele_counts_file,
  normal_allele_counts_file,
  output_file,
  imputeinfofile,
  is_male,
  problem_loci_file = NA,
  use_loci_file = NA,
  heterozygous_filter = 0.1
) {
  # Load impute reference info
  impute_info <- parse_imputeinfofile(imputeinfofile, is_male, chrom = chrom)
  chrom_name <- unique(impute_info$chrom)

  # Load and combine known SNPs from legend files
  known_SNPs <- data.table::rbindlist(
    lapply(impute_info$impute_legend, data.table::fread, sep = " "),
    use.names = TRUE
  )
  data.table::setkeyv(known_SNPs, "position")

  # Filter problem loci (anti-join using base-style logic or setkey)
  if (!is.na(problem_loci_file) && problem_loci_file != "NA") {
    problemSNPs <- data.table::fread(
      problem_loci_file,
      sep = "\t",
      select = c("Chr", "Pos")
    )
    # Subset using standard logical indexing to avoid NSE warnings
    problemSNPs <- problemSNPs[problemSNPs[["Chr"]] == chrom_name, ]

    data.table::setkeyv(problemSNPs, "Pos")
    known_SNPs <- known_SNPs[!problemSNPs, on = c(position = "Pos")]
  }

  # Filter to explicitly allowed loci
  if (!is.na(use_loci_file) && use_loci_file != "NA") {
    use_loci <- data.table::fread(use_loci_file, sep = "\t")
    # Using standard column access
    goodSNPs <- use_loci[use_loci[["chr"]] == chrom_name, "pos", with = FALSE][[1]]
    known_SNPs <- known_SNPs[known_SNPs[["position"]] %in% goodSNPs, ]
  }

  # Load allele counts
  cnt_names <- c("chr", "position", "ref_base", "A", "C", "G", "T")

  snp_data <- data.table::fread(
    germline_allele_counts_file,
    sep = "\t",
    header = FALSE,
    comment.char = "#"
  )
  data.table::setnames(snp_data, cnt_names)

  normal_snp_data <- data.table::fread(
    normal_allele_counts_file,
    sep = "\t",
    header = FALSE,
    comment.char = "#"
  )
  data.table::setnames(normal_snp_data, cnt_names)

  if (nrow(snp_data) == 0) {
    log_failure("Germline allele counts file is empty: {germline_allele_counts_file}")
  }
  if (ncol(snp_data) < 7) {
    log_failure("Germline allele counts file has fewer than 7 columns: {germline_allele_counts_file}")
  }

  if (nrow(normal_snp_data) == 0) {
    log_failure("Normal allele counts file is empty: {normal_allele_counts_file}")
  }

  data.table::setkeyv(snp_data, "position")
  data.table::setkeyv(normal_snp_data, "position")

  # Join reference SNPs to observed data
  found_data <- known_SNPs[snp_data, nomatch = NULL][normal_snp_data, nomatch = NULL]

  n <- nrow(found_data)
  if (n == 0L) {
    log_failure("No SNPs matched between reference and allele counts")
  }

  # Compute BAF
  # Accessing columns by character strings to avoid NSE
  ref_cols <- paste0("i.", found_data[["a0"]])
  alt_cols <- paste0("i.", found_data[["a1"]])
  rows <- seq_len(n)

  # Column indexing via match ensures no variable binding issues
  ref_counts <- as.numeric(found_data[cbind(rows, match(ref_cols, names(found_data)))])
  alt_counts <- as.numeric(found_data[cbind(rows, match(alt_cols, names(found_data)))])

  BAFs <- alt_counts / (alt_counts + ref_counts)
  BAFs[is.nan(BAFs)] <- 0

  # Generate genotypes
  minBaf <- min(heterozygous_filter, 1 - heterozygous_filter)
  maxBaf <- max(heterozygous_filter, 1 - heterozygous_filter)

  genotypes <- matrix(0L, nrow = n, ncol = 3)
  genotypes[BAFs <= minBaf, 1] <- 1
  genotypes[BAFs > minBaf & BAFs < maxBaf, 2] <- 1
  genotypes[BAFs >= maxBaf, 3] <- 1

  genotype_dt <- data.table::as.data.table(genotypes)
  data.table::setnames(genotype_dt, c("G1", "G2", "G3"))

  # Assemble output
  # Use set() to modify by reference using a character string for the column name
  data.table::set(found_data, j = "snp.names", value = paste0("snp", seq_len(n)))
  found_data <- data.table::as.data.table(cbind(found_data, genotype_dt))

  output_cols <- c("snp.names", "id", "position", "a0", "a1", "G1", "G2", "G3")

  # Use with = FALSE to select columns by character vector
  data.table::fwrite(
    found_data[, output_cols, with = FALSE],
    file = output_file,
    sep = " ",
    col.names = FALSE,
    quote = FALSE
  )

  # Write sample_g.txt for sex chromosomes
  if (is.na(as.numeric(chrom_name))) {
    sample_g_file <- file.path(dirname(output_file), "sample_g.txt")
    sample_g_data <- data.table::data.table(
      ID_1 = c("0", "INDIVI1"),
      ID_2 = c("0", "INDIVI1"),
      missing = c("0", "0"),
      sex = c("D", "2")
    )
    data.table::fwrite(sample_g_data, sample_g_file, sep = " ")
  }

  invisible(NULL)
}

# Helper function to replicate stats::cor(matrix, vector) using collapse speed
# This performs C++ based scaling and a cross-product
fast_cor_vec <- function(X, y) {
  # Handle missing values equivalent to use = "complete.obs"
  keep <- stats::complete.cases(X, y)

  # Standardize using collapse (extremely fast)
  X_std <- collapse::fscale(as.matrix(X[keep, ]))
  y_std <- collapse::fscale(y[keep])

  # Correlation = (X'y) / (n - 1)
  n_obs <- sum(keep)
  res <- (crossprod(X_std, y_std) / (n_obs - 1))[, 1]
  return(res)
}

#' Function to correct LogR for waivyness that correlates with GC content
#' @param germline_LogR_file String pointing to the germline LogR output
#' @param outfile String pointing to where the GC corrected LogR should be written
#' @param correlations_outfile File where correlations are to be saved
#' @param gc_content_file_prefix String pointing to where GC windows for this reference genome can be
#' found. These files should be split per chromosome and this prefix must contain the full path until
#' chr in its name. The .txt extension is automatically added.
#' @param replic_timing_file_prefix Like the gc_content_file_prefix, containing replication timing info (supply NULL if no replication timing correction is to be applied)
#' @param chrom_names A vector containing chromosome names to be considered
#' @param recalc_corr_afterwards Set to TRUE to recalculate correlations after correction
#' @author jonas demeulemeester, sd11, Naser Ansari-Pour (BDI, Oxford)
#' @export
gc_correct_wgs_germline <- function(germline_LogR_file, outfile, correlations_outfile,
                                    gc_content_file_prefix, replic_timing_file_prefix,
                                    chrom_names, recalc_corr_afterwards = FALSE) {
  if (is.null(gc_content_file_prefix)) {
    log_failure("GC content reference files must be supplied to WGS GC content correction")
  }

  # Fast reading of LogR
  Germline_LogR <- read_logr(germline_LogR_file)

  log_info("Processing GC content data")
  chrom_idx <- seq_along(chrom_names)

  # Efficiently reading and binding GC data
  gc_files <- paste0(gc_content_file_prefix, chrom_idx, ".txt.gz")
  GC_data <- data.table::rbindlist(lapply(gc_files, read_gccontent))
  colnames(GC_data) <- c(
    "chr", "Position", paste0(c(25, 50, 100, 200, 500), "bp"),
    paste0(c(1, 2, 5, 10, 20, 50, 100), "kb")
  )

  # Optional Replication timing data
  if (!is.null(replic_timing_file_prefix)) {
    log_info("Processing replication timing data")
    replic_files <- paste0(replic_timing_file_prefix, chrom_idx, ".txt.gz")
    replic_data <- data.table::rbindlist(lapply(replic_files, read_replication))
  }

  # Efficient Loci Synchronization
  key_logr <- paste0(Germline_LogR$Chromosome, "_", Germline_LogR$Position)
  key_gc <- paste0(GC_data$chr, "_", GC_data$Position)

  locimatches <- collapse::fmatch(key_logr, key_gc)

  valid_idx <- which(!is.na(locimatches))
  matched_gc_idx <- locimatches[valid_idx]

  Germline_LogR <- Germline_LogR[valid_idx, ]
  GC_data <- GC_data[matched_gc_idx, ]

  if (!is.null(replic_timing_file_prefix)) {
    replic_data <- replic_data[matched_gc_idx, ]
  }

  rm(key_logr, key_gc, locimatches, valid_idx, matched_gc_idx)

  # Fast Correlation calculation
  # Replaced stats::cor and non-existent fcor with helper
  corr <- abs(
    fast_cor_vec(GC_data[, 3:ncol(GC_data)], Germline_LogR[[3]])
  )

  if (!is.null(replic_timing_file_prefix)) {
    corr_rep <- abs(
      fast_cor_vec(replic_data[, 3:ncol(replic_data)], Germline_LogR[[3]])
    )
  }

  # Identify best window sizes
  index_1kb <- which(names(corr) == "1kb")
  maxGCcol_insert <- names(which.max(corr[1:index_1kb]))
  index_100kb <- which(names(corr) == "100kb")
  maxGCcol_amplic <- names(which.max(corr[(index_1kb + 2):index_100kb]))

  if (!is.null(replic_timing_file_prefix)) {
    maxreplic <- names(which.max(corr_rep))
    log_info("Replication timing correlation: {paste(names(corr_rep), format(corr_rep, digits = 2), collapse = '; ')}")
    log_info("Replication dataset: {maxreplic}")
  }

  log_info("GC correlation: {paste(names(corr), format(corr, digits = 2), collapse = '; ')}")
  log_info("Short window size: {maxGCcol_insert}")
  log_info("Long window size: {maxGCcol_amplic}")

  logr_vec <- Germline_LogR[[3]]

  # Create spline design matrices
  X_ins <- splines::ns(GC_data[[maxGCcol_insert]], df = 5, intercept = TRUE)
  X_amp <- splines::ns(GC_data[[maxGCcol_amplic]], df = 5, intercept = TRUE)

  if (!is.null(replic_timing_file_prefix)) {
    X_rep <- splines::ns(replic_data[[maxreplic]], df = 5, intercept = TRUE)
    X_design <- cbind(X_ins, X_amp, X_rep)

    before_corr_df <- data.frame(
      windowsize = c(names(corr), names(corr_rep)),
      correlation = c(as.numeric(corr), as.numeric(corr_rep))
    )
  } else {
    X_design <- cbind(X_ins, X_amp)
    before_corr_df <- data.frame(
      windowsize = names(corr),
      correlation = as.numeric(corr)
    )
  }

  # Fast Linear Model via collapse
  coeffs <- collapse::flm(logr_vec, X_design)

  # Calculate residuals (Corrected LogR)
  Germline_LogR[, 3] <- logr_vec - (X_design %*% coeffs)

  rm(X_ins, X_amp, X_design, coeffs)
  if (!is.null(replic_timing_file_prefix)) rm(X_rep)

  data.table::fwrite(before_corr_df,
    file = gsub(".txt", "_beforeCorrection.txt", correlations_outfile),
    sep = "\t", quote = FALSE
  )

  if (!recalc_corr_afterwards) {
    rm(GC_data)
    if (exists("replic_data")) rm(replic_data)
  }

  data.table::fwrite(
    Germline_LogR[!is.na(Germline_LogR[[3]]), ],
    file = outfile, sep = "\t"
  )

  # Optional Post-correction Analysis
  if (recalc_corr_afterwards) {
    # Re-using the helper for consistency and speed
    post_corr <- abs(
      fast_cor_vec(GC_data[, 3:ncol(GC_data)], Germline_LogR[[3]])
    )

    if (!is.null(replic_timing_file_prefix)) {
      post_corr_rep <- abs(
        fast_cor_vec(
          replic_data[, 3:ncol(replic_data)], Germline_LogR[[3]]
        )
      )

      log_info("Replication timing correlation post correction: {paste(names(post_corr_rep), format(post_corr_rep, digits = 2), collapse = '; ')}")

      after_corr_df <- data.frame(
        windowsize = c(
          names(post_corr),
          names(post_corr_rep)
        ),
        correlation = c(
          as.numeric(post_corr),
          as.numeric(post_corr_rep)
        )
      )
    } else {
      after_corr_df <- data.frame(
        windowsize = names(post_corr),
        correlation = as.numeric(post_corr)
      )
    }

    log_info("GC correlation post correction: {paste(names(post_corr), format(post_corr, digits = 2), collapse = '; ')}")
    data.table::fwrite(
      after_corr_df,
      file = gsub(
        ".txt", "_afterCorrection.txt",
        correlations_outfile
      ),
      sep = "\t",
      quote = FALSE
    )
  }
}


#' Prepare WGS data of germline for haplotype construction
#'
#' This function performs part of the Battenberg WGS pipeline: Counting alleles, generating BAF and logR,
#' reconstructing normal-pair allele counts for the germline and performing GC content correction.
#'
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param chrom_coord Full path to the file with chromosome coordinates including start, end and left/right centromere positions
#' @param germlinebam Full path to the germline BAM file
#' @param germlinename Identifier to be used for germline output files (i.e. the germline BAM file name without the '.bam' extension).
#' @param g1000lociprefix Prefix path to the 1000 Genomes loci reference files
#' @param g1000allelesprefix Prefix path to the 1000 Genomes SNP allele reference files
#' @param gamma_ivd The PCF gamma value for segmentation of 1000G hetSNP IVD values (Default 1e5).
#' @param kmin_ivd The min number of SNPs to support a segment in PCF of 1000G hetSNP IVD values (Default 50)
#' @param centromere_noise_seg_size The maximum size of PCF segment to be removed as noise when it overlaps with the centromere due to the noisy nature of data (Default 1e6)
#' @param centromere_dist The minimum distance from the centromere to ignore in analysis due to the noisy nature of data in the vicinity of centromeres (Default 5e5)
#' @param min_het_dist The minimum distance for detecting higher resolution inter-hetSNP regions with potential LOH while accounting for inherent homozygote stretches (Default 1e5)

#' @param allele_counts_dir Directory containing the allele counts files
#' @param min_normal_depth Minimum depth required in the normal for a SNP to be included
#' @param libs Path to the R libraries to be used by parallel workers
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
prepare_wgs_germline <- function(
  chrom_names, chrom_coord, germlinebam,
  germlinename, g1000lociprefix, g1000allelesprefix,
  gamma_ivd = 1e5, kmin_ivd = 50,
  centromere_noise_seg_size = 1e6,
  centromere_dist = 5e5, min_het_dist = 2e3,
  gamma_logr = 100, length_adjacent = 5e4,
  gccorrectprefix, repliccorrectprefix,
  min_base_qual, min_map_qual,
  allele_counts_dir, min_normal_depth,
  libs
) {
  germline_prefix <- file.path(allele_counts_dir, germlinename)

  # Check existence of at least one file
  first_file <- paste0(germline_prefix, "_alleleFrequencies_chr", chrom_names[1], ".txt")
  if (!file.exists(first_file)) {
    log_failure("Expected allele counts file not found: {first_file}")
    log_failure("Missing allele counts file: {first_file}")
  }

  # Standardise Chr notation (removes 'chr' string if present)
  # Skipping modification of external files. Assuming files are correct or handled in R reading.
  # standardise_chr_notation_germline(GERMLINENAME = germlinename)

  # Obtain BAF and LogR from the raw allele counts of the germline
  cl_data <- germline_baf_logR(
    GERMLINENAME = germline_prefix,
    g1000alleles_prefix = g1000allelesprefix,
    chrom_names = chrom_names
  )

  run_parallel_or_serial(
    iterator = seq_along(chrom_names),
    func = function(i) {
      germline_reconstruct_normal(
        GERMLINENAME = germlinename,
        NORMALNAME = paste(germlinename, "_normal", sep = ""),
        chrom_coord = chrom_coord,
        chrom = i,
        GL_OHET = cl_data$OHET,
        GL_AL = cl_data$AL,
        GL_AC = cl_data$AC,
        GL_LogR = cl_data$LogR,
        GAMMA_IVD = gamma_ivd,
        KMIN_IVD = kmin_ivd,
        CENTROMERE_NOISE_SEG_SIZE = centromere_noise_seg_size,
        CENTROMERE_DIST = centromere_dist,
        MIN_HET_DIST = min_het_dist,
        GAMMA_LOGR = gamma_logr,
        LENGTH_ADJACENT = length_adjacent
      )
    }, libs
  )

  if (length(list.files(pattern = "normal_alleleFrequencies")) == length(chrom_names)) {
    log_info("STEP 2 - Normal allelecounts reconstruction - completed")
  } else {
    log_failure("Missing 'normal' allelecount files - all chromosomes NOT reconstructed")
  }

  # Perform GC correction
  gc_correct_wgs_germline(
    germline_LogR_file = paste(germlinename, "_mutantLogR.tab", sep = ""),
    outfile = paste(germlinename, "_mutantLogR_gcCorrected.tab", sep = ""),
    correlations_outfile = paste(germlinename, "_GCwindowCorrelations.txt", sep = ""),
    gc_content_file_prefix = gccorrectprefix,
    replic_timing_file_prefix = repliccorrectprefix,
    chrom_names = chrom_names
  )
}
