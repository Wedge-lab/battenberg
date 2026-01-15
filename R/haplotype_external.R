#' Split a single vcf into separate vcfs for each chromosome
#' @param chrom_names Names of the chromosomes
#' @param externalHaplotypeFile Full path of the external vcf containing phased haplotypes (Default: NA)
#' @param outprefix Full path and prefix of the output files
#' @author jdemeul
#' @export
split_input_haplotypes <- function(chrom_names, externalhaplotypefile = NA, outprefix) {
  if (is.na(externalhaplotypefile)) {
    return(NULL)
  }

  hetsnps <- VariantAnnotation::readVcf(
    file = externalhaplotypefile,
    param = VariantAnnotation::ScanVcfParam(fixed = "ALT", info = NA, geno = c("GT", "PS"), trimEmpty = TRUE)
  )

  hetsnps <- split(x = hetsnps, f = GenomicRanges::seqnames(hetsnps))
  hetsnps <- hetsnps[chrom_names]

  lapply(X = chrom_names, FUN = function(chrom, chrom_names, snps, outbase) {
    VariantAnnotation::writeVcf(obj = snps[[chrom]], filename = paste0(outbase, chrom, ".vcf"))
  }, snps = hetsnps, outbase = outprefix, chrom_names = chrom_names)

  return(NULL)
}


#' Combine imputation results with external haplotype blocks
#' @param chrom_names Names of the chromosomes
#' @param chrom  chromosome for which to reconstruct haplotypes
#' @param imputedHaplotypeFile Full path to the imputed haplotyope file for the indexed chromosome
#' @param externalHaplotypeFile Full path to the vcf containing haplotype blocks for the indexed chromosome (Default: NA)
#' @param oldfilesuffix Suffix to be added to the original imputedHaplotypeFile (Default: _noExt.txt)
#' @author jdemeul
#' @export
input_known_haplotypes <- function(chrom_names, chrom, imputedHaplotypeFile, externalHaplotypeFile = NA, oldfilesuffix = "_noExt.txt") {
  if (is.na(externalHaplotypeFile)) {
    return(NULL)
  }

  # read BB phasing input
  bbphasin <- read_imputed_output(filename = imputedHaplotypeFile)

  # turn into GRanges and subset for het SNPs
  bbphasingr <- GenomicRanges::GRanges(seqnames = rep(chrom, nrow(bbphasin)), ranges = IRanges::IRanges(start = bbphasin$pos, width = 1))
  S4Vectors::mcols(bbphasingr) <- bbphasin[, c("alt", "hap1", "hap2")]
  bbphasingr <- bbphasingr[which(xor(bbphasingr$hap1 == 1, bbphasingr$hap2 == 1))]

  # load vcf containing external haplotyped variants
  hetsnps <- suppressWarnings(VariantAnnotation::readVcf(
    file = externalHaplotypeFile,
    param = VariantAnnotation::ScanVcfParam(fixed = "ALT", info = NA, geno = c("GT", "PS"), trimEmpty = TRUE)
  ))

  # subset to phased het SNPs on chrom & drop any multiallelic var & indels if present
  hetsnps <- hetsnps[which(VariantAnnotation::geno(hetsnps)$GT %in% c("0|1", "1|0"))]
  hetsnps <- hetsnps[which(lengths(VariantAnnotation::alt(hetsnps)) == 1)]
  hetsnps <- hetsnps[which(S4Vectors::nchar(VariantAnnotation::ref(hetsnps)) == 1 & unlist(S4Vectors::nchar(VariantAnnotation::alt(hetsnps))) == 1)]

  # e.g. if no phasing on X, no need to continue
  if (length(hetsnps) == 0) {
    return(NULL)
  }

  # match Battenberg het SNPs with those in external file, take only ranges to avoid chrom names mismatch
  snvoverlaps <- IRanges::findOverlaps(query = IRanges::ranges(bbphasingr), subject = IRanges::ranges(hetsnps), type = "equal")

  # add the corresponding phaseblocks (PS) and genotypes (GT)
  bbphasingr$PS <- vector(mode = "integer", length = length(bbphasingr))
  bbphasingr$GT <- vector(mode = "character", length = length(bbphasingr))
  bbphasingr[S4Vectors::queryHits(snvoverlaps)]$PS <- VariantAnnotation::geno(hetsnps[S4Vectors::subjectHits(snvoverlaps)])$PS
  bbphasingr[S4Vectors::queryHits(snvoverlaps)]$GT <- VariantAnnotation::geno(hetsnps[S4Vectors::subjectHits(snvoverlaps)])$GT

  # extract external haplotype 1 and match to imputed haplotypes
  bbphasingr$hap1_10X <- substr(bbphasingr$GT, start = 1, stop = 1)
  bbphasingr$isH1 <- ifelse(bbphasingr$hap1_10X == "", NA, bbphasingr$hap1_10X == bbphasingr$hap1)

  # complete and extend the known haplotype blocks
  # by transfering imputed haplotypes to nearest non-phased het SNPs
  # bbphasingr <- GenomicRanges::GRangesList(split(x = bbphasingr, f = bbphasingr$hap1_10X != ""), compress = FALSE)
  bbphasingr <- methods::as(object = split(x = bbphasingr, f = bbphasingr$hap1_10X != ""), Class = "GRangesList")
  if (length(bbphasingr$"FALSE") > 0) {
    nearestidxs <- GenomicRanges::nearest(x = bbphasingr$"FALSE", subject = bbphasingr$"TRUE", select = "arbitrary")
    bbphasingr$"FALSE"$isH1 <- bbphasingr$"TRUE"$isH1[nearestidxs]
    bbphasingr$"FALSE"$PS <- bbphasingr$"TRUE"$PS[nearestidxs]
  }
  bbphasingr <- GenomicRanges::sort(unlist(bbphasingr, use.names = FALSE))

  # build final haplotypes by flipping blocks according to imputation
  # last haplotype assignment of first block must match first haplotype assignment of second block
  psrle <- S4Vectors::Rle(bbphasingr$PS)
  flip <- cumsum(c(FALSE, bbphasingr$isH1[S4Vectors::start(psrle)[-1]] == bbphasingr$isH1[S4Vectors::end(psrle)[-S4Vectors::nrun(psrle)]])) %% 2
  S4Vectors::runValue(psrle) <- flip
  bbphasingr$isH1 <- ifelse(as.vector(psrle, mode = "logical"), !bbphasingr$isH1, bbphasingr$isH1)
  bbphasingr$hapfinal <- ifelse(bbphasingr$isH1, bbphasingr$hap1, bbphasingr$hap2)

  # reinsert the phased het SNP haplotypes into the total chromosomal haplotypes
  matchidxs <- match(x = GenomicRanges::start(bbphasingr), table = bbphasin$pos)
  bbphasin[matchidxs, "hap1"] <- bbphasingr$hapfinal
  bbphasin[matchidxs, "hap2"] <- abs(bbphasin[matchidxs, "hap1"] - 1)

  # backup original imputedHaplotypeFile
  if (file.exists(imputedHaplotypeFile)) {
    file.copy(from = imputedHaplotypeFile, to = gsub(pattern = "\\.txt$", replacement = oldfilesuffix, x = imputedHaplotypeFile), overwrite = TRUE)
  }

  # and write new version
  data.table::fwrite(x = bbphasin, file = imputedHaplotypeFile, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  return(NULL)
}

#' Writes the imputation and copy number phased haplotypes to a VCF
#' @param tumourname Sample name
#' @param SNPfiles Character vector of alleleFrequency files (per chromosome)
#' @param imputedHaplotypeFiles Character vector of impute2 haplotype files
#' @param bafsegmented_file Path to BAFSegmented file
#' @param outprefix Output VCF prefix
#' @param chrom_names Chromosome names
#' @param include_homozygous Include homozygous SNPs (default FALSE)
#' @export
write_battenberg_phasing <- function(
  tumourname,
  SNPfiles,
  imputedHaplotypeFiles,
  bafsegmented_file,
  outprefix,
  chrom_names,
  include_homozygous = FALSE
) {
  ## ---- Load & standardize BAF segments ----
  baf_dt <- read_bafsegmented(bafsegmented_file)
  data.table::setDT(baf_dt)
  data.table::setnames(
    baf_dt,
    old = c("Chromosome", "chrom", "chr"),
    new = c("CHR", "CHR", "CHR"),
    skip_absent = TRUE
  )
  baf_dt[, Position := as.integer(Position)]
  data.table::setkey(baf_dt, CHR, Position)

  ## ---- Impute2 schema ----
  impute_cols <- c("index", "rsid", "Position", "ref", "alt", "hap1", "hap2")

  for (idx in seq_along(chrom_names)) {
    chrom <- chrom_names[[idx]]

    ## ---- SNP / allele frequency ----
    snp_dt <- data.table::fread(SNPfiles[[idx]])
    data.table::setDT(snp_dt)
    data.table::setnames(
      snp_dt,
      old = c("Chromosome", "Chr", "POS"),
      new = c("CHR", "CHR", "Position"),
      skip_absent = TRUE
    )
    snp_dt[, Position := as.integer(Position)]
    snp_dt[, CHR := chrom]
    data.table::setkey(snp_dt, Position)

    ## ---- Imputed haplotypes ----
    hap_dt <- data.table::fread(
      imputedHaplotypeFiles[[idx]],
      header = FALSE,
      col.names = impute_cols
    )
    hap_dt <- hap_dt[, .(Position, ref, alt, hap1, hap2)]
    hap_dt[, Position := as.integer(Position)]
    data.table::setkey(hap_dt, Position)

    ## ---- Merge SNP + haplotypes ----
    dt <- snp_dt[hap_dt, nomatch = NULL]
    if (nrow(dt) == 0L) next

    ## ---- Allele counts ----
    dt[, ref_count := data.table::fcase(
      dt[["ref"]] == "A", dt[["Count_A"]],
      dt[["ref"]] == "C", dt[["Count_C"]],
      dt[["ref"]] == "G", dt[["Count_G"]],
      dt[["ref"]] == "T", dt[["Count_T"]],
      default = NA_real_
    )]

    dt[, alt_count := data.table::fcase(
      dt[["alt"]] == "A", dt[["Count_A"]],
      dt[["alt"]] == "C", dt[["Count_C"]],
      dt[["alt"]] == "G", dt[["Count_G"]],
      dt[["alt"]] == "T", dt[["Count_T"]],
      default = NA_real_
    )]

    dt[, BAF := dt[["alt_count"]] / (dt[["ref_count"]] + dt[["alt_count"]])]

    ## ---- Merge BAF segments ----
    baf_chr <- baf_dt[CHR == chrom, .(Position, BAFphased, BAFseg)]
    data.table::setkey(baf_chr, Position)

    if (include_homozygous) {
      dt <- baf_chr[dt]
    } else {
      dt <- dt[baf_chr, nomatch = NULL]
    }
    if (nrow(dt) == 0L) next

    ## ---- Build VRanges ----
    vr <- VariantAnnotation::VRanges(
      seqnames = dt[["CHR"]],
      ranges = IRanges::IRanges(start = dt[["Position"]], width = 1),
      ref = dt[["ref"]],
      alt = dt[["alt"]],
      totalDepth = dt[["ref_count"]] + dt[["alt_count"]],
      refDepth = dt[["ref_count"]],
      altDepth = dt[["alt_count"]]
    )

    ## ---- Genotype logic ----
    gt_vec <- data.table::fcase(
      is.na(dt[["BAFphased"]]),
      paste0(dt[["hap1"]], "|", dt[["hap2"]]),
      dt[["BAFseg"]] > 0.525 | is.na(dt[["BAFseg"]]),
      ifelse(abs(dt[["BAFphased"]] - dt[["BAF"]]) < 1e-5, "1|0", "0|1"),
      default =
        ifelse(abs(dt[["BAFphased"]] - dt[["BAF"]]) < 1e-5, "1/0", "0/1")
    )

    ## ---- Phase set (PS) ----
    n <- nrow(dt)
    ps <- rep(NA_integer_, n)
    phased_idx <- which(dt[["BAFseg"]] > 0.525)

    if (length(phased_idx) > 0) {
      rle_seg <- S4Vectors::Rle(dt[["BAFseg"]][phased_idx])
      ps[phased_idx] <- rep(
        dt[["Position"]][phased_idx][S4Vectors::start(rle_seg)],
        S4Vectors::runLength(rle_seg)
      )

      unphased <- setdiff(seq_len(n), phased_idx)
      if (length(unphased) > 0) {
        nearest <- GenomicRanges::nearest(
          vr[unphased],
          vr[phased_idx],
          select = "arbitrary"
        )
        ps[unphased] <- ps[phased_idx][nearest]
      }
    } else {
      ps[] <- dt[["Position"]][1]
    }

    ## ---- Attach metadata ----
    S4Vectors::mcols(vr)$GT <- gt_vec
    S4Vectors::mcols(vr)$PS <- ps
    VariantAnnotation::sampleNames(vr) <- tumourname

    ## ---- Write VCF ----
    VariantAnnotation::writeVcf(
      vr,
      filename = paste0(outprefix, chrom, ".vcf"),
      index = FALSE
    )
  }

  invisible(NULL)
}

#' @param chrom chromosome for which to obtain haplotypes
#' @param bbphasingprefixes Vector containing prefixes of the Battenberg_phased_chr files for the multiple samples
#' @param maxlag Maximal number of upstream SNPs used to inform the haplotype at another SNPs
#' @param relative_weight_balanced Relative weight to give to haplotype info from a sample without allelic imbalance in the region (default 0.25)
#' @param outprefix Prefix of the ouput multisample phasing files
#' @author jdemeul
#' @export
get_multisample_phasing <- function(chrom, bbphasingprefixes, maxlag = 90, relative_weight_balanced = .25, outprefix) {
  vcfs <- lapply(X = paste0(bbphasingprefixes, chrom, ".vcf"), FUN = VariantAnnotation::readVcf)

  # get common hetSNP loci
  temp <- do.call(c, lapply(X = vcfs, FUN = SummarizedExperiment::rowRanges))
  commonloci <- unique(names(which(GenomicRanges::countOverlaps(query = temp, type = "equal", drop.self = FALSE, drop.redundant = FALSE) == length(vcfs))))
  vcfs_common <- lapply(X = vcfs, FUN = function(x, commonloci) GenomicRanges::sort(x[commonloci]), commonloci = commonloci)

  # clean up
  rm(vcfs, temp, commonloci)

  # go through each vcf and add relevant columns as appropriate
  loci <- SummarizedExperiment::rowRanges(vcfs_common[[1]])
  for (vcfidx in seq_along(vcfs_common)) {
    # add the genotype, BAF and phaseblock info for each sample to all common loci
    singlevcf <- vcfs_common[[vcfidx]]
    sid <- VariantAnnotation::samples(VariantAnnotation::header(singlevcf))
    adddf <- S4Vectors::DataFrame(
      Major = VariantAnnotation::geno(singlevcf)$GT[, 1], # Major = as.integer(ifelse(test = grepl(pattern = "|", x = geno(singlevcf)$GT, fixed = TRUE), substr(x = geno(singlevcf)$GT, 1, 1), NA)),
      # BAF = VariantAnnotation::geno(singlevcf)$AD[,1,2]/BiocGenerics::rowSums(VariantAnnotation::geno(singlevcf)$AD[,1,]),
      BAF = VariantAnnotation::geno(singlevcf)$AD[, 1, 2] / rowSums(VariantAnnotation::geno(singlevcf)$AD[, 1, ]),
      PS = VariantAnnotation::geno(singlevcf)$PS[, 1]
    )
    colnames(adddf) <- paste0(sid, "_", colnames(adddf))
    S4Vectors::mcols(loci) <- cbind(S4Vectors::mcols(loci), adddf)
  }


  # get call for alt-ref switches at different lag intervals 1:maxlag
  # also keep track of which are evidenced by allelic imbalance in >= 1 sample and downweight the inference contribution from the other samples to relative_weight_balanced
  gtswitcheslist <- list()
  evidencelist <- list()

  for (lag in 1:maxlag) {
    # lag <- 1
    gtswitcheslist[[lag]] <- rbind(matrix(NA, nrow = lag - 1, ncol = length(vcfs_common)), apply(
      MARGIN = 2, X = S4Vectors::mcols(loci)[, grep(pattern = "Major", x = colnames(S4Vectors::mcols(loci)))],
      FUN = function(x, lag) abs(diff(as.integer(substr(x, 1, 1)), lag = lag)), lag = lag
    ))

    # check whether all are phased, note that the filter takes into account past values only here! So needs to be shifted in next step
    evidencelist[[lag]] <- apply(
      MARGIN = 2,
      X = S4Vectors::mcols(loci)[, grep(pattern = "Major", x = colnames(S4Vectors::mcols(loci)))],
      FUN = function(x, lag) {
        # First, find positions where the pattern "|" exists
        logical_vector <- grepl(pattern = "|", x = x, fixed = TRUE)
        numeric_vector <- as.numeric(logical_vector)
        result <- rep(FALSE, length(numeric_vector))

        if (length(numeric_vector) > lag) {
          # Then apply time series smoothing using stats::filter
          smoothed <- stats::filter(x = numeric_vector, filter = rep(1, lag + 1), sides = 1)
          smoothed[is.na(smoothed)] <- 0

          # Check where the smoothed values equal lag+1
          result[seq_along(smoothed)] <- (smoothed == lag + 1)
        }
        return(result)
      },
      lag = lag
    )
    # and they have the same PS
    # evidencelist[[lag]] <- (evidencelist[[lag]][-1,] * rbind(matrix(NA, nrow = lag-1, ncol = length(vcfs_common)), apply(MARGIN = 2, X = mcols(loci)[,grep(pattern = "PS", x = colnames(mcols(loci)))],
    #                                                 FUN = function(x, lag) diff(x = x, lag = lag) == 0, lag = lag))) == 1
    evidencelist[[lag]] <- evidencelist[[lag]][-1, ] * rbind(matrix(NA, nrow = lag - 1, ncol = length(vcfs_common)), apply(
      MARGIN = 2, X = S4Vectors::mcols(loci)[, grep(pattern = "PS", x = colnames(S4Vectors::mcols(loci)))],
      FUN = function(x, lag) diff(x = x, lag = lag) == 0, lag = lag
    ))
    evidencelist[[lag]][evidencelist[[lag]] == 0] <- relative_weight_balanced
    evidencelist[[lag]] <- evidencelist[[lag]] / rowSums(evidencelist[[lag]])
  }

  # initiate the vector which will cntain the combined phased haplotype
  haplovect <- as.integer(rep(NA, length(loci)))

  # start with a simple majorty call for the first hetSNP
  haplovect[1] <- as.integer(names(sort(table(substr(unlist(S4Vectors::mcols(loci)[1, grep(pattern = "Major", x = colnames(S4Vectors::mcols(loci))), drop = T]), 1, 1)), decreasing = TRUE)[1]))

  # votes for next positions integrate more laged inferences
  for (pos in 2:length(loci)) {
    nvotesalt <- 0
    if (pos - 1 > maxlag) maxlag_used <- maxlag else maxlag_used <- pos - 1
    lagwsum <- sum(1:maxlag_used) # used to downweight larger distances
    for (lag in 1:maxlag_used) {
      nvotesalt <- nvotesalt + sum(abs(haplovect[pos - lag] - gtswitcheslist[[lag]][pos - 1, ]) * evidencelist[[lag]][pos - 1, ]) * (maxlag_used + 1 - lag) / lagwsum
    }
    haplovect[pos] <- round(nvotesalt)
    # haplovect[pos] <- round(nvotesalt/maxlag_used)
  }

  # write out the joint phasing
  jointphasing_vr <- VariantAnnotation::VRanges(seqnames = GenomicRanges::seqnames(loci), ranges = GenomicRanges::ranges(loci), ref = loci$REF, alt = unlist(loci$ALT))

  # assign the genotypes based on flipping of individual BAF values in regions of allelic imbalance according to BAFseg
  S4Vectors::mcols(jointphasing_vr)$GT <- paste0(haplovect, "|", ifelse(haplovect == 0, 1, 0))

  # add phase set annotation based on segmented BAF: every segment = phase set
  S4Vectors::mcols(jointphasing_vr)$PS <- GenomicRanges::start(loci)[1]

  # write out vcf
  VariantAnnotation::sampleNames(jointphasing_vr) <- "multisample"
  VariantAnnotation::writeVcf(obj = jointphasing_vr, filename = paste0(outprefix, chrom, ".vcf"), index = FALSE)

  # write out loci + haplovect to do MSAI detection and plotting after final multisample CN calling
  S4Vectors::mcols(loci)$multisample_haplo <- haplovect
  saveRDS(object = loci, file = paste0(outprefix, chrom, "_loci.RDS"))

  return(NULL)
}


#' Generates haplotype blocks, MSAI results, and plots from phasing information contained in multisample Battenberg runs
#' @param rdsprefix Prefix of the RDS files containing the multisample haplotypes and BAF
#' @param subclonesfiles Vectors containing the paths to the different subclones.txt files
#' @param chrom_names Names of the chromosomes
#' @param tumournames Vector of sample names
#' @param plotting Should the multisample phasing plots be made? (Default: TRUE)
#' @author jdemeul
#' @export
call_multisample_MSAI <- function(
  rdsprefix,
  subclonesfiles,
  chrom_names,
  tumournames,
  plotting = TRUE
) {
  # compile all CN results
  subclonescat <- lapply(
    X = subclonesfiles, FUN = function(x) utils::read.delim(file = x, as.is = TRUE)
  )
  imbalancedregions <- do.call(rbind, subclonescat)
  # add sample identifiers
  imbalancedregions$sampleid <- rep(
    x = tumournames, sapply(X = subclonescat, FUN = nrow)
  )
  # subset to regions which are imbalanced in at least 2 samples
  imbalancedregions <- imbalancedregions[which(imbalancedregions$nMaj1_A != imbalancedregions$nMin1_A | imbalancedregions$nMaj2_A != imbalancedregions$nMin2_A), ]
  imbalancedregions <- GenomicRanges::GRanges(
    seqnames = imbalancedregions$chr,
    ranges = IRanges::IRanges(
      start = imbalancedregions$startpos,
      end = imbalancedregions$endpos
    ),
    sampleid = imbalancedregions$sampleid
  )
  imbalancedregions_disj <- GenomicRanges::disjoin(imbalancedregions)
  imbalancedregions_disj <- imbalancedregions_disj[GenomicRanges::countOverlaps(query = imbalancedregions_disj, subject = imbalancedregions) > 1]

  # if nothing remains, stop here
  if (length(imbalancedregions_disj) == 0) {
    log_info("No recurrently copy number imbalanced regions")
    return(NULL)
  }

  # add the identifiers of aberrated samples to each region
  samplehits <- GenomicRanges::findOverlaps(query = imbalancedregions_disj, subject = imbalancedregions)
  S4Vectors::mcols(imbalancedregions_disj)$sampleids <- split(x = imbalancedregions$sampleid[S4Vectors::subjectHits(samplehits)], f = S4Vectors::queryHits(samplehits))

  # split per chromosome, keeping only the imbalanced ones
  imbalancedregions_disj <- methods::as(object = split(x = imbalancedregions_disj, f = GenomicRanges::seqnames(imbalancedregions_disj)), Class = "GRangesList")

  # for every chromosome with imbalance
  for (i in seq_along(chrom_names)) {
    chrom <- as.character(chrom_names[i])
    # load loci.RDS file and simplify genotype formatting
    loci <- readRDS(file = paste0(rdsprefix, chrom, "_loci.RDS"))
    S4Vectors::mcols(loci)[, paste0(tumournames, "_Major")] <- S4Vectors::DataFrame(apply(
      X = S4Vectors::mcols(loci)[, paste0(tumournames, "_Major")],
      MARGIN = 2, FUN = function(x) as.numeric(substr(x = x, start = 1, stop = 1))
    ))

    # if (length(imbalancedregions_disj[[chrom]]) > 0) {
    if (chrom %in% names(imbalancedregions_disj)) {
      # split loci by abberrated region, compare only ranges to avoid chr naming scheme mismatch
      locioverlaps <- IRanges::findOverlaps(query = IRanges::ranges(imbalancedregions_disj[[chrom]]), subject = IRanges::ranges(loci))
      imballoci <- split(x = loci[S4Vectors::subjectHits(locioverlaps)], f = S4Vectors::queryHits(locioverlaps), drop = FALSE)

      # now check for each region the GT of major allele (in imbalanced samples)
      imbalancedregions_disj[[chrom]] <- imbalancedregions_disj[[chrom]][unique(S4Vectors::queryHits(locioverlaps))]

      frac_consensus <- mapply(haps = imballoci, samples = imbalancedregions_disj[[chrom]]$sampleids, FUN = function(haps, samples) {
        colSums(x = S4Vectors::as.matrix(S4Vectors::mcols(haps)[, paste0(samples, "_Major")]) == S4Vectors::mcols(haps)[, "multisample_haplo"], na.rm = TRUE) / length(haps)
      }, SIMPLIFY = FALSE)

      # simplify notation and call MSAI
      imbalancedregions_disj[[chrom]]$frac_consensus <- sapply(X = frac_consensus, FUN = function(x) paste0(names(x), "=", round(x, digits = 2), collapse = ";"))
      imbalancedregions_disj[[chrom]]$msai <- sapply(X = frac_consensus, FUN = function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE) > .9)

      if (length(GenomicRanges::mcols(imbalancedregions_disj[[chrom]])$msai) > 0) {
        msaidf <- GenomicRanges::as.data.frame(imbalancedregions_disj[[chrom]][GenomicRanges::mcols(imbalancedregions_disj[[chrom]])$msai])
      } else {
        msaidf <- data.frame()
      }
    } else {
      msaidf <- data.frame()
    }

    if (plotting) {
      # Plot the resulting data
      df1 <- data.frame(
        pos = GenomicRanges::start(loci),
        haplo = S4Vectors::mcols(loci)$multisample_haplo,
        BAF = as.numeric(rep(NA, length(loci)))
      )

      # visualise the haplotypes for the different samples
      for (tumour in tumournames) {
        df1$BAF <- ifelse(df1$haplo == 1, S4Vectors::mcols(
          loci
        )[, paste0(tumour, "_BAF")],
        1 - S4Vectors::mcols(loci)[, paste0(tumour, "_BAF")]
        )

        p1 <- ggplot2::ggplot()
        if (nrow(msaidf) > 0) {
          p1 <- p1 + ggplot2::geom_rect(
            data = msaidf, mapping = ggplot2::aes(
              xmin = rlang::.data$start,
              xmax = rlang::.data$end,
              ymin = 0, ymax = 1
            ),
            alpha = .05, color = "gray", size = 0
          )
        }
        p1 <- p1 + ggplot2::geom_point(data = df1, mapping = ggplot2::aes(
          x = rlang::.data$pos, y = 1 - rlang::.data$BAF
        ), alpha = .6, colour = "#67a9cf", shape = 46, show.legend = FALSE)
        p1 <- p1 + ggplot2::geom_point(
          data = df1, mapping = ggplot2::aes(x = rlang::.data$pos, y = rlang::.data$BAF),
          alpha = .6, colour = "#ef8a62", shape = 46, show.legend = FALSE
        ) + ggplot2::theme_minimal()
        p1 <- p1 + ggplot2::labs(
          x = "Position", y = "BAF",
          title = paste0(tumour, ": multisample phasing chr", chrom)
        )

        ggplot2::ggsave(
          filename = paste0(tumour, "_multisample_phasing_chr", chrom, ".png"),
          plot = p1, width = 20, height = 5
        )
      }
    }
  }

  # write out final MSAI dataframe
  msaiout <- GenomicRanges::as.data.frame(unlist(imbalancedregions_disj, use.names = FALSE))
  list_cols <- sapply(msaiout, is.list)
  for (col in names(msaiout)[list_cols]) {
    msaiout[[col]] <- sapply(msaiout[[col]], function(x) paste(x, collapse = ","))
  }
  data.table::fwrite(x = msaiout[, -c(4:6)], file = paste0("multisample_MSAI.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  return(NULL)
}
