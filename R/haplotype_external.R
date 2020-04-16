
#' Split a single vcf into separate vcfs for each chromosome
#' @param chrom_names Names of the chromosomes
#' @param externalHaplotypeFile Full path of the external vcf containing phased haplotypes
#' @param outprefix Full path and prefix of the output files
#' @author jdemeul
#' @export
split_input_haplotypes <- function(chrom_names, externalhaplotypefile=NA, outprefix) {

  if (is.na(externalhaplotypefile)) return(NULL)
  
  hetsnps <- VariantAnnotation::readVcf(file = externalhaplotypefile,
                     param = VariantAnnotation::ScanVcfParam(fixed = "ALT", info = NA, geno = c("GT", "PS"), trimEmpty = T))
  
  hetsnps <- split(x = hetsnps, f = GenomicRanges::seqnames(hetsnps))
  hetsnps <- hetsnps[chrom_names]
  
  lapply(X = 1:length(chrom_names), FUN = function(idx, chrom_names, snps, outbase) {
    VariantAnnotation::writeVcf(obj = snps[[chrom_names[idx]]], filename = paste0(outbase, idx, ".vcf"))
  }, snps = hetsnps, outbase = outprefix, chrom_names = chrom_names)
  
  return(NULL)
}



#' Combine imputation results with external haplotype blocks 
#' @param chrom_names Names of the chromosomes
#' @param chrom Index of the chromosome for which to reconstruct haplotypes
#' @param imputedHaplotypeFile Full path to the imputed haplotyope file for the indexed chromosome
#' @param externalHaplotypeFile Full path to the vcf containing haplotype blocks for the indexed chromosome (default NA)
#' @author jdemeul
#' @export
input_known_haplotypes = function(chrom_names, chrom, imputedHaplotypeFile, externalHaplotypeFile=NA) {

  if (is.na(externalHaplotypeFile)) return(NULL)
  
  # read BB phasing input
  bbphasin <- read_imputed_output(file = imputedHaplotypeFile)
  
  # turn into GRanges and subset for het SNPs
  bbphasingr <- GenomicRanges::GRanges(seqnames = rep(chrom_names[chrom], nrow(bbphasin)), ranges = IRanges::IRanges(start = bbphasin$pos, width = 1))
  S4Vectors::mcols(bbphasingr) <- bbphasin[, c("alt", "hap1", "hap2")]
  bbphasingr <- bbphasingr[which(xor(bbphasingr$hap1 == 1, bbphasingr$hap2 == 1))]
  
  # load vcf containing external haplotyped variants
  hetsnps <- suppressWarnings(VariantAnnotation::readVcf(file = externalHaplotypeFile,
                                      param = VariantAnnotation::ScanVcfParam(fixed = "ALT", info = NA, geno = c("GT", "PS"), trimEmpty = T)))
  
  # subset to phased het SNPs on chrom & drop any multiallelic var & indels if present
  hetsnps <- hetsnps[which(VariantAnnotation::geno(hetsnps)$GT %in% c("0|1", "1|0"))]
  hetsnps <- hetsnps[which(lengths(VariantAnnotation::alt(hetsnps)) == 1)]
  hetsnps <- hetsnps[which(S4Vectors::nchar(VariantAnnotation::ref(hetsnps)) == 1 & unlist(S4Vectors::nchar(VariantAnnotation::alt(hetsnps))) == 1)]
  
  # e.g. if no phasing on X, no need to continue
  if (length(hetsnps) == 0) return(NULL)
  
  # match Battenberg het SNPs with those in external file
  snvoverlaps <- GenomicRanges::findOverlaps(query = bbphasingr, subject = hetsnps, type = "equal")
  # and make sure we're phasing the same REF/ALT alleles (ref will always be the same)
  snvoverlaps_sub <- snvoverlaps[which(bbphasingr[S4Vectors::queryHits(snvoverlaps)]$alt ==
                                         as.character(unlist(VariantAnnotation::alt(hetsnps[S4Vectors::subjectHits(snvoverlaps)]))))]
  
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
  bbphasingr <- GenomicRanges::GRangesList(split(x = bbphasingr, f = bbphasingr$hap1_10X != ""), compress = F)
  nearestidxs <- GenomicRanges::nearest(x = bbphasingr$'FALSE', subject = bbphasingr$'TRUE', select = "arbitrary")
  bbphasingr$'FALSE'$isH1 <- bbphasingr$'TRUE'$isH1[nearestidxs]
  bbphasingr$'FALSE'$PS <- bbphasingr$'TRUE'$PS[nearestidxs]
  bbphasingr <- GenomicRanges::sort(unlist(bbphasingr, use.names = F))
  
  # build final haplotypes by flipping blocks according to imputation
  # last haplotype assignment of first block must match first haplotype assignment of second block
  psrle <- S4Vectors::Rle(bbphasingr$PS)
  flip <- cumsum(c(F, bbphasingr$isH1[S4Vectors::start(psrle)[-1]] == bbphasingr$isH1[S4Vectors::end(psrle)[-S4Vectors::nrun(psrle)]])) %% 2
  S4Vectors::runValue(psrle) <- flip
  bbphasingr$isH1 <- ifelse(as.vector(psrle, mode = "logical"), !bbphasingr$isH1, bbphasingr$isH1)
  bbphasingr$hapfinal <- ifelse(bbphasingr$isH1, bbphasingr$hap1, bbphasingr$hap2)
  
  # reinsert the phased het SNP haplotypes into the total chromosomal haplotypes
  matchidxs <- match(x = GenomicRanges::start(bbphasingr), table = bbphasin$pos)
  bbphasin[matchidxs, "hap1"] <- bbphasingr$hapfinal
  bbphasin[matchidxs, "hap2"] <- abs(bbphasin[matchidxs, "hap1"] - 1)
  
  # backup original imputedHaplotypeFile
  file.rename(from = imputedHaplotypeFile, to = gsub(pattern = "\\.txt$", replacement = ".noExt.txt", x = imputedHaplotypeFile))
  
  # and write new version
  write.table(x = bbphasin, file=imputedHaplotypeFile, row.names=F, col.names=F, quote=F, sep="\t")
  return(NULL)
}



#' Writes the imputation and copy number phased haplotypes to a vcf
#' @param tumourname Sample name
#' @param SNPfiles Character vector of the paths to the alleleFrequencies files, ordered by chromosome index
#' @param imputedHaplotypeFiles Character vector of the paths to the impute_output files, ordered by chromosome index
#' @param outprefix Prefix to write the output vcf files to
#' @param chrom_names Names of the chromosomes
#' @param include_homozygous Include homozygous SNPs in the output vcf file
#' @author jdemeul
#' @export
write_battenberg_phasing <- function(tumourname, SNPfiles, imputedHaplotypeFiles, bafsegmented_file, outprefix, chrom_names, include_homozygous = F) {
  
  bafsegmented <- read_bafsegmented(bafsegmented_file)[, c("Chromosome", "Position", "BAFphased", "BAFseg")]
  bafsegmented <- split(x = bafsegmented[, c("Position", "BAFphased", "BAFseg")], f = bafsegmented$Chromosome)
  
  for (chrom in 1:length(chrom_names)) {
    
    # read allele counts and imputed haplotypes (for the actually used alleles & loci)
    snp_data <- read_alleleFrequencies(SNPfiles[chrom])
    allele_data <- read_imputed_output(imputedHaplotypeFiles[chrom])[, c("pos", "ref", "alt", "hap1", "hap2")]
    merge_data <- merge(x = allele_data, y = snp_data, by.x = "pos", by.y = "POS", sort = F)
    
    # map counts to ref/alt
    merge_data$ref_count <- ifelse(merge_data$ref == "A", merge_data$Count_A,
                                   ifelse(merge_data$ref == "C", merge_data$Count_C,
                                          ifelse(merge_data$ref == "G", merge_data$Count_G, merge_data$Count_T)))
    merge_data$alt_count <- ifelse(merge_data$alt == "A", merge_data$Count_A,
                                   ifelse(merge_data$alt == "C", merge_data$Count_C,
                                          ifelse(merge_data$alt == "G", merge_data$Count_G, merge_data$Count_T)))
    merge_data <- cbind(merge_data[, c("CHR", "pos", "ref", "alt", "ref_count", "alt_count", "hap1", "hap2")], BAF = merge_data$alt_count/(merge_data$ref_count+merge_data$alt_count))
    
    # add in the segmented BAF values and start creating output vcf
    merge_data <- merge(x = merge_data, y = bafsegmented[[chrom_names[chrom]]], by.x = "pos", by.y = "Position",
                        all.x = include_homozygous, sort = T)
    
    bbphasing_vr <- VariantAnnotation::VRanges(seqnames = merge_data$CHR, ranges = IRanges::IRanges(start = merge_data$pos, width = 1),
                            ref = merge_data$ref, alt = merge_data$alt, 
                            totalDepth = merge_data$ref_count+merge_data$alt_count,
                            refDepth = merge_data$ref_count, altDepth = merge_data$alt_count)
    
    # assign the genotypes based on flipping of individual BAF values in regions of allelic imbalance according to BAFseg
    S4Vectors::mcols(bbphasing_vr)$GT <- ifelse(is.na(merge_data$BAFphased), paste0(merge_data$hap1, "|", merge_data$hap2), 
                                     ifelse(merge_data$BAFseg > 0.525 | is.na(merge_data$BAFseg), 
                                            ifelse(abs(merge_data$BAFphased-merge_data$BAF) < 1e-5, "1|0", "0|1"),
                                            ifelse(abs(merge_data$BAFphased-merge_data$BAF) < 1e-5, "1/0", "0/1")))
    
    # add phase set annotation based on segmented BAF: every segment = phase set
    S4Vectors::mcols(bbphasing_vr)$PS <- as.integer(NA)
    phasedidx <- which(merge_data$BAFseg > 0.525)
    if (length(phasedidx) > 0) {
      hetsegrle <- Rle(merge_data$BAFseg[phasedidx])
      S4Vectors::mcols(bbphasing_vr)$PS[phasedidx] <- rep(start(bbphasing_vr)[phasedidx][start(hetsegrle)], S4Vectors::runLength(hetsegrle))
      
      if (length(phasedidx) < nrow(merge_data)) {
        S4Vectors::mcols(bbphasing_vr)$PS[-phasedidx] <- S4Vectors::mcols(bbphasing_vr)$PS[phasedidx][GenomicRanges::nearest(x = bbphasing_vr[-phasedidx], subject = bbphasing_vr[phasedidx], select = "arbitrary")]
      }
    }
    
    # write out vcf
    VariantAnnotation::sampleNames(bbphasing_vr) <- tumourname
    
    VariantAnnotation::writeVcf(obj = bbphasing_vr, filename = paste0(outprefix, chrom, ".vcf"), index = F)
    
  }
  return(NULL)
}





# combine Battenberg phasing from multiple samples, requires output of the Battenberg phasing as VCF files
# combine_battenberg_phasing <- function(bbdirs, chrom_names) {


#' Generates haplotype blocks, MSAI results, and plots from phasing information contained in multisample Battenberg runs 
#' @param chrom Index of the chromosome for which to obtain haplotypes
#' @param bbphasingprefixes Vector containing prefixes of the Battenberg_phased_chr files for the multiple samples
#' @param maxlag Maximal number of upstream SNPs used to inform the haplotype at another SNPs
#' @param relative_weight_balanced Relative weight to give to haplotype info from a sample without allelic imbalance in the region (default 0.25)
#' @param outdir Folder to write the output to
#' @param plotting Should the multisample phasing plots be made? (default TRUE)
#' @author jdemeul
#' @export
get_multisample_phasing <- function(chrom, bbphasingprefixes, maxlag = 100, relative_weight_balanced = .25, outdir, plotting = T) {

  vcfs <- lapply(X = paste0(bbphasingprefixes, chrom, ".vcf"), FUN = VariantAnnotation::readVcf)
  samplenames <- sapply(X = vcfs, FUN = function(x) samples(header(x)))
  
  # get common hetSNP loci
  temp <- do.call(c, lapply(X = vcfs, FUN = rowRanges))
  commonloci <- unique(names(which(GenomicRanges::countOverlaps(query = temp, type = "equal", drop.self = F, drop.redundant = F) == length(vcfs))))
  vcfs_common <- lapply(X = vcfs, FUN = function(x, commonloci) GenomicRanges::sort(x[commonloci]), commonloci = commonloci)
  
  # clean up
  rm(vcfs, temp, commonloci)
  
  # go through each vcf and add relevant columns as appropriate
  loci <- GenomicRan ::rowRanges(vcfs_common[[1]])
  for (vcfidx in 1:length(vcfs_common)) {
    # add the genotype, BAF and phaseblock info for each sample to all common loci
    singlevcf <- vcfs_common[[vcfidx]]
    sid <- samples(header(singlevcf))
    adddf <- DataFrame(Major = geno(singlevcf)$GT[,1], #Major = as.integer(ifelse(test = grepl(pattern = "|", x = geno(singlevcf)$GT, fixed = T), substr(x = geno(singlevcf)$GT, 1, 1), NA)),
                       BAF = geno(singlevcf)$AD[,1,2]/rowSums(geno(singlevcf)$AD[,1,]),
                       PS = geno(singlevcf)$PS[,1])
    colnames(adddf) <- paste0(sid, "_", colnames(adddf))
    S4Vectors::mcols(loci) <- cbind(S4Vectors::mcols(loci), adddf)
  }
  
  
  # get call for alt-ref switches at different lag intervals 1:maxlag
  # also keep track of which are evidenced by allelic imbalance in ≥ 1 sample and downweight the inference contribution from the other samples to relative_weight_balanced
  gtswitcheslist <- list()
  evidencelist <- list()
  
  for (lag in 1:maxlag) {
    # lag <- 1
    gtswitcheslist[[lag]] <- rbind(matrix(NA, nrow = lag - 1, ncol = length(vcfs_common)), apply(MARGIN = 2, X = S4Vectors::mcols(loci)[,grep(pattern = "Major", x = colnames(S4Vectors::mcols(loci)))],
                                                                                                 FUN = function(x, lag) abs(diff(as.integer(substr(x,1,1)), lag = lag)), lag = lag))
    
    # check whether all are phased, note that the filter takes into account past values only here! So needs to be shifted in next step
    evidencelist[[lag]] <- apply(MARGIN = 2, X = S4Vectors::mcols(loci)[,grep(pattern = "Major", x = colnames(S4Vectors::mcols(loci)))],
                                 FUN = function(x, lag) filter(x = grepl(pattern = "|", x = x, fixed = T), filter = rep(1, lag + 1), sides = 1) == lag+1, lag = lag)
    # and they have the same PS
    # evidencelist[[lag]] <- (evidencelist[[lag]][-1,] * rbind(matrix(NA, nrow = lag-1, ncol = length(vcfs_common)), apply(MARGIN = 2, X = mcols(loci)[,grep(pattern = "PS", x = colnames(mcols(loci)))],
    #                                                 FUN = function(x, lag) diff(x = x, lag = lag) == 0, lag = lag))) == 1
    evidencelist[[lag]] <- evidencelist[[lag]][-1,] * rbind(matrix(NA, nrow = lag-1, ncol = length(vcfs_common)), apply(MARGIN = 2, X = S4Vectors::mcols(loci)[,grep(pattern = "PS", x = colnames(S4Vectors::mcols(loci)))],
                                                                                                                        FUN = function(x, lag) diff(x = x, lag = lag) == 0, lag = lag))
    evidencelist[[lag]][evidencelist[[lag]] == 0] <- relative_weight_balanced
    evidencelist[[lag]] <- evidencelist[[lag]]/rowSums(evidencelist[[lag]])
  }
  
  # initiate the vector which will cntain the combined phased haplotype
  haplovect <- as.integer(rep(NA, length(loci)))
  
  # start with a simple majorty call for the first hetSNP
  haplovect[1] <- as.integer(names(sort(table(substr(unlist(S4Vectors::mcols(loci)[1,grep(pattern = "Major", x = colnames(S4Vectors::mcols(loci))), drop = T]),1,1)), decreasing = T)[1]))
  
  # votes for next positions integrate more laged inferences
  for (pos in 2:length(loci)) {
    nvotesalt <- 0
    if (pos - 1 > maxlag) maxlag_used <- maxlag else maxlag_used <- pos - 1
    lagwsum <- sum(1:maxlag_used) # used to downweight larger distances
    for (lag in 1:maxlag_used) {
      nvotesalt <- nvotesalt + sum(abs(haplovect[pos-lag] - gtswitcheslist[[lag]][pos-1, ])*evidencelist[[lag]][pos-1,]) * (maxlag_used + 1 - lag) / lagwsum
    }
    haplovect[pos] <- round(nvotesalt)
    # haplovect[pos] <- round(nvotesalt/maxlag_used)
  }
  
  # write out the joint phasing
  jointphasing_vr <- VariantAnnotation::VRanges(seqnames = GenomicRanges::seqnames(loci), ranges = GenomicRanges::ranges(loci), ref = loci$REF, alt = unlist(loci$ALT))
  
  # assign the genotypes based on flipping of individual BAF values in regions of allelic imbalance according to BAFseg
  S4Vectors::mcols(jointphasing_vr)$GT <- paste0(haplovect, "|", ifelse(haplovect == 0, 1, 0))
  
  # add phase set annotation based on segmented BAF: every segment = phase set
  S4Vectors::mcols(jointphasing_vr)$PS <- start(loci)[1]
  
  # write out vcf
  VariantAnnotation::sampleNames(jointphasing_vr) <- "multisample"
  VariantAnnotation::writeVcf(obj = jointphasing_vr, filename = paste0(outdir, "/multisample_phasing_chr", chrom, ".vcf"), index = F)
  
  
  #### added code for MSAI detection
  segrle <- S4Vectors::Rle(rowSums(as.data.frame(S4Vectors::mcols(loci)[, grep(pattern = "_PS", x = colnames(S4Vectors::mcols(loci)), value = T)])))
  jointsegments <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(loci[start(segrle)]), ranges = IRanges::IRanges(start = start(loci[start(segrle)]), end = start(loci[end(segrle)])))
  jointsegments$nhetsnps <- GenomicRanges::countOverlaps(query = jointsegments, subject = loci)
  # haplostr <- as.character(haplovect)
  for (samplename in samplenames) {
    # samplename <- samplenames[1]
    # browser()
    samplegt <- substr(x = S4Vectors::mcols(loci)[, paste0(samplename, "_Major")], 1,1) == haplovect
    # samplegt <- ifelse(grepl(pattern = "|", x = mcols(loci)[, paste0(samplename, "_Major")], fixed = T), substr(x = mcols(loci)[, paste0(samplename, "_Major")], 1,1), substr(x = mcols(loci)[, paste0(samplename, "_Major")], 1,1)) == haplovect
    S4Vectors::mcols(jointsegments)[, paste0(samplename, "_refmatch")] <- c(by(data = samplegt, INDICES = segrle, FUN = sum, simplify = T)) / jointsegments$nhetsnps
  }
  jointsegments$MSAI <- apply(X = S4Vectors::mcols(jointsegments)[, grep(pattern = "_refmatch", x = colnames(S4Vectors::mcols(jointsegments)))], MARGIN = 1, FUN = function(x) max(x, na.rm = T) - min(x, na.rm = T) > .9)
  # jointsegments$MSAI[jointsegments$nhetsnps < 100] <- NA
  
  write.table(x = as.data.frame(jointsegments), file = paste0(outdir, "/multisample_MSAI_chr", chrom, ".txt"), row.names = F, sep = "\t", quote = F)
  msaidf <- as.data.frame(jointsegments[which(jointsegments$MSAI), ])
  print(msaidf)
  
  # visualise the haplotypes for the different samples
  for (sample in samplenames) {
    df1 <- data.frame(pos = GenomicRanges::start(loci), newhap = haplovect, baf = ifelse(haplovect == 1, S4Vectors::mcols(loci)[ ,paste0(sample, "_BAF")], 1-S4Vectors::mcols(loci)[ ,paste0(sample, "_BAF")]))
    
    p1 <- ggplot()
    if (nrow(msaidf) > 0) {
      p1 <- p1 + geom_rect(data = msaidf, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 1), alpha = .1, color = "gray", size = 0)
    }
    p1 <- p1 + geom_point(data = df1, mapping = aes(x = pos, y = 1-baf), alpha = .6, colour = "#67a9cf", shape = 46, show.legend = F)
    p1 <- p1 + geom_point(data = df1, mapping = aes(x = pos, y = baf), alpha = .6, colour = "#ef8a62", shape = 46, show.legend = F) + theme_minimal()
    p1 <- p1 + labs(x = "Position", y = "BAF", title = paste0(sample, ": multisample phasing chr", chrom))
    
    ggsave(filename = paste0(outdir, "/", sample, "_multisample_phasing_chr", chrom, ".png"), plot = p1, width = 20, height = 5)
  }
  return(NULL)
}




