
#' Obtain allele counts for 1000 Genomes loci through external program alleleCount
#'
#' @param bam.file A BAM alignment file on which the counter should be run.
#' @param output.file The file where output should go.
#' @param g1000.loci A file with 1000 Genomes SNP loci.
#' @param min.base.qual The minimum base quality required for it to be counted (optional, default=20).
#' @param min.map.qual The minimum mapping quality required for it to be counted (optional, default=35).
#' @param allelecounter.exe A pointer to where the alleleCounter executable can be found (optional, default points to $PATH).
#' @author sd11
#' @export
getAlleleCounts = function(bam.file, output.file, g1000.loci, min.base.qual=20, min.map.qual=35, allelecounter.exe="alleleCounter") {
  cmd = paste(allelecounter.exe,
              "-b", bam.file,
              "-l", g1000.loci,
              "-o", output.file,
              "-m", min.base.qual,
              "-q", min.map.qual)


  # alleleCount >= v4.0.0 is sped up considerably on 1000G loci when run in dense-snp mode
  counter_version = system(paste(allelecounter.exe, "--version"), intern = T)
  if (as.integer(substr(x = counter_version, start = 1, stop = 1)) >= 4)
    cmd = paste(cmd, "--dense-snps")

  system(cmd, wait=T)
}


#' Obtain BAF and LogR from the allele counts
#'
#' @param tumourAlleleCountsFile.prefix Prefix of the allele counts files for the tumour.
#' @param normalAlleleCountsFile.prefix Prefix of the allele counts files for the normal.
#' @param figuresFile.prefix Prefix for output figures file names.
#' @param BAFnormalFile File where BAF from the normal will be written.
#' @param BAFmutantFile File where BAF from the tumour will be written.
#' @param logRnormalFile File where LogR from the normal will be written.
#' @param logRmutantFile File where LogR from the tumour will be written.
#' @param combinedAlleleCountsFile File where combined allele counts for tumour and normal will be written.
#' @param chr_names A vector with allowed chromosome names.
#' @param g1000file.prefix Prefix to where 1000 Genomes reference files can be found.
#' @param minCounts Integer, minimum depth required for a SNP to be included (optional, default=NA).
#' @param samplename String, name of the sample (optional, default=sample1).
#' @param seed A seed to be set for when randomising the alleles.
#' @author dw9, sd11
#' @export
getBAFsAndLogRs = function(tumourAlleleCountsFile.prefix, normalAlleleCountsFile.prefix, figuresFile.prefix, BAFnormalFile, BAFmutantFile, logRnormalFile, logRmutantFile, combinedAlleleCountsFile, chr_names, g1000file.prefix, minCounts=NA, samplename="sample1", seed=as.integer(Sys.time())) {

  set.seed(seed)

  input_data = concatenateAlleleCountFiles(tumourAlleleCountsFile.prefix, ".txt", chr_names)
  normal_input_data = concatenateAlleleCountFiles(normalAlleleCountsFile.prefix, ".txt", chr_names)
  allele_data = concatenateG1000SnpFiles(g1000file.prefix, ".txt", chr_names)
  
  # We're no longer stripping out the "chr", which is causing problems
  #allele_data[,1] = gsub("chr","",allele_data[,1])
  #normal_input_data[,1] = gsub("chr","",normal_input_data[,1])
  #input_data[,1] = gsub("chr","",input_data[,1])

  # Synchronise all the data frames
  chrpos_allele = paste(allele_data[,1], "_", allele_data[,2], sep="")
  chrpos_normal = paste(normal_input_data[,1], "_", normal_input_data[,2], sep="")
  chrpos_tumour = paste(input_data[,1], "_", input_data[,2], sep="")
  matched_data = Reduce(intersect, list(chrpos_allele, chrpos_normal, chrpos_tumour))

  allele_data = allele_data[chrpos_allele %in% matched_data,]
  normal_input_data = normal_input_data[chrpos_tumour %in% matched_data,]
  input_data = input_data[chrpos_tumour %in% matched_data,]

  # Clean up and reduce amount of unneeded data
  names(input_data)[1] = "CHR"
  names(normal_input_data)[1] = "CHR"

  normal_data = normal_input_data[,3:6]
  mutant_data = input_data[,3:6]

  # Obtain depth for both alleles for tumour and normal
  len = nrow(normal_data)
  normCount1 = normal_data[cbind(1:len,allele_data[,3])]
  normCount2 = normal_data[cbind(1:len,allele_data[,4])]
  totalNormal = normCount1 + normCount2
  mutCount1 = mutant_data[cbind(1:len,allele_data[,3])]
  mutCount2 = mutant_data[cbind(1:len,allele_data[,4])]
  totalMutant = mutCount1 + mutCount2

  # Clean up a few unused variables to save some memory
  rm(normal_data, mutant_data, allele_data, normal_input_data)

  # Clear SNPs where there is not enough coverage
  indices = 1:nrow(input_data)
  if(!is.na(minCounts)){
    print(paste("minCount=", minCounts,sep=""))
    # Only normal has to have min coverage, mutant must have at least 1 read to prevent division by zero
    indices = which(totalNormal>=minCounts & totalMutant>=1)

    totalNormal = totalNormal[indices]
    totalMutant = totalMutant[indices]
    normCount1 = normCount1[indices]
    normCount2 = normCount2[indices]
    mutCount1 = mutCount1[indices]
    mutCount2 = mutCount2[indices]
  }
  n = length(indices)

  normalBAF = vector(length=n, mode="numeric")
  mutantBAF = vector(length=n, mode="numeric")
  normalLogR = vector(length=n, mode="numeric")
  mutantLogR = vector(length=n, mode="numeric")

  # randomise A and B alleles
  selector = round(runif(n))
  normalBAF[which(selector==0)] = normCount1[which(selector==0)] / totalNormal[which(selector==0)]
  normalBAF[which(selector==1)] = normCount2[which(selector==1)] / totalNormal[which(selector==1)]
  mutantBAF[which(selector==0)] = mutCount1[which(selector==0)] / totalMutant[which(selector==0)]
  mutantBAF[which(selector==1)] = mutCount2[which(selector==1)] / totalMutant[which(selector==1)]

  normalLogR = vector(length=n, mode="integer") #assume that normallogR is 0, and normalise mutantLogR to normalLogR
  mutantLogR = totalMutant/totalNormal
  rm(selector)

  # Create the output data.frames
  germline.BAF = data.frame(Chromosome=input_data$CHR[indices], Position=input_data$POS[indices], baf=normalBAF)
  germline.LogR = data.frame(Chromosome=input_data$CHR[indices], Position=input_data$POS[indices], samplename=normalLogR)
  tumor.BAF = data.frame(Chromosome=input_data$CHR[indices], Position=input_data$POS[indices], baf=mutantBAF)
  tumor.LogR = data.frame(Chromosome=input_data$CHR[indices], Position=input_data$POS[indices], samplename=log2(mutantLogR/mean(mutantLogR, na.rm=T)))
  alleleCounts = data.frame(Chromosome=input_data$CHR[indices], Position=input_data$POS[indices], mutCountT1=mutCount1, mutCountT2=mutCount2, mutCountN1=normCount1, mutCountN2=normCount2)

  # Save data.frames to disk
  write.table(germline.BAF,file=BAFnormalFile, row.names=F, quote=F, sep="\t", col.names=c("Chromosome","Position",samplename))
  write.table(tumor.BAF,file=BAFmutantFile, row.names=F, quote=F, sep="\t", col.names=c("Chromosome","Position",samplename))
  write.table(germline.LogR,file=logRnormalFile, row.names=F, quote=F, sep="\t", col.names=c("Chromosome","Position",samplename))
  write.table(tumor.LogR,file=logRmutantFile, row.names=F, quote=F, sep="\t", col.names=c("Chromosome","Position",samplename))
  write.table(alleleCounts, file=combinedAlleleCountsFile, row.names=F, quote=F, sep="\t")

  # Plot the raw data using ASCAT
  # Manually create an ASCAT object, which saves reading in the above files again
  SNPpos = germline.BAF[,c("Chromosome", "Position")]
  ch = list()
  for (i in 1:length(chr_names)) {
    temp = which(SNPpos$Chromosome==chr_names[i])
    if (length(temp) == 0) {
	ch[[i]] = 0
    } else {
    	ch[[i]] = temp[1]:temp[length(temp)]
    }
  }

  ascat.bc = list(Tumor_LogR=as.data.frame(tumor.LogR[,3]), Tumor_BAF=as.data.frame(tumor.BAF[,3]),
                  Germline_LogR=as.data.frame(germline.LogR[,3]), Germline_BAF=as.data.frame(germline.BAF[,3]),
                  Tumor_LogR_segmented=NULL, Tumor_BAF_segmented=NULL, Tumor_counts=NULL, Germline_counts=NULL,
                  SNPpos=tumor.LogR[,1:2], chrs=chr_names, samples=c(samplename), chrom=split_genome(tumor.LogR[,1:2]),
                  ch=ch)

  ASCAT::ascat.plotRawData(ascat.bc) #, parentDir=figuresFile.prefix)
}

#' Prepare data for impute
#'
#' @param chrom The chromosome for which impute input should be generated.
#' @param tumour.allele.counts.file Output from the allele counter on the matched tumour for this chromosome.
#' @param normal.allele.counts.file Output from the allele counter on the matched normal for this chromosome.
#' @param output.file File where the impute input for this chromosome will be written.
#' @param imputeinfofile Info file with impute reference information.
#' @param is.male Boolean denoting whether this sample is male (TRUE), or female (FALSE).
#' @param problemLociFile A file containing genomic locations that must be discarded (optional).
#' @param useLociFile A file containing genomic locations that must be included (optional).
#' @param heterozygousFilter The cutoff where a SNP will be considered as heterozygous (default 0.01).
#' @author dw9, sd11
#' @export
generate.impute.input.wgs = function(chrom, tumour.allele.counts.file, normal.allele.counts.file, output.file, imputeinfofile, is.male, problemLociFile=NA, useLociFile=NA, heterozygousFilter=0.1) {

  # Read in the 1000 genomes reference file paths for the specified chrom
  impute.info = parse.imputeinfofile(imputeinfofile, is.male, chrom=chrom)
  chr_names = unique(impute.info$chrom)
  chrom_name = chrom
  
  #print(paste("GenerateImputeInput is.male? ", is.male,sep=""))
  #print(paste("GenerateImputeInput #impute files? ", nrow(impute.info),sep=""))

  # Read in the known SNP locations from the 1000 genomes reference files
  known_SNPs = read.table(impute.info$impute_legend[1], sep=" ", header=T, stringsAsFactors=F)
  if(nrow(impute.info)>1){
    for(r in 2:nrow(impute.info)){
      known_SNPs = rbind(known_SNPs, read.table(impute.info$impute_legend[r], sep=" ", header=T, stringsAsFactors=F))
    }
  }

  # filter out bad SNPs (streaks in BAF)
  if((problemLociFile != "NA") & (!is.na(problemLociFile))) {
    problemSNPs = read.table(problemLociFile, header=T, sep="\t", stringsAsFactors=F)
    problemSNPs = problemSNPs$Pos[problemSNPs$Chr==chrom_name]
    badIndices = match(known_SNPs$position, problemSNPs)
    known_SNPs = known_SNPs[is.na(badIndices),]
    rm(problemSNPs, badIndices)
  }

  # filter 'good' SNPs (e.g. SNP6 positions)
  if((useLociFile != "NA") & (!is.na(useLociFile))) {
    goodSNPs = read.table(useLociFile, header=T, sep="\t", stringsAsFactors=F)
    goodSNPs = goodSNPs$pos[goodSNPs$chr==chrom_name]
    len = length(goodSNPs)
    goodIndices = match(known_SNPs$position, goodSNPs)
    known_SNPs = known_SNPs[!is.na(goodIndices),]
    rm(goodSNPs, goodIndices)
  }

  # Read in the allele counts and see which known SNPs are covered
  snp_data = read.table(tumour.allele.counts.file, comment.char="#", sep="\t", header=F, stringsAsFactors=F)
  normal_snp_data = read.table(normal.allele.counts.file, comment.char="#", sep="\t", header=F, stringsAsFactors=F)
  snp_data = cbind(snp_data, normal_snp_data)
  indices = match(known_SNPs$position, snp_data[,2])
  found_snp_data = snp_data[indices[!is.na(indices)],]
  rm(snp_data)

  # Obtain BAF for this chromosome (note: this is quicker than reading in the whole genome BAF file generated in the earlier step)
  nucleotides = c("A","C","G","T")
  ref_indices = match(known_SNPs[!is.na(indices),3], nucleotides)+ncol(normal_snp_data)+2
  alt_indices = match(known_SNPs[!is.na(indices),4], nucleotides)+ncol(normal_snp_data)+2
  BAFs = as.numeric(found_snp_data[cbind(1:nrow(found_snp_data),alt_indices)])/(as.numeric(found_snp_data[cbind(1:nrow(found_snp_data),alt_indices)])+as.numeric(found_snp_data[cbind(1:nrow(found_snp_data),ref_indices)]))
  BAFs[is.nan(BAFs)] = 0
  rm(nucleotides, ref_indices, alt_indices, found_snp_data, normal_snp_data)

  # Set the minimum level to use for obtaining genotypes
  minBaf = min(heterozygousFilter, 1.0-heterozygousFilter)
  maxBaf = max(heterozygousFilter, 1.0-heterozygousFilter)

  # Obtain genotypes that impute2 is able to understand
  genotypes = array(0,c(sum(!is.na(indices)),3))
  genotypes[BAFs<=minBaf,1] = 1
  genotypes[BAFs>minBaf & BAFs<maxBaf,2] = 1
  genotypes[BAFs>=maxBaf,3] = 1

  # Create the output
  snp.names = paste("snp",1:sum(!is.na(indices)), sep="")
  out.data = cbind(snp.names, known_SNPs[!is.na(indices),1:4], genotypes)

  write.table(out.data, file=output.file, row.names=F, col.names=F, quote=F)
  if(is.na(chrom_name)) {
    sample.g.file = paste(dirname(output.file), "/sample_g.txt", sep="")
    #not sure this is necessary, because only the PAR regions are used for males
    #if(is.male){
    #	sample_g_data=data.frame(ID_1=c(0,"INDIVI1"),ID_2=c(0,"INDIVI1"),missing=c(0,0),sex=c("D",1))
    #}else{
    sample_g_data = data.frame(ID_1=c(0,"INDIVI1"), ID_2=c(0,"INDIVI1"), missing=c(0,0), sex=c("D",2))
    #}
    write.table(sample_g_data, file=sample.g.file, row.names=F, col.names=T, quote=F)
  }
}

#' Function to correct LogR for waivyness that correlates with GC content
#' @param Tumour_LogR_file String pointing to the tumour LogR output
#' @param outfile String pointing to where the GC corrected LogR should be written
#' @param correlations_outfile File where correlations are to be saved
#' @param gc_content_file_prefix String pointing to where GC windows for this reference genome can be
#' found. These files should be split per chromosome and this prefix must contain the full path until
#' chr in its name. The .txt extension is automatically added.
#' @param replic_timing_file_prefix Like the gc_content_file_prefix, containing replication timing info (supply NULL if no replication timing correction is to be applied)
#' @param chrom_names A vector containing chromosome names to be considered
#' @param recalc_corr_afterwards Set to TRUE to recalculate correlations after correction
#' @author jdemeul, sd11
#' @export
gc.correct.wgs = function(Tumour_LogR_file, outfile, correlations_outfile, gc_content_file_prefix, replic_timing_file_prefix, chrom_names, recalc_corr_afterwards=F) {

  if (is.null(gc_content_file_prefix)) {
    stop("GC content reference files must be supplied to WGS GC content correction")
  }

  Tumor_LogR = read_logr(Tumour_LogR_file)

  print("Processing GC content data")
  gc_files = paste0(gc_content_file_prefix, chrom_names, ".txt.gz")
  GC_data = do.call(rbind, lapply(gc_files, read_gccontent))
  colnames(GC_data) = c("chr", "Position", paste0(c(25,50,100,200,500), "bp"),
                        paste0(c(1,2,5,10,20,50,100), "kb"))#,200,500), "kb"),
                        # paste0(c(1,2,5,10), "Mb"))

  if (!is.null(replic_timing_file_prefix)) {
    print("Processing replication timing data")
    replic_files = paste0(replic_timing_file_prefix, chrom_names, ".txt.gz")
    replic_data = do.call(rbind, lapply(replic_files, read_replication))
  }

  # omit non-matching loci, replication data generated at exactly same GC loci
  locimatches = match(x = paste0(Tumor_LogR$Chromosome, "_", Tumor_LogR$Position),
                       table = paste0(GC_data$chr, "_", GC_data$Position))
  Tumor_LogR = Tumor_LogR[which(!is.na(locimatches)), ]
  GC_data = GC_data[na.omit(locimatches), ]
  if (!is.null(replic_timing_file_prefix)) {
    replic_data = replic_data[na.omit(locimatches), ]
  }
  rm(locimatches)

  corr = abs(cor(GC_data[, 3:ncol(GC_data)], Tumor_LogR[,3], use="complete.obs")[,1])
  if (!is.null(replic_timing_file_prefix)) {
    corr_rep = abs(cor(replic_data[, 3:ncol(replic_data)], Tumor_LogR[,3], use="complete.obs")[,1])
  }

  index_1kb = which(names(corr)=="1kb")
  maxGCcol_insert = names(which.max(corr[1:index_1kb]))
  index_100kb = which(names(corr)=="100kb")
  # start large window sizes at 5kb rather than 2kb to avoid overly correlated expl variables
  maxGCcol_amplic = names(which.max(corr[(index_1kb+2):index_100kb]))
  if (!is.null(replic_timing_file_prefix)) {
    maxreplic = names(which.max(corr_rep))
  }

  if (!is.null(replic_timing_file_prefix)) {
    cat("Replication timing correlation: ",paste(names(corr_rep),format(corr_rep,digits=2), ";"),"\n")
    cat("Replication dataset: " ,maxreplic,"\n")
  }
  cat("GC correlation: ",paste(names(corr),format(corr,digits=2), ";"),"\n")
  cat("Short window size: ",maxGCcol_insert,"\n")
  cat("Long window size: ",maxGCcol_amplic,"\n")

  if (!is.null(replic_timing_file_prefix)) {
    # Multiple regression - with replication timing
    corrdata = data.frame(logr = Tumor_LogR[,3, drop = T],
                          GC_insert = GC_data[,maxGCcol_insert, drop = T],
                          GC_amplic = GC_data[,maxGCcol_amplic, drop = T],
                          replic = replic_data[, maxreplic, drop = T])
    colnames(corrdata) = c("logr", "GC_insert", "GC_amplic", "replic")
    if (!recalc_corr_afterwards)
      rm(GC_data, replic_data)

    model = lm(logr ~ splines::ns(x = GC_insert, df = 5, intercept = T) + splines::ns(x = GC_amplic, df = 5, intercept = T) + splines::ns(x = replic, df = 5, intercept = T), y=F, model = F, data = corrdata, na.action="na.exclude")

    corr = data.frame(windowsize=c(names(corr), names(corr_rep)), correlation=c(corr, corr_rep))
    write.table(corr, file=gsub(".txt", "_beforeCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)

  } else {
    # Multiple regression  - without replication timing
    corrdata = data.frame(logr = Tumor_LogR[,3, drop = T],
                          GC_insert = GC_data[,maxGCcol_insert, drop = T],
                          GC_amplic = GC_data[,maxGCcol_amplic, drop = T])
    colnames(corrdata) = c("logr", "GC_insert", "GC_amplic")
    if (!recalc_corr_afterwards)
      rm(GC_data)

    model = lm(logr ~ splines::ns(x = GC_insert, df = 5, intercept = T) + splines::ns(x = GC_amplic, df = 5, intercept = T), y=F, model = F, data = corrdata, na.action="na.exclude")

    corr = data.frame(windowsize=names(corr), correlation=corr)
    write.table(corr, file=gsub(".txt", "_beforeCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)
  }

  Tumor_LogR[,3] = residuals(model)
  rm(model, corrdata)

  readr::write_tsv(x=Tumor_LogR[which(!is.na(Tumor_LogR[,3])), ], path=outfile)

  if (recalc_corr_afterwards) {
    # Recalculate the correlations to see how much there is left
    corr = abs(cor(GC_data[, 3:ncol(GC_data)], Tumor_LogR[,3], use="complete.obs")[,1])
    if (!is.null(replic_timing_file_prefix)) {
      corr_rep = abs(cor(replic_data[, 3:ncol(replic_data)], Tumor_LogR[,3], use="complete.obs")[,1])
      cat("Replication timing correlation post correction: ",paste(names(corr_rep),format(corr_rep,digits=2), ";"),"\n")
    }
    cat("GC correlation post correction: ",paste(names(corr),format(corr,digits=2), ";"),"\n")

    if (!is.null(replic_timing_file_prefix)) {
      corr = data.frame(windowsize=c(names(corr), names(corr_rep)), correlation=c(corr, corr_rep))
      write.table(corr, file=gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)
    } else {
      corr = data.frame(windowsize=c(names(corr)), correlation=corr)
      write.table(corr, file=gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)
    }
  } else {
    corr$correlation = NA
    write.table(corr, file=gsub(".txt", "_afterCorrection.txt", correlations_outfile), sep="\t", quote=F, row.names=F)
  }
}


#' Prepare WGS data for haplotype construction
#'
#' This function performs part of the Battenberg WGS pipeline: Counting alleles, constructing BAF and logR
#' and performing GC content correction.
#'
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param tumourbam Full path to the tumour BAM file
#' @param normalbam Full path to the normal BAM file
#' @param tumourname Identifier to be used for tumour output files
#' @param normalname Identifier to be used for normal output files
#' @param g1000allelesprefix Prefix path to the 1000 Genomes alleles reference files
#' @param g1000prefix Prefix path to the 1000 Genomes SNP reference files
#' @param gccorrectprefix Prefix path to GC content reference data
#' @param repliccorrectprefix Prefix path to replication timing reference data (supply NULL if no replication timing correction is to be applied)
#' @param min_base_qual Minimum base quality required for a read to be counted
#' @param min_map_qual Minimum mapping quality required for a read to be counted
#' @param allelecounter_exe Path to the allele counter executable (can be found in $PATH)
#' @param min_normal_depth Minimum depth required in the normal for a SNP to be included
#' @param nthreads The number of paralel processes to run
#' @param skip_allele_counting Flag, set to TRUE if allele counting is already complete (files are expected in the working directory on disk)
#' @param skip_allele_counting_normal Flag, set to TRUE from the second sample onwards for multisample case (Default: FALSE)
#' @author sd11
#' @export
prepare_wgs = function(chrom_names, tumourbam, normalbam, tumourname, normalname, g1000allelesprefix, g1000prefix, gccorrectprefix,
                       repliccorrectprefix, min_base_qual, min_map_qual, allelecounter_exe, min_normal_depth, nthreads, skip_allele_counting, skip_allele_counting_normal = F) {

  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")

  if (!skip_allele_counting) {
    # Obtain allele counts for 1000 Genomes locations for both tumour and normal
    foreach::foreach(i=1:length(chrom_names)) %dopar% {
      getAlleleCounts(bam.file=tumourbam,
                      output.file=paste(tumourname,"_alleleFrequencies_chr", chrom_names[i], ".txt", sep=""),
                      g1000.loci=paste(g1000prefix, chrom_names[i], ".txt", sep=""),
                      min.base.qual=min_base_qual,
                      min.map.qual=min_map_qual,
                      allelecounter.exe=allelecounter_exe)
      
      if (!skip_allele_counting_normal) {
        getAlleleCounts(bam.file=normalbam,
                        output.file=paste(normalname,"_alleleFrequencies_chr", chrom_names[i], ".txt",  sep=""),
                        g1000.loci=paste(g1000prefix, chrom_names[i], ".txt", sep=""),
                        min.base.qual=min_base_qual,
                        min.map.qual=min_map_qual,
                        allelecounter.exe=allelecounter_exe)
      }
    }
  }

  # Obtain BAF and LogR from the raw allele counts
  getBAFsAndLogRs(tumourAlleleCountsFile.prefix=paste(tumourname,"_alleleFrequencies_chr", sep=""),
                  normalAlleleCountsFile.prefix=paste(normalname,"_alleleFrequencies_chr", sep=""),
                  figuresFile.prefix=paste(tumourname, "_", sep=''),
                  BAFnormalFile=paste(tumourname,"_normalBAF.tab", sep=""),
                  BAFmutantFile=paste(tumourname,"_mutantBAF.tab", sep=""),
                  logRnormalFile=paste(tumourname,"_normalLogR.tab", sep=""),
                  logRmutantFile=paste(tumourname,"_mutantLogR.tab", sep=""),
                  combinedAlleleCountsFile=paste(tumourname,"_alleleCounts.tab", sep=""),
                  chr_names=chrom_names,
                  g1000file.prefix=g1000allelesprefix,
                  minCounts=min_normal_depth,
                  samplename=tumourname)
  # Perform GC correction
  gc.correct.wgs(Tumour_LogR_file=paste(tumourname,"_mutantLogR.tab", sep=""),
                 outfile=paste(tumourname,"_mutantLogR_gcCorrected.tab", sep=""),
                 correlations_outfile=paste(tumourname, "_GCwindowCorrelations.txt", sep=""),
                 gc_content_file_prefix=gccorrectprefix,
                 replic_timing_file_prefix=repliccorrectprefix,
                 chrom_names=chrom_names)
}
