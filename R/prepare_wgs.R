
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
  
  system(cmd, wait=T)
}


#' Obtain BAF and LogR from the allele counts
#'
#' @param outputdir Directory where ASCAT plots should be written.
#' @param tumourAlleleCountsFile.prefix Prefix of the allele counts files for the tumour.
#' @param normalAlleleCountsFile.prefix Prefix of the allele counts files for the normal.
#' @param BAFnormalFile File where BAF from the normal will be written.
#' @param BAFmutantFile File where BAF from the tumour will be written.
#' @param logRnormalFile File where LogR from the normal will be written.
#' @param logRmutantFile File where LogR from the tumour will be written.
#' @param combinedAlleleCountsFile File where combined allele counts for tumour and normal will be written. 
#' @param chr_names A vector with allowed chromosome names.
#' @param g1000file.prefix Prefix to where 1000 Genomes reference files can be found.
#' @param minCounts Integer, minimum depth required for a SNP to be included (optional, default=NA).
#' @param samplename String, name of the sample (optional, default=sample1).
#' @author dw9, sd11
#' @export
getBAFsAndLogRs = function(tumourAlleleCountsFile.prefix, normalAlleleCountsFile.prefix, figuresFile.prefix, BAFnormalFile, BAFmutantFile, logRnormalFile, logRmutantFile, combinedAlleleCountsFile, chr_names, g1000file.prefix, minCounts=NA, samplename="sample1") {
  
  input_data = concatenateAlleleCountFiles(tumourAlleleCountsFile.prefix, ".txt", length(chr_names))
  normal_input_data = concatenateAlleleCountFiles(normalAlleleCountsFile.prefix, ".txt", length(chr_names))
  allele_data = concatenateG1000SnpFiles(g1000file.prefix, ".txt", length(chr_names))
  
  #save(file=paste(samplename, "_BAFLogR_concatenated_input.RData", sep=""), input_data, normal_input_data, allele_data)
  #load(paste(samplename, "_BAFLogR_concatenated_input.RData", sep=""))
  
  names(input_data)[1] = "CHR"
  names(normal_input_data)[1] = "CHR"
  
  normal_data = normal_input_data[,3:6]
  mutant_data = input_data[,3:6]
  
  # Obtain depth for both alleles for tumour and normal
  len = nrow(normal_data)
  normCount1 = normal_data[cbind(1:len,allele_data[,2])]
  normCount2 = normal_data[cbind(1:len,allele_data[,3])]
  totalNormal = normCount1 + normCount2
  mutCount1 = mutant_data[cbind(1:len,allele_data[,2])]
  mutCount2 = mutant_data[cbind(1:len,allele_data[,3])]
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
    # These variables are cleaned
    #normal_data = normal_data[indices,]
    #mutant_data = mutant_data[indices,]
    
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
  #ascat.bc = ascat.loadData(logRmutantFile, BAFmutantFile, logRnormalFile, BAFnormalFile)
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
  # TODO: On PD7422a this produces different plots than before (see streak on chrom 14)
  ascat.plotRawData(ascat.bc) #, parentDir=figuresFile.prefix)
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
  chrom_name = parse.imputeinfofile(imputeinfofile, is.male)$chrom[chrom]
  
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
  if(is.na(as.numeric(chrom_name))) {
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
#' @param Tumour_LogR_file: String pointing to the tumour LogR output
#' @param outfile: String pointing to where the GC corrected LogR should be written
#' @param gc_content_file_prefix: String pointing to where GC windows for this reference genome can be 
#' found. These files should be split per chromosome and this prefix must contain the full path until
#' chr in its name. The .txt extension is automatically added.
#' @param chrom_names: A vector containing chromosome names to be considered
#' @author tjm
gc.correct.wgs = function(Tumour_LogR_file, outfile, gc_content_file_prefix, chrom_names) {

  Tumor_LogR_new = NULL
  Tumor_LogR = read.table(Tumour_LogR_file, stringsAsFactors=F, header=T)

  GC_data = NULL
  for (chrindex in 1:length(chrom_names)) {
    print(paste("chr =", chrindex))
    Tumor_LogR_chr = Tumor_LogR[Tumor_LogR$Chromosome %in% chrindex,]
    GC_newlist = read.table(paste(gc_content_file_prefix, chrindex, ".txt", sep=""), header=T, stringsAsFactors=F)
    colnames(GC_newlist)[c(1,2)] = c("Chr","Position")
    GC_newlist = GC_newlist[GC_newlist$Position %in% Tumor_LogR_chr$Position,]
    GC_data = rbind(GC_data, GC_newlist)
    Tumor_LogR_new = rbind(Tumor_LogR_new, Tumor_LogR_chr[!is.na(match(Tumor_LogR_chr$Position, GC_newlist$Position)),])
  }
  Tumor_LogR = Tumor_LogR_new  
  rm(Tumor_LogR_new)
  
#     ovl = which(Tumor_LogR_chr$Position %in% GC_newlist$Position)
#     Tumor_LogR_chr = Tumor_LogR_chr[ovl,,drop=F]
  
  flag_nona = is.finite(Tumor_LogR[,3])
  corr = cor(GC_data[flag_nona, 3:ncol(GC_data)], Tumor_LogR[flag_nona,3], use="complete.obs")
  length = nrow(Tumor_LogR)
  
  corr = apply(corr, 1, function(x) sum(abs(x*length))/sum(length))
  index_1M = c(which(names(corr)=="X1M"), which(names(corr)=="X1Mb"))
  maxGCcol_short = which(corr[1:(index_1M-1)]==max(corr[1:(index_1M-1)]))
  maxGCcol_long = which(corr[index_1M:length(corr)]==max(corr[index_1M:length(corr)]))
  maxGCcol_long = maxGCcol_long+(index_1M-1)
  
  cat("weighted correlation: ",paste(names(corr),format(corr,digits=2), ";"),"\n")   
  cat("Short window size: ",names(GC_data)[maxGCcol_short+2],"\n")
  cat("Long window size: ",names(GC_data)[maxGCcol_long+2],"\n")
  
  # Multiple regression 
  flag_NA = (is.infinite(Tumor_LogR[,3])) | (is.na(GC_data[,2+maxGCcol_short])) | (is.na(GC_data[,2+maxGCcol_long]))
  td_select = Tumor_LogR[!flag_NA,3]
  GC_data_select = GC_data[!flag_NA,]
  x1 = GC_data_select[,2+maxGCcol_short]
  x2 = GC_data_select[,2+maxGCcol_long]
  x3 = (x1)^2
  x4 = (x2)^2
  model = lm(td_select ~ x1+x2+x3+x4, y=T, na.action="na.omit")
  
  GCcorrected = Tumor_LogR[!flag_NA,]
  GCcorrected[,3] = model$residuals
    
    #Tumor_LogR_new = rbind(Tumor_LogR_new,GCcorrected)
  #write.table(Tumor_LogR_new, file=outfile, sep="\t", quote=F)
  write.table(GCcorrected, file=outfile, sep="\t", quote=F)
}
