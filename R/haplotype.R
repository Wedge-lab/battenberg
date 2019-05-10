#' Morphs phased SNPs from SNP6 input into haplotype blocks
#' 
#' This function matches allele frequencies and halplotype info, reverses frequencies by haplotype, combines the output and saves it to disk.
#' @param chrom The chromosome number for which this function should run.
#' @param alleleFreqFile File containing allele frequency information.
#' @param haplotypeFile File containing haplotype information.
#' @param samplename Identifier for the sample (used in header of output file).
#' @param outputfile Full path pointing to where output should be written.
#' @param chr_names Vector of chromosome names
#' @author dw9
#' @export
GetChromosomeBAFs_SNP6 = function(chrom, alleleFreqFile, haplotypeFile, samplename, outputfile, chr_names) {
  # Read in the allele frequencies and variant info
  alleleFreqData = read.csv(alleleFreqFile, header=T)
  variant_data = read.table(haplotypeFile, header=F)
  
  # TODO: Check columns input
  
  # Match the two
  alleleFreqData = alleleFreqData[alleleFreqData[,1] %in% variant_data[,3],]
  select = match(alleleFreqData[,1], variant_data[,3])
  variant_data = variant_data[select,]
  
  chr_name = chr_names[as.numeric(chrom)]
  print(chr_name)

  # Switch the haplotypes where required
  alleleFreqs = alleleFreqData$allele.frequency
  reversedHaplotypes = variant_data[,6]==1
  alleleFreqs[reversedHaplotypes] = 1.0-alleleFreqs[reversedHaplotypes]
  
  print(paste(nrow(variant_data),length(alleleFreqs),sep=","))
  # Combine the allele frequencies and variant info and save output
  knownMutBAFs = cbind(chr_name,variant_data[,3],alleleFreqs)
  write.table(knownMutBAFs, outputfile, sep="\t", row.names=F, col.names=c("Chromosome", "Position", samplename), quote=F)	
}

#' Morphs phased SNPs from WGS input into haplotype blocks
#' 
#' @param chrom The chromosome number for which this function is called.
#' @param SNP_file File containing allele counts for each SNP location.
#' @param haplotypeFile File containing impute phasing output.
#' @param samplename Name of the sample (used in header of output file).
#' @param outfile Full path to where the output will be written.
#' @param chr_names Names of all allowed chromosomes as a Vector.
#' @param minCounts An integer describing the minimum number of reads covering this position to be included in the output.
#' @author dw9
#' @export
GetChromosomeBAFs = function(chrom, SNP_file, haplotypeFile, samplename, outfile, chr_names, minCounts=1) {
  # Read in the SNP and haplotype info
  snp_data = read.table(SNP_file, comment.char="", sep="\t", header=T, stringsAsFactors=F)
  variant_data = read.table(haplotypeFile, header=F)
  
  # TODO: Check columns input
  
  print(snp_data[1:3,])
  print(chr_names)
  print(chrom)
  
  # Just select heterozygous SNPs
  het_variant_data = variant_data[variant_data[,6] != variant_data[,7],]
  
  chr_name = chr_names[as.numeric(chrom)]
  print(chr_name)
  
  # Match allele counts and phasing info
  indices = match(het_variant_data[,3],snp_data[,2])
  het_variant_data = het_variant_data[!is.na(indices),]
  snp_indices = indices[!is.na(indices)]
  filtered_snp_data = snp_data[snp_indices,]
  
  # No matches found, save empty file and quit
  if(nrow(het_variant_data)==0 | is.null(het_variant_data)) {
    write.table(array(NA,c(0,3)),outfile,sep="\t",col.names=c("Chromosome","Position",samplename),quote=F,row.names=F)		
    return()
  }
  print(filtered_snp_data[1:3,])
  
  # Decode 1,2,3,4 to A,C,G,T (encoding used in the variant_data input files)
  # TODO: place this in utils script? Isn't this also performed in GenerateImputeInputFromAlleleFrequencies.R?
  nucleotides=c("A","C","G","T")
  ref_indices = match(het_variant_data[cbind(1:nrow(het_variant_data),4+het_variant_data[,6])],nucleotides)
  alt_indices = match(het_variant_data[cbind(1:nrow(het_variant_data),4+het_variant_data[,7])],nucleotides)
  
  # Obtain counts for both alleles and the total
  ref.count = as.numeric(filtered_snp_data[cbind(1:nrow(filtered_snp_data),alt_indices+2)])
  alt.count = as.numeric(filtered_snp_data[cbind(1:nrow(filtered_snp_data),ref_indices+2)])
  denom = ref.count+alt.count

  # Filter out those SNPs that have less than minCounts reads
  min_indices = denom>=minCounts
  filtered_snp_data = filtered_snp_data[min_indices,]
  denom = denom[min_indices]
  alt.count = alt.count[min_indices]
  
  # No matches found, save empty file and quit
  if(nrow(filtered_snp_data)==0 | is.null(filtered_snp_data)) {
    write.table(array(NA,c(0,3)),outfile,sep="\t",col.names=c("Chromosome","Position",samplename),quote=F,row.names=F)  	
    return()
  }
  
  # Save all to disk
  hetMutBAFs = cbind(chr_name,filtered_snp_data[,2],alt.count/denom)
  write.table(hetMutBAFs,outfile,sep="\t",row.names=F,col.names=c("Chromosome","Position",samplename),quote=F)
}

#' Plot haplotyped SNPs
#' 
#' This function takes haplotyped SNPs and plots them to a png file.
#' @param haplotyped.baf.file File containing the haplotyped SNP info.
#' @param imageFileName Filename as which the png will be saved.
#' @param samplename Name of the sample to be used in image title.
#' @param chrom The chromosome that is plotted.
#' @param chr_names A list of allowed chromosome names.
#' @author dw9
#' @export
plot.haplotype.data = function(haplotyped.baf.file, imageFileName, samplename, chrom, chr_names) {
  chr_name = chr_names[chrom]
  mut_data = read.table(haplotyped.baf.file,sep="\t",header=T)
  
  if (nrow(mut_data) > 0) {
    x_min = min(mut_data$Position,na.rm=T)
    x_max = max(mut_data$Position,na.rm=T)
  } else {
    x_min = 1
    x_max = 2
  }

  png(filename = imageFileName, width = 10000, height = 2500, res = 500)
  create.haplotype.plot(chrom.position=mut_data$Position, 
                        points.blue=mut_data[,3], 
                        points.red=1-mut_data[,3], 
                        x.min=x_min, 
                        x.max=x_max, 
                        title=paste(samplename,", chromosome",chr_name, sep=" "), 
                        xlab="pos", 
                        ylab="BAF")
  dev.off()
}

#' Combines all separate BAF files per chromosome into a single file
#'
#' @param inputfile.prefix Prefix of the input files until the chromosome number. The chromosome number will be added internally
#' @param inputfile.postfix Postfix of the input files from the chromosome number
#' @param outputfile Full path to where the output will be written
#' @param no.chrs Number of chromosomes, i.e. the number of files that need to be concatenated
#' @author dw9
#' @export
combine.baf.files = function(inputfile.prefix, inputfile.postfix, outputfile, no.chrs) {
  concatenateBAFfiles(inputfile.prefix, inputfile.postfix, outputfile, no.chrs)
}
