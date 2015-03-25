#' Morphs phased SNPs from SNP6 input into haplotype blocks
#' 
#' This function matches allele frequencies and halplotype info, reverses frequencies by haplotype, combines the output and saves it to disk.
#' @param chrom The chromosome number for which this function should run.
#' @param alleleFreqFile File containing allele frequency information.
#' @param haplotypeFile File containing haplotype information.
#' @param samplename Identifier for the sample (used in header of output file).
#' @param outputfile Full path pointing to where output should be written.
#' @author dw9
#' @export
GetChromosomeBAFs_SNP6 = function(chrom, alleleFreqFile, haplotypeFile, samplename, outputfile, chr_names) {
  # Read in the allele frequencies and variant info
  alleleFreqData = read.csv(alleleFreqFile, header=T)
  variant_data = read.table(haplotypeFile, header=F)
  
  # TODO: Check columns input
  
  # Match the two
  alleleFreqData = alleleFreqData[alleleFreqData[,1] %in% variant_data[,3],]
  
#   chr_name = toString(chrom)
#   if(chr_name == "23"){
#     chr_name = "X"
#   }
  chr_name = chr_names[as.numeric(chrom)]
  print(chr_name)

  # Switch the haplotypes where required
  alleleFreqs = alleleFreqData$allele.frequency
  reversedHaplotypes = variant_data[,6]==1
  alleleFreqs[reversedHaplotypes] = 1.0-alleleFreqs[reversedHaplotypes]
  
  print(paste(nrow(variant_data),length(alleleFreqs),sep=","))
  # Combine the allele frequencies and variant info and save output
  knownMutBAFs = cbind(chr_name,variant_data[,3],alleleFreqs)
  write.table(knownMutBAFs,outputfile,sep="\t",row.names=paste("snp",1:nrow(knownMutBAFs),sep=""),col.names=c("Chromosome","Position",samplename),quote=F)	
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
GetChromosomeBAFs = function(chrom, SNP_file, haplotypeFile, samplename, outfile, chr_names, minCounts=10) {
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
  if(nrow(het_variant_data)==0){
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
  
  # Save all to disk
  hetMutBAFs = cbind(chr_name,filtered_snp_data[,2],alt.count/denom)
  write.table(hetMutBAFs,outfile,sep="\t",row.names=F,col.names=c("Chromosome","Position",samplename),quote=F)
}

#' Plot haplotyped SNPs
#' 
#' This function takes haplotyped SNPs and plots them to a png file.
#' @param mutFile File containing the haplotyped SNP info.
#' @param imageFileName Filename as which the png will be saved.
#' @param samplename Name of the sample to be used in image title.
#' @param chrom The chromosome that is plotted.
#' @param chr_names A list of allowed chromosome names.
#' @author dw9
#' @export
plot.haplotype.data = function(haplotyped.baf.file, imageFileName, samplename, chrom, chr_names) {
  #impute.info = read.table(imputeInfoFile,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
  # TODO: is this really required?
  #chr_names = unique(parse.imputeinfofile(imputeinfofile, FALSE)$chrom)
  chr_name = chr_names[chrom]
  #imageFileName = paste(samplename,"_chr",chr_name,"_heterozygousData.png",sep="")
  
  mut_data = read.table(haplotyped.baf.file,sep="\t",header=T)
  
  png(filename = imageFileName, width = 10000, height = 2500, res = 500)
  # TODO: This should move to the plotting script
#   par(pch = ".", cex=1, cex.main=0.8, cex.axis = 0.6, cex.lab=0.7,yaxp=c(-0.05,1.05,6))
#   plot(c(min(mut_data$Position,na.rm=T),max(mut_data$Position,na.rm=T)),c(0,1),type="n",,main=paste(samplename,", chromosome",mut_data[1,1],sep=" "),xlab="pos", ylab="BAF")
#   points(mut_data$Position,mut_data[,3],col="blue")
#   points(mut_data$Position,1-mut_data[,3],col="red")
  create.haplotype.plot(chrom.position=mut_data$Position, 
                        points.blue=mut_data[,3], 
                        points.red=1-mut_data[,3], 
                        x.min=min(mut_data$Position,na.rm=T), 
                        x.max=max(mut_data$Position,na.rm=T), 
                        title=paste(samplename,", chromosome",mut_data[1,1], sep=" "), 
                        xlab="pos", 
                        ylab="BAF")
  dev.off()
}

#' Combines all separate BAF files per chromosome into a single file
#'
#' @param inputfile.prefix Prefix of the input files until the chromosome number. The chromosome number will be added internally.
#' @param inputfile.postfix Postfix of the input files from the chromosome number.
#' @param outputfile Full path to where the output will be written.
#' @param no.chrs Number of chromosomes.
#' @author dw9
#' @export
combine.baf.files = function(inputfile.prefix, inputfile.postfix, outputfile, no.chrs) {
  # TODO: Possibly deal with ignored chromosomes here
#   impute.info = read.table(impute_info_file,header=F,row.names=NULL,sep="\t",stringsAsFactors=F)
#   if(is.male){
#     impute.info = impute.info[impute.info[,7]==1,]
#   }
#   no.chrs = length(unique(impute.info[,1]))
  #no.chrs = length(unique(parse.imputeinfofile(imputeinfofile, FALSE)$chrom))
  concatenateBAFfiles(inputfile.prefix, inputfile.postfix, outputfile, no.chrs)
}

#' Segment the haplotyped and phased data using fastPCF.
#' @param
#' @author dw9
segment.baf.phased = function(samplename, inputfile, outputfile, gamma=10, phasegamma=3, kmin=3, phasekmin=3) {
  #BAFraw = read.table(paste(sample,"_allChromosomes_heterozygousMutBAFs_haplotyped.txt",sep=""),sep="\t",header=T)
  BAFraw = read.table(inputfile,sep="\t",header=T, stringsAsFactors=F)
  
  BAFoutput = NULL
  
  for (chr in unique(BAFraw[,1])) {
    BAFrawchr = BAFraw[BAFraw[,1]==chr,c(2,3)]
    BAFrawchr = BAFrawchr[!is.na(BAFrawchr[,2]),]
    
    BAF = BAFrawchr[,2]
    pos = BAFrawchr[,1]
    names(BAF) = rownames(BAFrawchr)
    names(pos) = rownames(BAFrawchr)
    
    sdev <- getMad(ifelse(BAF<0.5,BAF,1-BAF),k=25)
    #DCW 250314
    #for cell lines, sdev goes to zero in regions of LOH, which causes problems.
    #0.09 is around the value expected for a binomial distribution around 0.5 with depth 30
    if(sdev<0.09){
      sdev = 0.09
    }
    
    print(paste("BAFlen=",length(BAF),sep=""))
    if(length(BAF)<50){
      BAFsegm = rep(mean(BAF),length(BAF))
    }else{
      #res= selectFastPcf(BAF,3,3*sdev,T)
      res= selectFastPcf(BAF,phasekmin,phasegamma*sdev,T)
      BAFsegm = res$yhat
    }
    
    png(filename = paste(samplename,"_RAFseg_chr",chr,".png",sep=""), width = 2000, height = 1000, res = 200)
    #     par(mar = c(5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, cex.lab = 2)
    #     plot(c(min(pos)/1000000,max(pos)/1000000),c(0,1),pch=".",type = "n", 
    #          main = paste(sample,", chromosome ", chr, sep=""), xlab = "Position (Mb)", ylab = "BAF (phased)")
    #     points(pos/1000000,BAF,pch=".",col="red",cex=2)
    #     points(pos/1000000,BAFsegm,pch=19,cex=0.5,col="green")
    create.segmented.plot(chrom.position=pos/1000000, 
                                  points.red=BAF, 
                                  points.green=BAFsegm, 
                                  x.min=min(pos)/1000000, 
                                  x.max=max(pos)/1000000, 
                                  title=paste(samplename,", chromosome ", chr, sep=""), 
                                  xlab="Position (Mb)", 
                                  ylab="BAF (phased)")
    dev.off()
    
    BAFphased = ifelse(BAFsegm>0.5,BAF,1-BAF)
    
    if(length(BAFphased)<50){
      BAFphseg = rep(mean(BAFphased),length(BAFphased))
    }else{
      res = selectFastPcf(BAFphased,kmin,gamma*sdev,T)
      BAFphseg = res$yhat
    }
    
    png(filename = paste(samplename,"_chr",chr,".png",sep=""), width = 2000, height = 1000, res = 200)
    #     par(mar = c(5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, cex.lab = 2)
    #     plot(c(min(pos)/1000000,max(pos)/1000000),c(0,1),pch=".",type = "n", 
    #          main = paste(sample,", chromosome ", chr, sep=""), xlab = "Position (Mb)", ylab = "BAF (phased)")
    #     points(pos/1000000,BAF,pch=".",col=ifelse(BAFsegm>0.5,"red","blue"),cex=2)
    #     points(pos/1000000,BAFphseg,pch=19,cex=0.5,col="darkred")
    #     points(pos/1000000,1-BAFphseg,pch=19,cex=0.5,col="darkblue")
    create.baf.plot(chrom.position=pos/1000000, 
                    points.red.blue=BAF, 
                    points.darkred=BAFphseg, 
                    points.darkblue=1-BAFphseg, 
                    x.min=min(pos)/1000000, 
                    x.max=max(pos)/1000000, 
                    title=paste(samplename,", chromosome ", chr, sep=""), 
                    xlab="Position (Mb)", 
                    ylab="BAF (phased)")
    dev.off()
    
    BAFphased = ifelse(BAFsegm>0.5, BAF, 1-BAF)
    BAFoutputchr = cbind(rep(chr, length(BAFphseg)), pos, BAF, BAFphased, BAFphseg)
    BAFoutput = rbind(BAFoutput, BAFoutputchr)
  }
  colnames(BAFoutput) = c("Chromosome","Position","BAF","BAFphased","BAFseg")
  #write.table(BAFoutput,paste(sample,".BAFsegmented.txt",sep=""),sep="\t",row.names=T,col.names=NA,quote=F)
  write.table(BAFoutput, outputfile, sep="\t", row.names=F, col.names=T, quote=F)
}
