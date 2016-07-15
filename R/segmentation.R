#' Segment the haplotyped and phased data using fastPCF.
#' 
#' This function performs segmentation. This is done in two steps. First a segmentation step
#' that aims to find short segments. These are used to find haplotype blocks that have been
#' switched. These blocks are switched into the correct order first after which the second
#' segmentation step is performed. This second step aims to segment the data that will go into
#' fit.copy.number. This function produces a BAF segmented file with 5 columns: chromosome, position,
#' original BAF, switched BAF and BAF segment. The BAF segment column should be used subsequently
#' @param samplename Name of the sample, which is used to name output figures
#' @param inputfile String that points to the output from the \code{combine.baf.files} function. This contains the phased SNPs with their BAF values
#' @param outputfile String where the segmentation output will be written
#' @param gamma The gamma parameter controls the size of the penalty of starting a new segment during segmentation. It is therefore the key parameter for controlling the number of segments (Default 10)
#' @param kmin Kmin represents the minimum number of probes/SNPs that a segment should consist of (Default 3)
#' @param phasegamma Gamma parameter used when correcting phasing mistakes (Default 3)
#' @param phasekmin Kmin parameter used when correcting phasing mistakes (Default 3)
#' @author dw9
#' @export
segment.baf.phased = function(samplename, inputfile, outputfile, gamma=10, phasegamma=3, kmin=3, phasekmin=3) {
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
    # Standard deviation is not defined for a single value
    if (is.na(sdev)) {
      sdev = 0
    }
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
      res= selectFastPcf(BAF,phasekmin,phasegamma*sdev,T)
      BAFsegm = res$yhat
    }
    
    png(filename = paste(samplename,"_RAFseg_chr",chr,".png",sep=""), width = 2000, height = 1000, res = 200)
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
    
    png(filename = paste(samplename,"_segment_chr",chr,".png",sep=""), width = 2000, height = 1000, res = 200)
    create.baf.plot(chrom.position=pos/1000000, 
                    points.red.blue=BAF, 
                    plot.red=BAFsegm>0.5,
                    points.darkred=BAFphseg, 
                    points.darkblue=1-BAFphseg, 
                    x.min=min(pos)/1000000, 
                    x.max=max(pos)/1000000, 
                    title=paste(samplename,", chromosome ", chr, sep=""), 
                    xlab="Position (Mb)", 
                    ylab="BAF (phased)")
    dev.off()
    
    BAFphased = ifelse(BAFsegm>0.5, BAF, 1-BAF)
    BAFoutputchr = data.frame(Chromosome=rep(chr, length(BAFphseg)), Position=pos, BAF=BAF, BAFphased=BAFphased, BAFseg=BAFphseg)
    BAFoutput = rbind(BAFoutput, BAFoutputchr)
  }
  colnames(BAFoutput) = c("Chromosome","Position","BAF","BAFphased","BAFseg")
  write.table(BAFoutput, outputfile, sep="\t", row.names=F, col.names=T, quote=F)
}

#' Segment BAF with the inclusion of structural variant breakpoints
#' 
#' This function takes the SV breakpoints as initial segments and runs PCF on each
#' of those independently. The SVs must be supplied as a simple data.frame with columns
#' chromosome and position
#' @param samplename Name of the sample, which is used to name output figures
#' @param inputfile String that points to the output from the \code{combine.baf.files} function. This contains the phased SNPs with their BAF values
#' @param outputfile String where the segmentation output will be written
#' @param svs Data.frame with chromosome and position columns
#' @param gamma The gamma parameter controls the size of the penalty of starting a new segment during segmentation. It is therefore the key parameter for controlling the number of segments (Default 10)
#' @param kmin Kmin represents the minimum number of probes/SNPs that a segment should consist of (Default 3)
#' @param phasegamma Gamma parameter used when correcting phasing mistakes (Default 3)
#' @param phasekmin Kmin parameter used when correcting phasing mistakes (Default 3)
#' @author sd11, dw9
#' @export
segment.baf.phased.sv = function(samplename, inputfile, outputfile, svs, gamma=10, phasegamma=3, kmin=3, phasekmin=3) {
  
  #' Helper function that creates segment breakpoints from SV calls
  #' @param svs_chrom Structural variant breakpoints for a single chromosome
  #' @param BAFrawchr Raw BAF values of germline heterozygous SNPs on a single chromosome
  #' @return A data.frame with chrom, start and end columns
  #' @author sd11
  svs_to_presegment_breakpoints = function(chrom, svs_chrom, BAFrawchr) {
    svs_breakpoints = svs_chrom$position
    
    # If there are no SVs we cannot insert any breakpoints
    if (length(svs_breakpoints) > 0) {
      breakpoints = data.frame()
      
      # check which comes first, the breakpoint or the first SNP
      if (BAFrawchr$Position[1] < svs_breakpoints[1]) {
        startpos = BAFrawchr$Position[1]
        startfromsv = 1 # We're starting from SNP data, so the first SV should be added first
      } else {
        startpos = svs_breakpoints[1]
        startfromsv = 2 # We've just added the first SV, don't use it again
      }
      
      for (svposition in svs_breakpoints[startfromsv:length(svs_breakpoints)]) {
        selectedsnps = BAFrawchr$Position >= startpos & BAFrawchr$Position <= svposition
        if (sum(selectedsnps, na.rm=T) > 0) {
          endindex = max(which(selectedsnps))
          breakpoints = rbind(breakpoints, data.frame(chrom=chrom, start=startpos, end=BAFrawchr$Position[endindex]))
          # Previous SV is the new starting point for the next segment
          startpos = BAFrawchr$Position[endindex + 1]
          print(paste0("Set ", startpos, " to ", BAFrawchr$Position[endindex]))
        }
      }
      
      # Add the remainder of the chromosome, if available
      if (BAFrawchr$Position[nrow(BAFrawchr)] > svs_breakpoints[length(svs_breakpoints)]) {
        endindex = nrow(BAFrawchr)
        breakpoints = rbind(breakpoints, data.frame(chrom=chrom, start=startpos, end=BAFrawchr$Position[endindex]))
        print(paste0("Set ", startpos, " to ", BAFrawchr$Position[endindex]))
      }
    } else {
      # There are no SVs, so create one big segment
      print("No SVs")
      breakpoints = data.frame(chrom=chrom, start=BAFrawchr$Position[1], end=BAFrawchr$Position[nrow(BAFrawchr)])
      print(paste0("Set ", BAFrawchr$Position[1], " to ", BAFrawchr$Position[nrow(BAFrawchr)]))
    }
    return(breakpoints)
  }
  
  #' Run PCF on presegmented data
  #' @param BAFrawchr
  #' @param presegment_chrom_start
  #' @param presegment_chrom_end
  #' @param phasekmin
  #' @param phasegamma
  #' @param kmin
  #' @param gamma
  #' @return A data.frame with columns Chromosome,Position,BAF,BAFphased,BAFseg
  run_pcf = function(BAFrawchr, presegment_chrom_start, presegment_chrom_end, phasekmin, phasegamma, kmin, gamma) {
    row.indices = which(BAFrawchr$Position >= presegment_chrom_start & 
                          BAFrawchr$Position <= presegment_chrom_end)
    
    BAF = BAFrawchr[row.indices,2]
    pos = BAFrawchr[row.indices,1]
    # names(BAF) = rownames(BAFrawchr[row.indices])
    # names(pos) = rownames(BAFrawchr[row.indices])
    
    sdev <- Battenberg:::getMad(ifelse(BAF<0.5,BAF,1-BAF),k=25)
    # Standard deviation is not defined for a single value
    if (is.na(sdev)) {
      sdev = 0
    }
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
      res= Battenberg:::selectFastPcf(BAF,phasekmin,phasegamma*sdev,T)
      BAFsegm = res$yhat
    }
    
    BAFphased = ifelse(BAFsegm>0.5,BAF,1-BAF)
    
    if(length(BAFphased)<50){
      BAFphseg = rep(mean(BAFphased),length(BAFphased))
    }else{
      res = Battenberg:::selectFastPcf(BAFphased,kmin,gamma*sdev,T)
      BAFphseg = res$yhat
    }
    
    return(data.frame(Chromosome=rep(chr, length(row.indices)), 
                      Position=BAFrawchr[row.indices,1], 
                      BAF=BAF, 
                      BAFphased=BAFphased, 
                      BAFseg=BAFphseg,
                      tempBAFsegm=BAFsegm)) # Keep track of BAFsegm for the plot below
  }
  
  # bafsegments, breakpoints, kmin, gamma_param, samplename, filename_suffix="jabba_aspcf"
  BAFraw = read.table(inputfile,sep="\t",header=T, stringsAsFactors=F)
  
  BAFoutput = NULL
  for (chr in unique(BAFraw[,1])) {
    print(paste0("Segmenting ", chr))
    BAFrawchr = BAFraw[BAFraw[,1]==chr,c(2,3)]
    # BAFrawchr = bafsegments[bafsegments$Chromosome==chr, c(2,3)]
    BAFrawchr = BAFrawchr[!is.na(BAFrawchr[,2]),]
    svs_chrom = svs[svs$chromosome==chr,]
    
    breakpoints_chrom = svs_to_presegment_breakpoints(chr, svs_chrom, BAFrawchr)
    print(breakpoints_chrom)
    BAFoutputchr = NULL
    
    for (r in 1:nrow(breakpoints_chrom)) {
      print(breakpoints_chrom[r,])
      BAFoutput_preseg = run_pcf(BAFrawchr, breakpoints_chrom$start[r], breakpoints_chrom$end[r], phasekmin, phasegamma, kmin, gamma)
      BAFoutputchr = rbind(BAFoutputchr, BAFoutput_preseg)
    }
    
    png(filename = paste(samplename,"_RAFseg_chr",chr,".png",sep=""), width = 2000, height = 1000, res = 200)
    create.segmented.plot(chrom.position=BAFoutputchr$Position/1000000, 
                          points.red=BAFoutputchr$BAF, 
                          points.green=BAFoutputchr$tempBAFsegm, 
                          x.min=min(BAFoutputchr$Position)/1000000, 
                          x.max=max(BAFoutputchr$Position)/1000000, 
                          title=paste(samplename,", chromosome ", chr, sep=""), 
                          xlab="Position (Mb)", 
                          ylab="BAF (phased)",
                          svs_pos=svs_chrom$position/1000000)
    dev.off()
    
    png(filename = paste(samplename,"_segment_chr",chr,".png",sep=""), width = 2000, height = 1000, res = 200)
    create.baf.plot(chrom.position=BAFoutputchr$Position/1000000, 
                    points.red.blue=BAFoutputchr$BAF, 
                    plot.red=BAFoutputchr$tempBAFsegm>0.5,
                    points.darkred=BAFoutputchr$BAFseg, 
                    points.darkblue=1-BAFoutputchr$BAFseg, 
                    x.min=min(BAFoutputchr$Position)/1000000, 
                    x.max=max(BAFoutputchr$Position)/1000000, 
                    title=paste(samplename,", chromosome ", chr, sep=""), 
                    xlab="Position (Mb)", 
                    ylab="BAF (phased)",
                    svs_pos=svs_chrom$position/1000000)
    dev.off()
    
    BAFoutputchr$BAFphased = ifelse(BAFoutputchr$tempBAFsegm>0.5, BAFoutputchr$BAF, 1-BAFoutputchr$BAF)
    # Remove the temp BAFsegm values as they are only needed for plotting
    BAFoutput = rbind(BAFoutput, BAFoutputchr[,c(1:5)])
  }
  colnames(BAFoutput) = c("Chromosome","Position","BAF","BAFphased","BAFseg")
  write.table(BAFoutput, outputfile, sep="\t", row.names=F, col.names=T, quote=F)
}