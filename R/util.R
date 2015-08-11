########################################################################################
# Concatenate files
########################################################################################
#' Function to concatenate Impute output
#' @noRd
concatenateImputeFiles<-function(inputStart, boundaries) { #outputFile, 
  infiles = c()
  for(i in 1:nrow(boundaries)) {
    filename = paste(inputStart,"_",boundaries[i,1]/1000,"K_",boundaries[i,2]/1000,"K.txt_haps",sep="")
    # Only add files that exist and have data
    if(file.exists(filename) && file.info(filename)$size>0) {
      infiles = c(infiles, filename)
    }
  }
  return(do.call(rbind, lapply(infiles, FUN=function(x) { read.table(x, sep=" ") })))
}

#' Function to concatenate haplotyped BAF output
#' @noRd
concatenateBAFfiles<-function(inputStart, inputEnd, outputFile, no.chrs) {
  all_data<-NULL
  colNames<-NULL
  for(i in 1:no.chrs)
  {
    filename = paste(inputStart,i,inputEnd,sep="")
    if(file.exists(filename) && file.info(filename)$size>0)
    {
      data<-read.table(filename, sep = "\t", header=T)
      all_data<-rbind(all_data,data)
      colNames<-names(data)
    }
  }
  #rnames=paste("snp",1:nrow(all_data),sep="")
  write.table(all_data,outputFile, row.names=F, col.names=colNames, quote=F, sep="\t")
}

#' Function to concatenate allele counter output
#' @noRd
concatenateAlleleCountFiles = function(inputStart, inputEnd, no.chrs) {
  infiles = c()
  for(chrom in 1:no.chrs) {
    filename = paste(inputStart, chrom, inputEnd, sep="")
    # Only add files that exist and have data
    if(file.exists(filename) && file.info(filename)$size>0) {
      infiles = c(infiles, filename)
    }
  }
  return(do.call(rbind, lapply(infiles, FUN=function(x) { read.table(x, header=T, sep="\t", comment.char="") })))
}

#' Function to concatenate 1000 Genomes SNP reference files
#' @noRd
concatenateG1000SnpFiles = function(inputStart, inputEnd, no.chrs) {
  infiles = c()
  for(chrom in 1:no.chrs) {
    filename = paste(inputStart, chrom, inputEnd, sep="")
    # Only add files that exist and have data
    if(file.exists(filename) && file.info(filename)$size>0) {
      infiles = c(infiles, filename)
    }
  }
  #infiles = paste(inputStart, 1:no.chrs, inputEnd, sep="")
  return(do.call(rbind, lapply(infiles, FUN=function(x) { read.table(x, sep="\t", header=T) })))
}

#' Check if a file exists, if it doesn't, exit non-clean
#' @noRd
assert.file.exists = function(filename) {
  if (!file.exists(filename)) {
    warning(paste("Supplied file does not exist: ", filename, sep=""))
    quit(save="no", status=1)
  }
}