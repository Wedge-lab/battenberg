########################################################################################
# Generic table reader
########################################################################################
#' Generic reading function using the readr R package, tailored for reading in genomic data
#' @param file Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @param row.names Whether the file contains row names (Default: FALSE)
#' @param stringsAsFactor Legacy parameter that is no longer used (Default: FALSE)
#' @param sep Column separator (Default: \\t)
#' @param chrom_col The column number that contains chromosome denominations. This column will automatically be cast as a character. Should be counted including the row.names (Default: 1)
#' @param skip The number of rows to skip before reading (Default: 0)
#' @return A data frame with contents of the file 
#' @export
read_table_generic = function(file, header=T, row.names=F, stringsAsFactor=F, sep="\t", chrom_col=1, skip=0) {
  # stringsAsFactor is not needed here, but kept for legacy purposes
  
  # Read in first line to obtain the header
  d = readr::read_delim(file=file, delim=sep, col_names=header, n_max=1, skip=skip)
  
  # fetch the name of the first column to set its col_type for reading in the whole file
  # this is needed as readr does not understand the chromosome column properly
  col_types = list()
  for (i in chrom_col) {
    first_colname = colnames(d)[i]
    col_types[[first_colname]] = readr::col_character()
  }
  d = readr::read_delim(file=file, delim=sep, col_names=header, col_types=col_types, skip=skip)
  
  # readr never reads row.names, so this needs to be manually corrected
  if (row.names) {
    row.names(d) = d[,1]
    d = d[,-1]
  }
  # Replace spaces with dots as is the standard with the regular read.table
  colnames(d) = gsub(" ", ".", colnames(d))
  return(d)
}

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
      data<-as.data.frame(read_table_generic(filename))
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
  return(as.data.frame(do.call(rbind, lapply(infiles, FUN=function(x) { read_table_generic(x) }))))
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
  return(as.data.frame(do.call(rbind, lapply(infiles, FUN=function(x) { read_table_generic(x) }))))
}

#' Check if a file exists, if it doesn't, exit non-clean
#' @noRd
assert.file.exists = function(filename) {
  if (!file.exists(filename)) {
    warning(paste("Supplied file does not exist: ", filename, sep=""))
    quit(save="no", status=1)
  }
}
