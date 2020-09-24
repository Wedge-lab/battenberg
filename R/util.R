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
  d = readr::read_delim(file=file, delim=sep, col_names=header, n_max=1, skip=skip, col_types = readr::cols())
  
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

#' Parser for logR data
#' @param filename Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @return A data frame with logR content
read_logr = function(filename, header=T) {
  return(readr::read_tsv(file = filename, col_names = header, col_types = "cin"))
}

#' Parser for BAF data
#' @param filename Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @return A data frame with BAF content
read_baf = function(filename, header=T) {
  return(readr::read_tsv(file = filename, col_names = header, col_types = "cin"))
}

#' Parser for GC content reference data
#' @param filename Filename of the file to read in
#' @return A data frame with GC content
read_gccontent = function(filename) {
  return(readr::read_tsv(file=filename, skip = 1, col_names = F, col_types = "-cinnnnnnnnnnnn------"))
}

#' Parser for replication timing reference data
#' @param filename Filename of the file to read in
#' @return A data frame with replication timing
read_replication = function(filename) {
  return(readr::read_tsv(file=filename, col_types = paste0("ci", paste0(rep("n", 15), collapse = ""))))
}

#' Parser for BAFsegmented data
#' @param filename Filename of the file to read in
#' @param header Whether the file contains a header (Default: TRUE)
#' @return A data frame with BAFsegmented content
read_bafsegmented = function(filename, header=T) {
  return(readr::read_tsv(file = filename, col_names = header, col_types = "cinnn"))
}

#' Parser for imputed genotype data
#' @param filename Filename of the file to read in
#' @return A data frame with the imputed genotype output
read_imputed_output = function(filename) {
  return(readr::read_tsv(file = filename, col_names = c("snpidx", "rsidx", "pos", "ref", "alt", "hap1", "hap2"), col_types = "cciccii"))
}

#' Parser for allele frequencies data
#' @param filename Filename of the file to read in
#' @return A data frame with the alleleCounter output
read_alleleFrequencies = function(filename) {
  return(readr::read_tsv(file = filename, col_names = c("CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth"), col_types = "ciiiiii", comment = "#"))
}

#' Parser for impute input data
#' @param filename Filename of the file to read in
#' @return A data frame with the input for impute
read_impute_input = function(filename) {
  return(readr::read_delim(file = filename, col_names = F, col_types = "ccicciii", delim = " "))
}

#' Parser for beagle5 output data
#' @param filename Filename of the file to read in
#' @return A data frame with the beagle5 output
read_beagle_output = function(filename) {
  return(readr::read_tsv(file = filename, col_names = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMP001"), col_types = "cicccccccc", comment = "#"))
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
concatenateBAFfiles<-function(inputStart, inputEnd, outputFile, chr_names) {
  all_data<-NULL
  colNames<-NULL
  for(i in chr_names)
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
concatenateAlleleCountFiles = function(inputStart, inputEnd, chr_names) {
  infiles = c()
  for(chrom in chr_names) {
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
concatenateG1000SnpFiles = function(inputStart, inputEnd, chr_names) {
  data = list()
  for(chrom in chr_names) {
    filename = paste(inputStart, chrom, inputEnd, sep="")
    # Only add files that exist and have data
    if(file.exists(filename) && file.info(filename)$size>0) {
      # infiles = c(infiles, filename)
      data[[chrom]] = cbind(chromosome=chrom, read_table_generic(filename))
    }
  }
  return(as.data.frame(do.call(rbind, data)))
}



########################################################################################
# Various functions for calculating from data
########################################################################################
#' Calc copy number of major allele per segment from a subclones data.frame
#' @noRd
calc_total_cn_major = function(bb) {
  return(bb$nMaj1_A*bb$frac1_A + ifelse(bb$frac1_A < 1, bb$nMaj2_A*bb$frac2_A, 0))
}

#' Calc copy number of minor allele per segment from a subclones data.frame
#' @noRd
calc_total_cn_minor = function(bb) {
  return(bb$nMin1_A*bb$frac1_A + ifelse(bb$frac1_A < 1, bb$nMin2_A*bb$frac2_A, 0))
}

#' Calc total copy number per segment from a subclones data.frame
#' @noRd
calculate_bb_total_cn = function(bb) {
  return((bb$nMaj1_A+bb$nMin1_A)*bb$frac1_A + ifelse(!is.na(bb$frac2_A), (bb$nMaj2_A+bb$nMin2_A)*bb$frac2_A, 0))
}

#' Calc ploidy from a subclones data.frame
#' @noRd
calc_ploidy = function(bb) {
  bb$len = bb$endpos/1000-bb$startpos/1000
  bb$total_cn = calculate_bb_total_cn(bb)
  ploidy = sum(bb$total_cn*bb$len) / sum(bb$len)
  return(ploidy)
}

#' Transform logR into an estimate of total copy number given purity and total ploidy (tumour+normal)
#' @noRd
logr2tumcn = function(cellularity, total_ploidy, logR) {
  return(((total_ploidy*(2^logR)) - 2*(1-cellularity)) / cellularity)
}

#' Calc psi from psi_t and rho
#' @noRd
psit2psi = function(rho, psi_t) {
  return(rho*psi_t + 2*(1-rho))
}

#' Calc psi_t from psi and rho
#' @noRd
psi2psit = function(rho, psi) {
  return((psi-2*(1-rho))/rho)
}

########################################################################################
# Refitting functions
########################################################################################
#' Calculate rho and psi values from a refit suggestion
#' 
#' Use this function to calculate the refit values from a refit suggestion.
#' @param refBAF BAF of the segment
#' @param refLogR logR of the segment
#' @param refMajor Major allele copy number
#' @param refMinor Minor allele copy number
#' @param rho Sample rho parameter
#' @param gamma_param Platform gamma parameter
#' @return A list with a field for rho and psi_t
#' @author sd11
#' @export
calc_rho_psi_refit = function(refBAF, refLogR, refMajor, refMinor, rho, gamma_param) {
  rho = (2*refBAF-1)/(2*refBAF-refBAF*(refMajor+refMinor)-1+refMajor)
  psi = (rho*(refMajor+refMinor)+2-2*rho)/(2^(refLogR/gamma_param))
  psi_t = psi2psit(rho, psi)
  return(list(rho=rho, psi_t=psi_t))
}

#' Calculate refit values from a refit suggestion
#' 
#' Use this function to calculate the refit values from a refit suggestion.
#' @param subclones_file A Battenberg subclones.txt file
#' @param segment_chrom Chromsome of the segment to use for refitting
#' @param segment_pos Position within the start/end coordinates of the segment to use for refitting
#' @param new_nMaj Major allele copy number
#' @param new_nMin Minor allele copy number
#' @param rho Sample rho parameter
#' @param gamma_param Platform gamma parameter
#' @return A list with a field for rho and psi_t
#' @author sd11
#' @export
suggest_refit = function(subclones_file, segment_chrom, segment_pos, new_nMaj, new_nMin, rho, gamma_param) {
  # segment_pos = as.numeric(gsub("M", "000000", segment_pos))
  subclones = read.table(subclones_file, header=T, stringsAsFactors=F)
  segment = subclones[subclones$chr==segment_chrom & subclones$startpos<=segment_pos & subclones$endpos>=segment_pos,]
  segment_BAF = segment$BAF
  segment_LogR = segment$LogR
  return(calc_rho_psi_refit(segment_BAF, segment_LogR, new_nMaj, new_nMin, rho, gamma_param))
}

#' Create refit suggestions for a fit copy number profile
#' 
#' This function takes a fit copy number profile and generates refit suggestions for a future rerun.
#' If there are clonal alterations above a specified size, then those written out as supplied as suggestions,
#' otherwise a refit suggestion of an external purity value will be saved.
#' @param samplename Samplename for the output file
#' @param subclones_file File containing a fit copy number profile
#' @param rho_psi_file File with rho and psi values
#' @param gamma_param Platform gamma parameter
#' @param min_segment_size_mb Minimum size of a segment in Mb to be considered for a refit suggestion (Default: 2)
#' @author sd11
#' @export
cnfit_to_refit_suggestions = function(samplename, subclones_file, rho_psi_file, gamma_param, min_segment_size_mb=2) {
  # samplename = "NASCR-0016"
  # subclones_file = "NASCR-0016_subclones.txt"
  subclones = Battenberg::read_table_generic(subclones_file)
  subclones$len = subclones$endpos/1000000-subclones$startpos/1000000
  subclones$is_cna = subclones$nMaj1_A!=subclones$nMin1_A
  
  if (any(subclones$len > min_segment_size_mb & subclones$is_cna)) {
    # There are large scale alterations, save the top couple as suggestions
    rho_psi = read.table(rho_psi_file, header=T, stringsAsFactors=F)
    rho = rho_psi["FRAC_GENOME", "rho"]
    psi_t = rho_psi["FRAC_GENOME", "psi"]
    
    # Take only segments that are clonal and are an alteration
    is_subclonal = subclones$frac1_A < 1
    subclones_clonal_cna = subset(subclones, !is_subclonal & subclones$is_cna)
    subclones_clonal_cna = subclones_clonal_cna[with(subclones_clonal_cna, order(len, decreasing=T)),]

    if (nrow(subclones_clonal_cna)==0) {
	output = data.frame(project=NA, samplename=samplename, qc=NA, cellularity_refit=T, chrom=NA, pos=NA, maj=NA, min=NA, baf=NA, logr=NA, rho_estimate=NA, psi_t_estimate=NA, rho_diff=NA, psi_t_diff=NA)
    } else {
    
    # Generate a couple of solutions, but not more than are possibly available
    max_solutions = ifelse(nrow(subclones_clonal_cna) >= 5, 5, nrow(subclones_clonal_cna))
    subclones_clonal_cna = subclones_clonal_cna[1:max_solutions, , drop=F]
    
    # Determine position in Mb within the segment
    position = subclones_clonal_cna$startpos + (subclones_clonal_cna$endpos - subclones_clonal_cna$startpos) / 2
    position = position / 1000000
    position_round_up = ceiling(position)
    position_round_down = floor(position)
    position = ifelse(position_round_up < subclones_clonal_cna$endpos, position_round_up, position_round_down)
    
    output = data.frame(project=rep(NA, max_solutions),
                        samplename=rep(samplename, max_solutions),
                        qc=rep(NA, max_solutions),
                        cellularity_refit=rep(F, max_solutions),
                        chrom=subclones_clonal_cna$chr[1:max_solutions],
                        pos=paste(position, "M", sep=""),
                        maj=subclones_clonal_cna$nMaj1_A[1:max_solutions],
                        min=subclones_clonal_cna$nMin1_A[1:max_solutions],
                        baf=subclones_clonal_cna$BAF[1:max_solutions],
                        logr=subclones_clonal_cna$LogR[1:max_solutions])
    
    #refBAF, refLogR, refMajor, refMinor, rho, gamma_param
    res = calc_rho_psi_refit(output$baf, output$logr, output$maj, output$min, rho, gamma_param)
    output$rho_estimate = res$rho
    output$psi_t_estimate = res$psi_t
    output$rho_diff = abs(rho-output$rho_estimate)
    output$psi_t_diff = abs(psi_t-output$psi_t_estimate)
    }
  } else {
    # No large clonal alteration, save a suggestion that should use an external purity value
    output = data.frame(project=NA, samplename=samplename, qc=NA, cellularity_refit=T, chrom=NA, pos=NA, maj=NA, min=NA, baf=NA, logr=NA, rho_estimate=NA, psi_t_estimate=NA, rho_diff=NA, psi_t_diff=NA)
  }
  write.table(output, file=paste0(samplename, "_refit_suggestion.txt"), quote=F, sep="\t", row.names=F)
}

########################################################################################
# Other
########################################################################################
#' Check if a file exists, if it doesn't, exit non-clean
#' @noRd
assert.file.exists = function(filename) {
  if (!file.exists(filename)) {
    warning(paste("Supplied file does not exist: ", filename, sep=""))
    quit(save="no", status=1)
  }
}
