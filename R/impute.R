#' Run impute on the specified inputfile
#'
#' This function runs impute across the input using the specified region.size.
#' @param inputfile Full path to a csv file with columns: Physical.Position, Allele.A, Allele.B, allele.frequency, id ,position, a0, a1
#' @param outputfile.prefix Prefix to the output file. Region boundaries are added as suffix.
#' @param is.male Boolean describing whether the sample is male (TRUE) or female (FALSE)
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @param impute.exe Pointer to where the impute2 executable can be found (optional).
#' @param region.size An integer describing the region size to be used by impute (optional).
#' @param chrom The name of a chromosome on which this function should run (names are used, supply X as 'X') (optional).
#' @param seed The seed to be set
#' @author dw9
#' @export
run.impute = function(inputfile, outputfile.prefix, is.male, imputeinfofile, impute.exe="impute2", region.size=5000000, chrom=NA, seed=as.integer(Sys.time())) {
  
  # Read in the impute file information
  impute.info = parse.imputeinfofile(imputeinfofile, is.male, chrom=chrom)
  
  # Run impute for each region of the size specified above
  for(r in 1:nrow(impute.info)){
    boundaries = seq(as.numeric(impute.info[r,]$start),as.numeric(impute.info[r,]$end),region.size)
    if(boundaries[length(boundaries)] != impute.info[r,]$end){
      boundaries = c(boundaries,impute.info[r,]$end)
    }

    # Take the start of the region+1 here to make sure there are no overlapping regions, wich causes a 
    # problem with SNPs on exactly the boundary. It does mean the first base on the first chromosome
    # cannot be phased
    for(b in 1:(length(boundaries)-1)){
      cmd = paste(impute.exe,
                  " -m ", impute.info[r,]$genetic_map,
                  " -h ", impute.info[r,]$impute_hap,
                  " -l ", impute.info[r,]$impute_legend,
                  " -g ", inputfile,
                  " -int ", boundaries[b]+1, " ", boundaries[b+1],
                  " -Ne 20000", # Authors of impute2 mention that this parameter works best on all population types, thus hardcoded.
                  " -o ", outputfile.prefix, "_", boundaries[b]/1000, "K_", boundaries[b+1]/1000, "K.txt", 
                  " -phase", 
                  " -seed ",
                  " -os 2", sep="") # lowers computational cost by not imputing reference only SNPs
      system(cmd, wait=T)
    }
  }
}

#' Read in the imputeinfofile.
#'
#' Reads in a file with the following columns: 
#'   chromosome : 1-X
#'   impute_legend : Legend file in IMPUTE -l format
#'   genetic_map : Genetic map file in IMPUTE -m format
#'   impute_hap : Phased haplotype file in IMPUTE -h format
#'   start : Start of the chromosome 
#'   end : End of the chromosome
#'   is_par : 1 when pseudo autosomal region, 0 when not
#'   
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @param is.male A boolean describing whether the sample under study is male.
#' @param chrom The name of a chromosome to subset the contents of the imputeinfofile with (optional)
#' @return A data.frame with 7 columns: Chromosome, impute_legend, genetic_map, impute_hap, start, end, is_par
#' @author sd11
#' @export
parse.imputeinfofile = function(imputeinfofile, is.male, chrom=NA) {
  impute.info = read.table(imputeinfofile, stringsAsFactors=F)
  colnames(impute.info) = c("chrom", "impute_legend", "genetic_map", "impute_hap", "start", "end", "is_par")
  # Remove the non-pseudo autosomal region (i.e. where not both men and woman are diploid)
  if(is.male){ impute.info = impute.info[impute.info$is_par==1,] }
  chr_names=unique(impute.info$chrom)
  # Subset for a particular chromosome
  if (!is.na(chrom)) {
    impute.info = impute.info[impute.info$chrom==chr_names[chrom],]
  }
  return(impute.info)
}

#' Check impute info file consistency
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @author sd11
check.imputeinfofile = function(imputeinfofile, is.male) {
  impute.info = parse.imputeinfofile(imputeinfofile, is.male)
  if (any(!file.exists(impute.info$impute_legend) | !file.exists(impute.info$genetic_map) | !file.exists(impute.info$impute_hap))) {
    print("Could not find reference files, make sure paths in impute_info.txt point to the correct location")  
    stop("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
  }
}

#' Returns the chromosome names that are supported
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @param is.male A boolean describing whether the sample under study is male.
#' @param chrom The name of a chromosome to subset the contents of the imputeinfofile with (optional)
#' @return A vector containing the supported chromosome names
#' @author sd11
#' @export
get.chrom.names = function(imputeinfofile, is.male, chrom=NA) {
  return(unique(parse.imputeinfofile(imputeinfofile, is.male, chrom=chrom)$chrom))
}

#' Concatenate the impute output generated for each of the regions.
#'
#' This function assembles the impute output generated.
#' @param inputfile.prefix Prefix of the input files (this is typically the outputfile.prefix option supplied when calling run.impute).
#' @param outputfile Where to store the output.
#' @param is.male Boolean describing whether the sample is male (TRUE) or female (FALSE).
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @param region.size An integer describing the region size to be used by impute (optional).
#' @param chrom The name of a chromosome on which this function should run (names are used, supply X as 'X').
#' @author dw9
#' @export
combine.impute.output = function(inputfile.prefix, outputfile, is.male, imputeinfofile, region.size=5000000, chrom=NA) {
  # Read in the impute file information
  impute.info = parse.imputeinfofile(imputeinfofile, is.male, chrom=chrom)
  
  # Assemble the start and end points of all regions
  all.boundaries = array(0,c(0,2))
  for(r in 1:nrow(impute.info)){
    boundaries = seq(as.numeric(impute.info[r,]$start),as.numeric(impute.info[r,]$end),region.size)
    if(boundaries[length(boundaries)] != impute.info[r,]$end){
      boundaries = c(boundaries,impute.info[r,]$end)
    }
    all.boundaries = rbind(all.boundaries,cbind(boundaries[-(length(boundaries))],boundaries[-1]))
  }
  # Concatenate all the regions  
  impute.output = concatenateImputeFiles(inputfile.prefix, all.boundaries)
  write.table(impute.output, file=outputfile, row.names=F, col.names=F, quote=F, sep=" ")
}


#' Construct haplotypes for a chromosome
#' 
#' This function takes preprocessed data and performs haplotype reconstruction.
#' 
#' @param chrom The chromosome for which to reconstruct haplotypes
#' @param tumourname Identifier of the tumour, used to match data files on disk
#' @param normalname Identifier of the normal, used to match data files on disk
#' @param ismale Boolean, set to TRUE if the sample is male
#' @param imputeinfofile Full path to the imputeinfo reference file
#' @param problemloci Full path to the problematic loci reference file
#' @param impute_exe Path to the impute executable (can be found if its in $PATH)
#' @param min_normal_depth Minimal depth in the matched normal required for a SNP to be used
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param snp6_reference_info_file SNP6 only parameter Default: NA
#' @param heterozygousFilter SNP6 only parameter Default: NA
#' @author sd11
#' @export
run_haplotyping = function(chrom, tumourname, normalname, ismale, imputeinfofile, problemloci, impute_exe, min_normal_depth, chrom_names,
                           snp6_reference_info_file=NA, heterozygousFilter=NA) {
  
  if (file.exists(paste(tumourname, "_alleleFrequencies_chr", chrom, ".txt", sep=""))) {
    generate.impute.input.wgs(chrom=chrom,
                              tumour.allele.counts.file=paste(tumourname,"_alleleFrequencies_chr", chrom, ".txt", sep=""),
                              normal.allele.counts.file=paste(normalname,"_alleleFrequencies_chr", chrom, ".txt", sep=""),
                              output.file=paste(tumourname, "_impute_input_chr", chrom, ".txt", sep=""),
                              imputeinfofile=imputeinfofile,
                              is.male=ismale,
                              problemLociFile=problemloci,
                              useLociFile=NA)
  } else {
    generate.impute.input.snp6(infile.germlineBAF=paste(tumourname, "_germlineBAF.tab", sep=""),
                               infile.tumourBAF=paste(tumourname, "_mutantBAF.tab", sep=""),
                               outFileStart=paste(tumourname, "_impute_input_chr", sep=""),
                               chrom=chrom,
                               chr_names=chrom_names,
                               problemLociFile=problemloci,
                               snp6_reference_info_file=snp6_reference_info_file,
                               imputeinfofile=imputeinfofile,
                               is.male=ismale,
                               heterozygousFilter=heterozygousFilter)
  }

  # Run impute on the files
  run.impute(inputfile=paste(tumourname, "_impute_input_chr", chrom, ".txt", sep=""),
             outputfile.prefix=paste(tumourname, "_impute_output_chr", chrom, ".txt", sep=""),
             is.male=ismale,
             imputeinfofile=imputeinfofile,
             impute.exe=impute_exe,
             region.size=5000000,
             chrom=chrom)

  # As impute runs in windows across a chromosome we need to assemble the output
  combine.impute.output(inputfile.prefix=paste(tumourname, "_impute_output_chr", chrom, ".txt", sep=""),
                        outputfile=paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep=""),
                        is.male=ismale,
                        imputeinfofile=imputeinfofile,
                        region.size=5000000,
                        chrom=chrom)

  # If an allele counts file exists we assume this is a WGS sample and run the corresponding step, otherwise it must be SNP6
  print(paste(tumourname, "_alleleFrequencies_chr", chrom, ".txt", sep=""))
  print(file.exists(paste(tumourname, "_alleleFrequencies_chr", chrom, ".txt", sep="")))
  if (file.exists(paste(tumourname, "_alleleFrequencies_chr", chrom, ".txt", sep=""))) {
    # WGS - Transform the impute output into haplotyped BAFs
    GetChromosomeBAFs(chrom=chrom,
                      SNP_file=paste(tumourname, "_alleleFrequencies_chr", chrom, ".txt", sep=""),
                      haplotypeFile=paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep=""),
                      samplename=tumourname,
                      outfile=paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                      chr_names=chrom_names,
                      minCounts=min_normal_depth)
  } else {
    print("SNP6 get BAFs")
    # SNP6 - Transform the impute output into haplotyped BAFs
    GetChromosomeBAFs_SNP6(chrom=chrom,
                           alleleFreqFile=paste(tumourname, "_impute_input_chr", chrom, "_withAlleleFreq.csv", sep=""),
                           haplotypeFile=paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep=""),
                           samplename=tumourname,
                           outputfile=paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                           chr_names=chrom_names)
  }

  # Plot what we have until this point
  plot.haplotype.data(haplotyped.baf.file=paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                      imageFileName=paste(tumourname,"_chr",chrom,"_heterozygousData.png",sep=""),
                      samplename=tumourname,
                      chrom=chrom,
                      chr_names=chrom_names)

  # Cleanup temp Impute output
  unlink(paste(tumourname, "_impute_output_chr", chrom, ".txt*K.txt*", sep=""))
}
