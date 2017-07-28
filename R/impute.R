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
