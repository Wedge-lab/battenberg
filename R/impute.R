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
      EXIT_CODE=system(cmd, wait=T)
      stopifnot(EXIT_CODE==0)
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
    impute.info = impute.info[impute.info$chrom==chrom,]
  }
  return(impute.info)
}

#' Check impute info file consistency
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @author sd11
check.imputeinfofile = function(imputeinfofile, is.male, usebeagle) {
  impute.info = parse.imputeinfofile(imputeinfofile, is.male)
  if (usebeagle){
    if (any(!file.exists(impute.info$impute_legend))) {
      print("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
      stop("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
    }
  } else {
    if (any(!file.exists(impute.info$impute_legend) | !file.exists(impute.info$genetic_map) | !file.exists(impute.info$impute_hap))) {
      print("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
      stop("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
    }
  }
}

#' Returns the chromosome names that are supported
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @param is.male A boolean describing whether the sample under study is male.
#' @param chrom The name of a chromosome to subset the contents of the imputeinfofile with (optional)
#' @param analaysis Depending on the type of analysis different sets of chromosomes are returned (Default:  paired)
#' @return A vector containing the supported chromosome names
#' @author sd11
#' @export
get.chrom.names = function(imputeinfofile, is.male, chrom=NA, analysis="paired") {
  chrom_names = unique(parse.imputeinfofile(imputeinfofile, is.male, chrom=chrom)$chrom)
  if (analysis=="cell_line" | analysis=="germline") {
    # Both cell line and germline analysis do not yield usable data on X and Y, so remove
    chrom_names = chrom_names[!chrom_names %in% c("X", "Y")]
  }
  return(chrom_names)
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



#' Converts impute input to a beagle input
#'
#' This function takes the impute input file and converts it to a beagle input
#'
#' @param imputeinput path to the impute input file
#' @param chrom chromosome
#' @author maxime.tarabichi
#' @export
convert.impute.input.to.beagle.input = function(imputeinput,
                                                chrom)
{
  chrom <- if(chrom=="23") "X" else chrom
  inp <- read_impute_input(imputeinput)
  coln <- c("#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "SAMP001")
  vcf <- cbind(rep(chrom,nrow(inp)),
               inp[,3],
               rep(".",nrow(inp)),
               inp[,4],
               inp[,5],
               rep(".",nrow(inp)),
               rep("PASS",nrow(inp)),
               rep(".",nrow(inp)),
               rep("GT",nrow(inp)),
               paste(inp$X6,inp$X7,inp$X8,sep="-"), stringsAsFactors = F)
  vcf[vcf[,10]=="1-0-0",10] <- "0/0"
  vcf[vcf[,10]=="0-1-0",10] <- "0/1"
  vcf[vcf[,10]=="0-0-1",10] <- "1/1"
  vcf <- vcf[vcf[,10]!="0-0-0",]
  colnames(vcf) <- coln
  vcf
}

#' Writes input file for beagle5
#'
#' This function writes a table formatted as a vcf to the drive for beagle5 to run on
#'
#' @param vcf data frame vcf-like for beagle
#' @param filepath character string for path to the file to write on disk
#' @param vcfversion character string for version for the vcf (default 4.2)
#' @param genomereference character string for genome build (default GRCh37)
#' @author maxime.tarabichi
#' @export
writevcf.beagle = function(vcf,
                           filepath,
                           vcfversion="4.2",
                           genomereference="GRCh37")
{
  cat(paste0('##fileformat=VCFv',vcfversion,
             '\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##reference=',
             genomereference,
             '\n'),
      file=filepath)
  suppressWarnings(write.table(vcf,
                               file=filepath,
                               sep="\t",col.names=T,row.names=F,quote=F,append=T))
}


#' Writes output of beagle as output from impute (interface bealge/impute for Battenberg)
#'
#' This function writes a table formatted as a vcf to the drive for beagle5 to run on
#'
#' @param vcf character string path for output from beagle
#' @param outfile character string path for impute-like outputfile
#' @author maxime.tarabichi
#' @export
writebeagle.as.impute = function(vcf,
                                 outfile)
{
  beagleout <- read_beagle_output(vcf)
  haplotypes <- strsplit(beagleout$SAMP001,split="\\|")
  dt <- cbind(paste0("snp_index",1:nrow(beagleout)),
              paste0("rs_index",1:nrow(beagleout)),
              beagleout[,2],
              beagleout[,4],
              beagleout[,5],
              sapply(haplotypes,"[",1),
              sapply(haplotypes,"[",2))
  write.table(dt,
              file=outfile,
              quote=F,
              col.names=F,
              row.names=F,
              sep="\t")
}


#' Command to run beagle5
#'
#' This runs beagle through a system call to the beagle java jar file.
#' It requires pre-formatted reference and plink files for the correct genome build.
#'
#' @param beaglejar character string path to Beagle5 java jar file
#' @param vcfpath character string path to the vcf input file to be phased
#' @param reffile character string path to the Beagle5 reference file
#' @param outpath character string path to Beagle's output vcf.gz file
#' @param plinkfile character string path to the plink file
#' @param nthreads integer number of threads
#' @param window integer max size of genomic window to be phased (cM; default 40; decrease for less memory usage; should be >1.1*overlap)
#' @param overlap integer overlap of windows (cM; default 4)
#' @param javajre Path to the Java JRE executable (default java, i.e. in $PATH)
#' @param maxheap.gb integer maximum heap size for the java process in gigabytes (default 10)
#' @author maxime.tarabichi
#' @export
run.beagle5 = function(beaglejar,
                       vcfpath,
                       reffile,
                       outpath,
                       plinkfile,
                       nthreads=1,
                       window=40,
                       overlap=4,
                       maxheap.gb=10,
		       javajre="java")
{
    cmd <- paste0(javajre,
		  " -Xmx",maxheap.gb,"g",
		  " -Xms", maxheap.gb, "g",
		  " -XX:+UseParallelOldGC",
                  " -jar ",beaglejar,
                  " gt=",vcfpath,
                  " ref=",reffile ,
                  " out=",outpath,
                  " map=",plinkfile,
                  " nthreads=",nthreads,
                  " window=",window,
                  " overlap=",overlap,
                  " impute=false")
    EXIT_CODE=system(cmd, wait=T)
    stopifnot(EXIT_CODE==0)
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
#' @param usebeagle Should use beagle5 instead of impute2 Default: FALSE
#' @param beaglejar Full path to Beagle java jar file Default: NA
#' @param beagleref Full path to Beagle reference file Default: NA
#' @param beagleplink Full path to Beagle plink file  Default: NA
#' @param beaglemaxmem Integer Beagle max heap size in Gb  Default: 10
#' @param beaglenthreads Integer number of threads used by beagle5 Default:1
#' @param beaglewindow Integer size of the genomic window for beagle5 (cM) Default:40
#' @param beagleoverlap Integer size of the overlap between windows beagle5 Default:4
#' @param javajre Path to the Java JRE executable (default java, i.e. in $PATH)
#' @author sd11, maxime.tarabichi, jdemeul
#' @export
run_haplotyping = function(chrom, tumourname, normalname, ismale, imputeinfofile, problemloci, impute_exe, min_normal_depth, chrom_names, 
                           externalhaplotypeprefix = NA,
                           use_previous_imputation=F,
                           snp6_reference_info_file=NA, heterozygousFilter=NA,
                           usebeagle=FALSE,
                           beaglejar=NA,
                           beagleref=NA,
                           beagleplink=NA,
                           beaglemaxmem=10,
                           beaglenthreads=1,
                           beaglewindow=40,
                           beagleoverlap=4,
			   javajre="java")
{
  
  previoushaplotypefile <- list.files(pattern = paste0("_impute_output_chr", chrom, "_allHaplotypeInfo.txt"))[1]
  if (use_previous_imputation & !is.na(previoushaplotypefile)) {
    
    print(paste0("Previous imputation results found, copying info from", previoushaplotypefile, " to flip alleles"))
    currenthaplotypefile <- paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep="")
    if (previoushaplotypefile != currenthaplotypefile) {
      file.copy(from = previoushaplotypefile, to = paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep=""))
    }
    
  } else {
    
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
    
    if(usebeagle){
      ## Convert input files for beagle5
      imputeinputfile <- paste(tumourname,
                               "_impute_input_chr",
                               chrom, ".txt", sep="")
      vcfbeagle <- convert.impute.input.to.beagle.input(imputeinput=imputeinputfile,
                                                        chrom=chrom)
      vcfbeagle_path <- paste(tumourname,"_beagle5_input_chr",chrom,".txt",sep="")
      outbeagle_path <- paste(tumourname,"_beagle5_output_chr",chrom,".txt",sep="")
      writevcf.beagle(vcfbeagle, filepath=vcfbeagle_path)
      ## Run beagle5 on the files
      run.beagle5(beaglejar=beaglejar,
                  vcfpath=vcfbeagle_path,
                  reffile=beagleref,
                  outpath=outbeagle_path,
                  plinkfile=beagleplink,
                  maxheap.gb=beaglemaxmem,
                  nthreads=beaglenthreads,
                  window=beaglewindow,
                  overlap=beagleoverlap,
                  javajre=javajre)
      outfile <- paste(tumourname,
                       "_impute_output_chr",
                       chrom, "_allHaplotypeInfo.txt", sep="")
      vcfout <- paste(outbeagle_path,".vcf.gz",sep="")
      ## Convert beagle output file to impute2-like file
      writebeagle.as.impute(vcf=vcfout,
                            outfile=outfile)
    }
    else {
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
      # Cleanup temp Impute output
      unlink(paste(tumourname, "_impute_output_chr", chrom, ".txt*K.txt*", sep=""))
    }
    
  }
  

  # If an allele counts file exists we assume this is a WGS sample and run the corresponding step, otherwise it must be SNP6
  allelefrequenciesfile <- paste0(tumourname, "_alleleFrequencies_chr", chrom, ".txt")
  print(allelefrequenciesfile)
  print(file.exists(allelefrequenciesfile))
  
  if (file.exists(allelefrequenciesfile)) {
    # WGS - Transform the impute output into haplotyped BAFs
    
    # if present, input external haplotype blocks
    if (!is.na(externalhaplotypeprefix) && file.exists(paste0(externalhaplotypeprefix, chrom, ".vcf"))) {
      print("Adding in the external haplotype blocks")
      
      # output BAFs to plot pre-external haplotyping
      GetChromosomeBAFs(chrom=chrom,
                        SNP_file=allelefrequenciesfile,
                        haplotypeFile=paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep=""),
                        samplename=tumourname,
                        outfile=paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt", sep=""),
                        chr_names=chrom_names,
                        minCounts=min_normal_depth)
      
      # Plot what we have before external haplotyping is incorporated
      plot.haplotype.data(haplotyped.baf.file=paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt", sep=""),
                          imageFileName=paste(tumourname,"_chr",chrom,"_heterozygousData_noExt.png",sep=""),
                          samplename=tumourname,
                          chrom=chrom,
                          chr_names=chrom_names)
      
      input_known_haplotypes(chrom = chrom,
                             chrom_names = chrom_names,
                             imputedHaplotypeFile = paste0(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"),
                             externalHaplotypeFile = paste0(externalhaplotypeprefix, chrom, ".vcf"))
      
    }
    
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
}
