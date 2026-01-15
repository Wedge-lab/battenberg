#' Run impute on the specified inputfile
#'
#' This function runs impute across the input using the specified region.size.
#' @param inputfile Full path to a csv file with columns: Physical.Position, Allele.A, Allele.B, allele.frequency, id ,position, a0, a1
#' @param outputfile_prefix Prefix to the output file. Region boundaries are added as suffix.
#' @param is_male Boolean describing whether the sample is male (TRUE) or female (FALSE)
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @param impute.exe Pointer to where the impute2 executable can be found (optional).
#' @param region.size An integer describing the region size to be used by impute (optional).
#' @param chrom The name of a chromosome on which this function should run (names are used, supply X as 'X') (optional).
#' @param seed The seed to be set
#' @author dw9
#' @export
run_impute <- function(
  inputfile, outputfile_prefix, is_male,
  imputeinfofile, impute.exe = "impute2",
  region.size = 5000000, chrom = NA,
  seed = as.integer(Sys.time())
) {
  # Read in the impute file information
  impute_info <- parse_imputeinfofile(imputeinfofile, is_male, chrom = chrom)

  # Run impute for each region of the size specified above
  for (r in seq_len(nrow(impute_info))) {
    boundaries <- seq(as.numeric(impute_info[r, ]$start), as.numeric(impute_info[r, ]$end), region.size)
    if (boundaries[length(boundaries)] != impute_info[r, ]$end) {
      boundaries <- c(boundaries, impute_info[r, ]$end)
    }

    # Take the start of the region+1 here to make sure there are no overlapping regions, wich causes a
    # problem with SNPs on exactly the boundary. It does mean the first base on the first chromosome
    # cannot be phased
    for (b in 1:(length(boundaries) - 1)) {
      cmd <- paste(impute.exe,
        " -m ", impute_info[r, ]$genetic_map,
        " -h ", impute_info[r, ]$impute_hap,
        " -l ", impute_info[r, ]$impute_legend,
        " -g ", inputfile,
        " -int ", boundaries[b] + 1, " ", boundaries[b + 1],
        " -Ne 20000", # Authors of impute2 mention that this parameter works best on all population types, thus hardcoded.
        " -o ", outputfile_prefix, "_", boundaries[b] / 1000, "K_", boundaries[b + 1] / 1000, "K.txt",
        " -phase",
        " -seed ",
        " -os 2",
        sep = ""
      ) # lowers computational cost by not imputing reference only SNPs
      exit_code <- system(cmd, wait = TRUE)
      stopifnot(exit_code == 0)
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
#' @param is_male A boolean describing whether the sample under study is male.
#' @param chrom The name of a chromosome to subset the contents of the imputeinfofile with (optional)
#' @return A data.frame with 7 columns: Chromosome, impute_legend, genetic_map, impute_hap, start, end, is_par
#' @author sd11
#' @export
parse_imputeinfofile <- function(imputeinfofile, is_male, chrom = NA) {
  is_par <- NULL
  # Use fread for high-speed reading.
  impute_info <- data.table::fread(
    imputeinfofile,
    col.names = c(
      "chrom", "impute_legend", "genetic_map",
      "impute_hap", "start", "end", "is_par"
    ),
    stringsAsFactors = FALSE
  )
  # Efficient filtering using data.table's internal optimization
  if (is_male) {
    impute_info <- impute_info[is_par == 1]
  }
  # Subset for a particular chromosome
  if (!is.na(chrom)) {
    target_chrom <- chrom
    impute_info <- impute_info[chrom == target_chrom]
  }
  return(impute_info)
}

#' Check impute info file consistency
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @author sd11
check_imputeinfofile <- function(imputeinfofile, is_male, usebeagle) {
  impute_info <- parse_imputeinfofile(imputeinfofile, is_male)
  if (usebeagle) {
    if (any(!file.exists(impute_info$impute_legend))) {
      log_failure("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
    }
  } else {
    if (any(!file.exists(impute_info$impute_legend) | !file.exists(impute_info$genetic_map) | !file.exists(impute_info$impute_hap))) {
      log_failure()("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
    }
  }
}

#' Returns the chromosome names that are supported
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @param is_male A boolean describing whether the sample under study is male.
#' @param chrom The name of a chromosome to subset the contents of the imputeinfofile with (optional)
#' @param analaysis Depending on the type of analysis different sets of chromosomes are returned (Default:  paired)
#' @return A vector containing the supported chromosome names
#' @author sd11
#' @export
get_chrom_names <- function(imputeinfofile, is_male, chrom = NA, analysis = "paired") {
  chrom_names <- unique(parse_imputeinfofile(imputeinfofile, is_male, chrom = chrom)$chrom)
  if (analysis == "cell_line" || analysis == "germline") {
    # Both cell line and germline analysis do not yield usable data on X and Y, so remove
    chrom_names <- chrom_names[!chrom_names %in% c("X", "Y")]
  }
  return(chrom_names)
}

#' Concatenate the impute output generated for each of the regions.
#'
#' This function assembles the impute output generated.
#' @param inputfile.prefix Prefix of the input files (this is typically the outputfile_prefix option supplied when calling run_impute).
#' @param outputfile Where to store the output.
#' @param is_male Boolean describing whether the sample is male (TRUE) or female (FALSE).
#' @param imputeinfofile Path to the imputeinfofile on disk.
#' @param region.size An integer describing the region size to be used by impute (optional).
#' @param chrom The name of a chromosome on which this function should run (names are used, supply X as 'X').
#' @author dw9
#' @export
combine_impute_output <- function(inputfile.prefix, outputfile, is_male, imputeinfofile, region.size = 5000000, chrom = NA) {
  # Read in the impute file information
  impute_info <- parse_imputeinfofile(imputeinfofile, is_male, chrom = chrom)

  # Assemble the start and end points of all regions
  all.boundaries <- array(0, c(0, 2))
  for (r in seq_len(nrow(impute_info))) {
    boundaries <- seq(as.numeric(impute_info[r, ]$start), as.numeric(impute_info[r, ]$end), region.size)
    if (boundaries[length(boundaries)] != impute_info[r, ]$end) {
      boundaries <- c(boundaries, impute_info[r, ]$end)
    }
    all.boundaries <- rbind(all.boundaries, cbind(boundaries[-(length(boundaries))], boundaries[-1]))
  }
  # Concatenate all the regions
  impute.output <- concatenateImputeFiles(inputfile.prefix, all.boundaries)
  data.table::fwrite(
    impute.output,
    file = outputfile,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = " "
  )
}


#' Converts impute input to a beagle input
#'
#' This function takes the impute input file and converts it to a beagle input
#'
#' @param imputeinput path to the impute input file
#' @param chrom chromosome
#' @author maxime.tarabichi
#' @export
convert_impute_input_to_beagle_input <- function(imputeinput, chrom) {
  # :: syntax and pure comments
  chrom_str <- ifelse(chrom == "23" | chrom == "chr23", "X", as.character(chrom))

  # inp is a data.frame from read_impute_input
  inp <- read_impute_input(imputeinput)
  coln <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMP001")

  # Pure comments instead of numbering
  # Handle the case where the impute input is empty
  if (nrow(inp) == 0) {
    empty_vcf <- matrix(character(), nrow = 0, ncol = 10)
    colnames(empty_vcf) <- coln
    return(empty_vcf)
  }

  clean_pos <- as.integer(as.numeric(trimws(inp[, 3])))

  # Build the genotype string using the X6, X7, X8 naming we forced
  # If the columns don't exist, paste will return "NA-NA-NA" which we handle
  gt_raw <- paste(inp$X6, inp$X7, inp$X8, sep = "-")

  # Use a data.frame to prevent vector collapsing
  vcf_df <- data.frame(
    CHROM = rep(chrom_str, nrow(inp)),
    POS = clean_pos,
    ID = rep(".", nrow(inp)),
    REF = inp$X4,
    ALT = inp$X5,
    QUAL = rep(".", nrow(inp)),
    FILTER = rep("PASS", nrow(inp)),
    INFO = rep(".", nrow(inp)),
    FORMAT = rep("GT", nrow(inp)),
    GT = gt_raw,
    stringsAsFactors = FALSE
  )

  # Standardize genotypes
  vcf_df$GT[vcf_df$GT == "1-0-0"] <- "0/0"
  vcf_df$GT[vcf_df$GT == "0-1-0"] <- "0/1"
  vcf_df$GT[vcf_df$GT == "0-0-1"] <- "1/1"

  vcf_df <- vcf_df[vcf_df$GT %in% c("0/0", "0/1", "1/1"), ]
  colnames(vcf_df) <- coln

  return(vcf_df)
}
#' Writes input file for beagle5
#'
#' @param vcf data frame vcf-like for beagle
#' @param filepath character string for path (e.g., "data.vcf")
#' @param vcfversion character string (default 4.2)
#' @param genomereference character string (default GRCh37)
#' @importFrom data.table fwrite
#' @export
writevcf_beagle <- function(vcf,
                            filepath,
                            vcfversion = "4.2",
                            genomereference = "GRCh37") {
  # :: syntax used
  # Pure comments instead of numbering

  vcf_df <- base::as.data.frame(vcf, stringsAsFactors = FALSE)

  header <- base::paste0(
    "##fileformat=VCFv", vcfversion, "\n",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    "##reference=", genomereference, "\n"
  )

  # Write header first
  base::cat(header, file = filepath)

  # Safely handle the #CHROM column name requirement
  actual_names <- base::colnames(vcf_df)
  if (base::length(actual_names) > 0) {
    actual_names[1] <- base::paste0("#", base::gsub("^#", "", actual_names[1]))
    base::colnames(vcf_df) <- actual_names
  }
  data.table::fwrite(
    x = vcf_df,
    file = filepath,
    sep = "\t",
    append = TRUE,
    col.names = TRUE,
    quote = FALSE,
    nThread = 2,
  )
}

#' Writes output of beagle as output from impute (interface bealge/impute for Battenberg)
#'
#' This function writes a table formatted as a vcf to the drive for beagle5 to run on
#'
#' @param vcf character string path for output from beagle
#' @param outfile character string path for impute-like outputfile
#' @author maxime.tarabichi
#' @export
writebeagle_as_impute <- function(vcf,
                                  outfile) {
  beagleout <- read_beagle_output(vcf)
  haplotypes <- strsplit(beagleout$SAMP001, split = "\\|")
  dt <- cbind(
    paste0("snp_index", seq_len(nrow(beagleout))),
    paste0("rs_index", seq_len(nrow(beagleout))),
    beagleout[, 2],
    beagleout[, 4],
    beagleout[, 5],
    sapply(haplotypes, "[", 1),
    sapply(haplotypes, "[", 2)
  )
  data.table::fwrite(dt,
    file = outfile,
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
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
#' @param maxheap_gb integer maximum heap size for the java process in gigabytes (default 10)
#' @author maxime.tarabichi
#' @export
run_beagle5 <- function(beaglejar,
                        vcfpath,
                        reffile,
                        outpath,
                        plinkfile,
                        nthreads = 1,
                        window = 40,
                        overlap = 4,
                        maxheap_gb = 10,
                        javajre = "java") {
  cmd <- paste0(
    javajre,
    " -Xmx", maxheap_gb, "g",
    " -Xms", maxheap_gb, "g",
    " -XX:+UseParallelGC",
    " -jar ", beaglejar,
    " gt=", vcfpath,
    " ref=", reffile,
    " out=", outpath,
    " map=", plinkfile,
    " nthreads=", nthreads,
    " window=", window,
    " overlap=", overlap,
    " impute=false"
  )
  exit_code <- system(cmd, wait = TRUE)
  stopifnot(exit_code == 0)
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
#' @param heterozygous_filter SNP6 only parameter Default: NA
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
run_haplotyping <- function(
  chrom, tumourname, normalname,
  ismale, imputeinfofile, problemloci,
  impute_exe, min_normal_depth, chrom_names,
  externalhaplotypeprefix = NA,
  use_previous_imputation = FALSE,
  snp6_reference_info_file = NA,
  heterozygous_filter = NA,
  usebeagle = FALSE,
  beaglejar = NA,
  beagleref = NA,
  beagleplink = NA,
  beaglemaxmem = 10,
  beaglenthreads = 1,
  beaglewindow = 40,
  beagleoverlap = 4,
  javajre = "java"
) {
  previoushaplotypefile <- list.files(pattern = paste0("_impute_output_chr", chrom, "_allHaplotypeInfo.txt"))[1]
  if (use_previous_imputation && !is.na(previoushaplotypefile)) {
    print(paste0("Previous imputation results found, copying info from", previoushaplotypefile, " to flip alleles"))
    currenthaplotypefile <- paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = "")
    if (previoushaplotypefile != currenthaplotypefile) {
      file.copy(from = previoushaplotypefile, to = paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = ""))
    }
  } else {
    if (file.exists(paste(tumourname, "_alleleFrequencies_chr", chrom, ".txt", sep = ""))) {
      log_info("Generating WGS impute input for chr{chrom}")
      generate_impute_input_wgs(
        chrom = chrom,
        tumour_allele_counts_file = paste(tumourname, "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
        normal_allele_counts_file = paste(normalname, "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
        output_file = paste(tumourname, "_impute_input_chr", chrom, ".txt", sep = ""),
        imputeinfofile = imputeinfofile,
        is_male = ismale,
        problem_loci_file = problemloci,
        use_loci_file = NA
      )
    } else {
      log_info("Generating SNP6 impute input for chr{chrom}")
      generate_impute_input_snp6(
        infile_germlineBAF = paste(tumourname, "_germlineBAF.tab", sep = ""),
        infile_tumourBAF = paste(tumourname, "_mutantBAF.tab", sep = ""),
        outFileStart = paste(tumourname, "_impute_input_chr", sep = ""),
        chrom = chrom,
        chr_names = chrom_names,
        problem_loci_file = problemloci,
        snp6_reference_info_file = snp6_reference_info_file,
        imputeinfofile = imputeinfofile,
        is_male = ismale,
        heterozygous_filter = heterozygous_filter
      )
    }

    if (usebeagle) {
      log_info("Mode: Beagle5 for chr{chrom}")
      ## Convert input files for beagle5
      imputeinputfile <- paste(tumourname,
        "_impute_input_chr",
        chrom, ".txt",
        sep = ""
      )
      vcfbeagle <- convert_impute_input_to_beagle_input(
        imputeinput = imputeinputfile,
        chrom = chrom
      )
      log_info("successfully converted impute input to beagle input")
      vcfbeagle_path <- paste(tumourname, "_beagle5_input_chr", chrom, ".txt", sep = "")
      outbeagle_path <- paste(tumourname, "_beagle5_output_chr", chrom, ".txt", sep = "")
      writevcf_beagle(vcfbeagle, filepath = vcfbeagle_path)
      ## Run beagle5 on the files
      log_info("Calling run_beagle5 for chr{chrom}")
      run_beagle5(
        beaglejar = beaglejar,
        vcfpath = vcfbeagle_path,
        reffile = beagleref,
        outpath = outbeagle_path,
        plinkfile = beagleplink,
        maxheap_gb = beaglemaxmem,
        nthreads = beaglenthreads,
        window = beaglewindow,
        overlap = beagleoverlap,
        javajre = javajre
      )
      outfile <- paste(tumourname,
        "_impute_output_chr",
        chrom, "_allHaplotypeInfo.txt",
        sep = ""
      )
      vcfout <- paste(outbeagle_path, ".vcf.gz", sep = "")
      log_info("Converting Beagle VCF back to Impute format for chr{chrom}")
      writebeagle_as_impute(
        vcf = vcfout,
        outfile = outfile
      )
    } else {
      # Run impute on the files
      run_impute(
        inputfile = paste(tumourname, "_impute_input_chr", chrom, ".txt", sep = ""),
        outputfile_prefix = paste(tumourname, "_impute_output_chr", chrom, ".txt", sep = ""),
        is_male = ismale,
        imputeinfofile = imputeinfofile,
        impute.exe = impute_exe,
        region.size = 5000000,
        chrom = chrom
      )

      # As impute runs in windows across a chromosome we need to assemble the output
      combine_impute_output(
        inputfile.prefix = paste(tumourname, "_impute_output_chr", chrom, ".txt", sep = ""),
        outputfile = paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = ""),
        is_male = ismale,
        imputeinfofile = imputeinfofile,
        region.size = 5000000,
        chrom = chrom
      )
      # Cleanup temp Impute output
      unlink(paste(tumourname, "_impute_output_chr", chrom, ".txt*K.txt*", sep = ""))
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
      log_info("Adding in the external haplotype blocks")

      # output BAFs to plot pre-external haplotyping
      GetChromosomeBAFs(
        chrom = chrom,
        SNP_file = allelefrequenciesfile,
        haplotypeFile = paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = ""),
        samplename = tumourname,
        outfile = paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt", sep = ""),
        chr_names = chrom_names,
        minCounts = min_normal_depth
      )

      # Plot what we have before external haplotyping is incorporated
      plot_haplotype_data(
        haplotyped_baf_file = paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt", sep = ""),
        image_file_name = paste(tumourname, "_chr", chrom, "_heterozygousData_noExt.png", sep = ""),
        samplename = tumourname,
        chrom = chrom
      )

      input_known_haplotypes(
        chrom = chrom,
        chrom_names = chrom_names,
        imputedHaplotypeFile = paste0(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"),
        externalHaplotypeFile = paste0(externalhaplotypeprefix, chrom, ".vcf")
      )
    }

    GetChromosomeBAFs(
      chrom = chrom,
      SNP_file = paste(tumourname, "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
      haplotypeFile = paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = ""),
      samplename = tumourname,
      outfile = paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
      chr_names = chrom_names,
      minCounts = min_normal_depth
    )
  } else {
    log_info("SNP6 get BAFs")
    # SNP6 - Transform the impute output into haplotyped BAFs
    GetChromosomeBAFs_SNP6(
      chrom = chrom,
      alleleFreqFile = paste(tumourname, "_impute_input_chr", chrom, "_withAlleleFreq.csv", sep = ""),
      haplotypeFile = paste(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = ""),
      samplename = tumourname,
      outputfile = paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
      chr_names = chrom_names
    )
  }

  # Plot what we have until this point
  plot_haplotype_data(
    haplotyped_baf_file = paste(tumourname, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
    image_file_name = paste(tumourname, "_chr", chrom, "_heterozygousData.png", sep = ""),
    samplename = tumourname,
    chrom = chrom
  )
}

#' Construct haplotypes for a chromosome - germline WGS version
#'
#' This function takes preprocessed data and performs haplotype reconstruction.
#'
#' @param chrom The chromosome for which to reconstruct haplotypes
#' @param germlinename Identifier of the germline sample, used to match data files on disk
#' @param normalname Identifier of the reconstructed normal, used to match data files on disk
#' @param ismale Boolean, set to TRUE if the sample is male
#' @param imputeinfofile Full path to the imputeinfo reference file
#' @param problemloci Full path to the problematic loci reference file
#' @param impute_exe Path to the impute executable (can be found if its in $PATH)
#' @param min_normal_depth Minimal depth in the matched normal required for a SNP to be used
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param snp6_reference_info_file SNP6 only parameter Default: NA
#' @param heterozygous_filter SNP6 only parameter Default: NA
#' @param usebeagle Should use beagle5 instead of impute2 Default: FALSE
#' @param beaglejar Full path to Beagle java jar file Default: NA
#' @param beagleref Full path to Beagle reference file Default: NA
#' @param beagleplink Full path to Beagle plink file  Default: NA
#' @param beaglemaxmem Integer Beagle max heap size in Gb  Default: 10
#' @param beaglenthreads Integer number of threads used by beagle5 Default:1
#' @param beaglewindow Integer size of the genomic window for beagle5 (cM) Default:40
#' @param beagleoverlap Integer size of the overlap between windows beagle5 Default:4
#' @param javajre Path to the Java JRE executable (default java, i.e. in $PATH)
#' @author sd11, maxime.tarabichi, jdemeul, Naser Ansari-Pour (BDI, Oxford)
#' @export

run_haplotyping_germline <- function(chrom, germlinename, normalname, ismale, imputeinfofile, problemloci, impute_exe, min_normal_depth, chrom_names,
                                     externalhaplotypeprefix = NA,
                                     use_previous_imputation = FALSE,
                                     snp6_reference_info_file = NA, heterozygous_filter = NA,
                                     usebeagle = FALSE,
                                     beaglejar = NA,
                                     beagleref = NA,
                                     beagleplink = NA,
                                     beaglemaxmem = 10,
                                     beaglenthreads = 1,
                                     beaglewindow = 40,
                                     beagleoverlap = 4,
                                     javajre = "java") {
  previoushaplotypefile <- list.files(pattern = paste0("_impute_output_chr", chrom, "_allHaplotypeInfo.txt"))[1]
  if (use_previous_imputation && !is.na(previoushaplotypefile)) {
    print(paste0("Previous imputation results found, copying info from", previoushaplotypefile, " to flip alleles"))
    currenthaplotypefile <- paste(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = "")
    if (previoushaplotypefile != currenthaplotypefile) {
      file.copy(from = previoushaplotypefile, to = paste(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = ""))
    }
  } else {
    if (file.exists(paste(germlinename, "_alleleFrequencies_chr", chrom, ".txt", sep = ""))) {
      generate_impute_input_wgs_germline(
        chrom = chrom,
        germline_allele_counts_file = paste(germlinename, "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
        normal_allele_counts_file = paste(normalname, "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
        output_file = paste(germlinename, "_impute_input_chr", chrom, ".txt", sep = ""),
        imputeinfofile = imputeinfofile,
        is_male = ismale,
        problem_loci_file = problemloci,
        use_loci_file = NA
      )
    } else {
      stop("Germline calling is currently on WGS data only - SNP array data is not sufficiently dense to detect all germline CNVs")
    }

    if (usebeagle) {
      ## Convert input files for beagle5
      imputeinputfile <- paste(germlinename,
        "_impute_input_chr",
        chrom, ".txt",
        sep = ""
      )
      vcfbeagle <- convert_impute_input_to_beagle_input(
        imputeinput = imputeinputfile,
        chrom = chrom
      )
      vcfbeagle_path <- paste(germlinename, "_beagle5_input_chr", chrom, ".txt", sep = "")
      outbeagle_path <- paste(germlinename, "_beagle5_output_chr", chrom, ".txt", sep = "")
      writevcf_beagle(vcfbeagle, filepath = vcfbeagle_path)
      ## Run beagle5 on the files
      run_beagle5(
        beaglejar = beaglejar,
        vcfpath = vcfbeagle_path,
        reffile = beagleref,
        outpath = outbeagle_path,
        plinkfile = beagleplink,
        maxheap_gb = beaglemaxmem,
        nthreads = beaglenthreads,
        window = beaglewindow,
        overlap = beagleoverlap,
        javajre = javajre
      )
      outfile <- paste(germlinename,
        "_impute_output_chr",
        chrom, "_allHaplotypeInfo.txt",
        sep = ""
      )
      vcfout <- paste(outbeagle_path, ".vcf.gz", sep = "")
      ## Convert beagle output file to impute2-like file
      writebeagle_as_impute(
        vcf = vcfout,
        outfile = outfile
      )
    } else {
      # Run impute on the files
      run_impute(
        inputfile = paste(germlinename, "_impute_input_chr", chrom, ".txt", sep = ""),
        outputfile_prefix = paste(germlinename, "_impute_output_chr", chrom, ".txt", sep = ""),
        is_male = ismale,
        imputeinfofile = imputeinfofile,
        impute.exe = impute_exe,
        region.size = 5000000,
        chrom = chrom
      )

      # As impute runs in windows across a chromosome we need to assemble the output
      combine_impute_output(
        inputfile.prefix = paste(germlinename, "_impute_output_chr", chrom, ".txt", sep = ""),
        outputfile = paste(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = ""),
        is_male = ismale,
        imputeinfofile = imputeinfofile,
        region.size = 5000000,
        chrom = chrom
      )
      # Cleanup temp Impute output
      unlink(paste(germlinename, "_impute_output_chr", chrom, ".txt*K.txt*", sep = ""))
    }
  }


  # If an allele counts file exists we assume this is a WGS sample and run the corresponding step, otherwise it must be SNP6
  allelefrequenciesfile <- paste0(germlinename, "_alleleFrequencies_chr", chrom, ".txt")
  print(allelefrequenciesfile)
  print(file.exists(allelefrequenciesfile))

  if (file.exists(allelefrequenciesfile)) {
    # WGS - Transform the impute output into haplotyped BAFs

    # if present, input external haplotype blocks
    if (!is.na(externalhaplotypeprefix) && file.exists(paste0(externalhaplotypeprefix, chrom, ".vcf"))) {
      log_info("Adding in the external haplotype blocks")

      # output BAFs to plot pre-external haplotyping
      GetChromosomeBAFs(
        chrom = chrom,
        SNP_file = allelefrequenciesfile,
        haplotypeFile = paste(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = ""),
        samplename = germlinename,
        outfile = paste(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt", sep = ""),
        chr_names = chrom_names,
        minCounts = min_normal_depth
      )

      # Plot what we have before external haplotyping is incorporated
      plot_haplotype_data(
        haplotyped_baf_file = paste(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped_noExt.txt", sep = ""),
        image_file_name = paste(germlinename, "_chr", chrom, "_heterozygousData_noExt.png", sep = ""),
        samplename = germlinename,
        chrom = chrom
      )

      input_known_haplotypes(
        chrom = chrom,
        chrom_names = chrom_names,
        imputedHaplotypeFile = paste0(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"),
        externalHaplotypeFile = paste0(externalhaplotypeprefix, chrom, ".vcf")
      )
    }

    GetChromosomeBAFs(
      chrom = chrom,
      SNP_file = paste(germlinename, "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
      haplotypeFile = paste(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep = ""),
      samplename = germlinename,
      outfile = paste(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
      chr_names = chrom_names,
      minCounts = min_normal_depth
    )
  } else {
    log_failure("Germline calling is only on WGS data - SNParray data not sufficiently dense")
  }

  # Plot what we have until this point
  plot_haplotype_data(
    haplotyped_baf_file = paste(germlinename, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
    image_file_name = paste(germlinename, "_chr", chrom, "_heterozygousData.png", sep = ""),
    samplename = germlinename,
    chrom = chrom
  )
}
