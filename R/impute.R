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
  # Use fread for high-speed reading.
  impute_info <- data.table::fread(
    imputeinfofile,
    col.names = c(
      "chrom", "impute_legend", "genetic_map",
      "impute_hap", "start", "end", "is_par"
    ),
    stringsAsFactors = FALSE
  )

  expected_cols <- c("chrom", "impute_legend", "genetic_map", "impute_hap", "start", "end", "is_par")
  if (!all(expected_cols %in% names(impute_info))) {
    # If columns are missing, try to assign them if possible, or fail
    if (ncol(impute_info) == length(expected_cols)) {
      names(impute_info) <- expected_cols
    } else {
      log_failure("Impute info file does not have the expected number of columns (7). Found: {ncol(impute_info)}")
    }
  }

  # Filter based on gender
  if (!is.na(is_male) && !is_male) {
    # If female, we exclude Y chromosome regions
    # and we might want to handle PAR specifically if the pipeline requires it.
    # But generally, we just want to ensure we don't return Y.
    impute_info <- impute_info[impute_info[["chrom"]] != "Y", ]
  }
  # Subset for a particular chromosome
  if (!is.na(chrom)) {
    impute_info <- impute_info[impute_info[["chrom"]] == chrom, ]
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
      log_failure("Could not find reference files, make sure paths in impute_info.txt point to the correct location")
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
#' @author sd11, maxime.tarabichi, jdemeul
#' @export
convert_beagle_to_impute <- function(beagle_file, output_file) {
  # Read VCF (skip metadata lines starting with ##)
  # We assume VCF has a header line starting with #CHROM
  vcf <- data.table::fread(beagle_file, skip = "#CHROM", header = TRUE)

  # Check if we have enough columns (standard VCF: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE...)
  if (ncol(vcf) < 10) {
    log_failure("Beagle VCF file does not have enough columns: {beagle_file}")
  }

  # Extract GT (Genotype)
  # We assume the last column is the sample genotype (or 10th column)
  # If multisample, this simple converter might need adjustment, but Battenberg usually runs per-sample or tumor/normal
  # For Battenberg pipeline, we typically process one sample's haplotypes here.
  # Let's assume the sample of interest is the first sample column (column 10).
  # If the user provides a multisample VCF, they might need to split it or we pick the first.
  # Given the context of filenames (tumourname_...), it's likely single sample.

  gt_col <- names(vcf)[10]
  gt_data <- vcf[[gt_col]]

  # Split GT string "0|1" -> "0" "1"
  # Beagle output is phased, so pipe | separator
  # We use tstrsplit for efficiency
  haplo <- data.table::tstrsplit(gt_data, "[|/]")

  if (length(haplo) != 2) {
    log_failure("Could not parse genotypes from Beagle VCF. Expected '0|1' format.")
  }

  # Construct IMPUTE2 format
  # 1: "---" (SNP ID placeholder)
  # 2: ID (rsID from VCF) - validation: IMPUTE format often expects non-empty
  # 3: POS
  # 4: REF
  # 5: ALT
  # 6: Hap1
  # 7: Hap2

  impute_dt <- data.table::data.table(
    V1 = "---",
    V2 = vcf$`ID`,
    V3 = vcf$`POS`,
    V4 = vcf$`REF`,
    V5 = vcf$`ALT`,
    V6 = haplo[[1]],
    V7 = haplo[[2]]
  )

  # Write out space-separated, no header (as expected by GetChromosomeBAFs read logic 'header=FALSE')
  data.table::fwrite(impute_dt, file = output_file, sep = " ", col.names = FALSE, quote = FALSE)
}

#' @param impute_results_dir Directory containing the impute/beagle output files
#' @author sd11, maxime.tarabichi, jdemeul
#' @export
run_haplotyping <- function(
  chrom, tumourname, normalname,
  ismale, imputeinfofile, problemloci,
  impute_results_dir, min_normal_depth, chrom_names,
  externalhaplotypeprefix = NA,
  use_previous_imputation = FALSE,
  snp6_reference_info_file = NA,
  heterozygous_filter = NA,
  usebeagle = FALSE
) {
  # Point to the existing haplotype file in the external directory
  if (usebeagle) {
    # Expected Beagle VCF file name
    # We try patterns: .vcf.gz, .vcf
    # Try multiple common naming patterns for Beagle VCFs
    beagle_patterns <- c(
      paste0(tumourname, "_beagle5_output_chr", chrom, ".txt.vcf.gz"),
      paste0(tumourname, "_beagle5_output_chr", chrom, ".txt.vcf"),
      paste0(tumourname, "_beagle_output_chr", chrom, ".vcf.gz"),
      paste0(tumourname, "_beagle_output_chr", chrom, ".vcf")
    )

    beagle_vcf <- NA
    for (pat in beagle_patterns) {
      temp_path <- file.path(impute_results_dir, pat)
      if (file.exists(temp_path)) {
        beagle_vcf <- temp_path
        break
      }
    }

    if (is.na(beagle_vcf)) {
      log_failure("Expected Beagle VCF file not found in {impute_results_dir}. Tried patterns: {paste(beagle_patterns, collapse=', ')}")
    }

    # We need to convert this to IMPUTE format for Battenberg to use
    # We'll create a temporary file or a converted file in the same dir?
    # Ideally in the same dir but we might not have write perms?
    # Let's accept that we write to the same dir or tempdir.
    # To avoid permission issues if impute_results_dir is read-only, we write to tempdir() or current work dir.
    # Current work dir is safer for persistence/debugging.
    haplotype_file <- paste0(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt")
    log_info("Converting Beagle VCF to IMPUTE format: {beagle_vcf} -> {haplotype_file}")
    convert_beagle_to_impute(beagle_vcf, haplotype_file)
  } else {
    haplotype_file <- file.path(impute_results_dir, paste0(tumourname, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"))
    if (!file.exists(haplotype_file)) {
      log_failure("Expected haplotype file not found: {haplotype_file}")
    }
    log_info("Using existing haplotype file: {haplotype_file}")
  }


  # If an allele counts file exists we assume this is a WGS sample and run the corresponding step, otherwise it must be SNP6
  allelefrequenciesfile <- paste0(tumourname, "_alleleFrequencies_chr", chrom, ".txt")

  if (file.exists(allelefrequenciesfile)) {
    # WGS - Transform the impute output into haplotyped BAFs

    # if present, input external haplotype blocks
    if (!is.na(externalhaplotypeprefix) && file.exists(paste0(externalhaplotypeprefix, chrom, ".vcf"))) {
      log_info("Adding in the external haplotype blocks")

      # output BAFs to plot pre-external haplotyping
      GetChromosomeBAFs(
        chrom = chrom,
        SNP_file = allelefrequenciesfile,
        haplotypeFile = haplotype_file,
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
        imputedHaplotypeFile = haplotype_file,
        externalHaplotypeFile = paste0(externalhaplotypeprefix, chrom, ".vcf")
      )
    }

    GetChromosomeBAFs(
      chrom = chrom,
      SNP_file = paste(tumourname, "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
      haplotypeFile = haplotype_file,
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
      haplotypeFile = haplotype_file,
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

#' @param usebeagle Logical, if TRUE expects Beagle VCF output and converts to IMPUTE format.
#' @author sd11, maxime.tarabichi, jdemeul, Naser Ansari-Pour (BDI, Oxford)
#' @export
run_haplotyping_germline <- function(
  chrom, germlinename, normalname, ismale, imputeinfofile, problemloci,
  impute_results_dir, min_normal_depth, chrom_names,
  externalhaplotypeprefix = NA,
  use_previous_imputation = FALSE,
  snp6_reference_info_file = NA, heterozygous_filter = NA,
  usebeagle = FALSE
) {
  # Point to the existing haplotype file in the external directory
  if (usebeagle) {
    # Try multiple common naming patterns for Beagle VCFs
    beagle_patterns <- c(
      paste0(germlinename, "_beagle5_output_chr", chrom, ".txt.vcf.gz"),
      paste0(germlinename, "_beagle5_output_chr", chrom, ".txt.vcf"),
      paste0(germlinename, "_beagle_output_chr", chrom, ".vcf.gz"),
      paste0(germlinename, "_beagle_output_chr", chrom, ".vcf")
    )

    beagle_vcf <- NA
    for (pat in beagle_patterns) {
      temp_path <- file.path(impute_results_dir, pat)
      if (file.exists(temp_path)) {
        beagle_vcf <- temp_path
        break
      }
    }

    if (is.na(beagle_vcf)) {
      log_failure("Expected Beagle VCF file not found in {impute_results_dir}. Tried patterns: {paste(beagle_patterns, collapse=', ')}")
    }

    haplotype_file <- paste0(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt")
    log_info("Converting Beagle VCF to IMPUTE format: {beagle_vcf} -> {haplotype_file}")
    convert_beagle_to_impute(beagle_vcf, haplotype_file)
  } else {
    haplotype_file <- file.path(impute_results_dir, paste0(germlinename, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"))
    if (!file.exists(haplotype_file)) {
      log_failure("Expected haplotype file not found: {haplotype_file}")
    }
    log_info("Using existing haplotype file: {haplotype_file}")
  }

  allelefrequenciesfile <- paste0(germlinename, "_alleleFrequencies_chr", chrom, ".txt")

  if (file.exists(allelefrequenciesfile)) {
    # WGS - Transform the impute output into haplotyped BAFs

    # if present, input external haplotype blocks
    if (!is.na(externalhaplotypeprefix) && file.exists(paste0(externalhaplotypeprefix, chrom, ".vcf"))) {
      log_info("Adding in the external haplotype blocks")

      # output BAFs to plot pre-external haplotyping
      GetChromosomeBAFs(
        chrom = chrom,
        SNP_file = allelefrequenciesfile,
        haplotypeFile = haplotype_file,
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
        imputedHaplotypeFile = haplotype_file,
        externalHaplotypeFile = paste0(externalhaplotypeprefix, chrom, ".vcf")
      )
    }

    GetChromosomeBAFs(
      chrom = chrom,
      SNP_file = paste(germlinename, "_alleleFrequencies_chr", chrom, ".txt", sep = ""),
      haplotypeFile = haplotype_file,
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
