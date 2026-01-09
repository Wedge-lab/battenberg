#' Transform cel files into BAF and LogR
#'
#' This function takes a cel file from a tumour and a matched normal and
#' extracts the BAF and LogR, which is saved into a single file. The \code{gc_correct}
#' function can read that file and transforms it into separate BAF and LogR files that
#' both Battenberg and ASCAT can use.
#' @param normal_cel_file String that points to the cel file containing the matched normal data
#' @param tumour_cel_file String that points to the cel file containing the tumour data
#' @param output_file String where the BAF and LogR should be written
#' @param snp6_reference_info_file String to the SNP6 reference info file that comes with Battenberg SNP6
#' @param apt_probeset_genotype_exe Path to the apt.probeset.genotype executable (Default $PATH)
#' @param apt_probeset_summarize_exe Path to the apt.probeset.summarize executable (Default $PATH)
#' @param norm_geno_clust_exe Path to the normalize_affy_geno_cluster.pl script (Default $PATH)
#' @author sd11
#' @export
cel2baf_logr <- function(
  normal_cel_file,
  tumour_cel_file,
  output_file,
  snp6_reference_info_file,
  apt_probeset_genotype_exe = "apt-probeset-genotype",
  apt_probeset_summarize_exe = "apt-probeset-summarize",
  norm_geno_clust_exe = "normalize_affy_geno_cluster.pl"
) {
  # Unpack pointers to reference files required during this step
  ref_files <- parse_snp6_ref_file(snp6_reference_info_file)
  GW_SNP6 <- ref_files[ref_files$variable == "GW_SNP6", ]$reference_file
  SNP6_BIRDSEED_MODELS <- ref_files[ref_files$variable == "SNP6_BIRDSEED_MODELS", ]$reference_file
  SNP6_SPECIALSNPS <- ref_files[ref_files$variable == "SNP6_SPECIALSNPS", ]$reference_file
  QUANT_NORM_TARGET <- ref_files[ref_files$variable == "QUANT_NORM_TARGET", ]$reference_file
  LOCFILE <- ref_files[ref_files$variable == "LOCFILE", ]$reference_file
  UNM_NORMALS <- ref_files[ref_files$variable == "UNM_NORMALS", ]$reference_file

  # Unpack the normal cel file
  cmd <- paste(apt_probeset_genotype_exe, "-c", GW_SNP6, "-a birdseed", "--read-models-birdseed", SNP6_BIRDSEED_MODELS, "--special-snps", SNP6_SPECIALSNPS, "--cels", normal_cel_file)
  print(cmd)
  exit_code <- system(cmd, wait = TRUE)
  stopifnot(exit_code == 0)
  # Unpack the tumour cel file
  cmd <- paste(apt_probeset_summarize_exe, "--cdf-file", GW_SNP6, "--analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true", "--target-sketch", QUANT_NORM_TARGET, normal_cel_file, tumour_cel_file)
  print(cmd)
  exit_code <- system(cmd, wait = TRUE)
  stopifnot(exit_code == 0)
  # Construct the LogR and BAF and push that to
  cmd <- paste(norm_geno_clust_exe, UNM_NORMALS, "quant-norm.pm-only.med-polish.expr.summary.txt", "-locfile", LOCFILE, "-out", output_file)
  print(cmd)
  exit_code <- system(cmd, wait = TRUE)
  stopifnot(exit_code == 0)
}

#' Correct the LogR estimates for GC content
#'
#' This function performs GC correction of the LogR
#' data. Sometimes a wave pattern is observed there
#' that correlates with GC content. Internally it uses
#' the ASCAT gc correction function.
#' @param samplename Name of the sample to be used to name columns
#' @param infile.logr.baf String that points to the raw combined BAF and LogR file that is the result of \code{cel2baf_logr}
#' @param outfile.tumor.LogR The filename of the file where the tumour LogR will be written
#' @param outfile.tumor.BAF The filename of the file where the tumour BAF will be written
#' @param outfile.normal.LogR The filename of the file where the normal LogR will be written
#' @param outfile.normal.BAF The filename of the file where the normal BAF will be written
#' @param outfile.probeBAF The filename of the file where the probe ids and their BAF will be saved
#' @param snp6_reference_info_file String to the SNP6 reference info file that comes with Battenberg SNP6
#' @param chr_names A vector of chromosome names that are to be used
#' @param birdseed_report_file Name of the birdseed output file. This is a temp output file of one of the internally called functions of which the name cannot be defined. Don't change this parameter. (Default birdseed.report.txt)
#' @author sd11
#' @export
gc_correct <- function(samplename, infile.logr.baf, outfile.tumor.LogR, outfile.tumor.BAF, outfile.normal.LogR, outfile.normal.BAF, outfile.probeBAF, snp6_reference_info_file, chr_names, birdseed_report_file = "birdseed.report.txt", genomebuild = "hg19") {
  # Read in needed reference files
  ref_files <- parse_snp6_ref_file(snp6_reference_info_file)
  SNP_POS_REF <- ref_files[ref_files$variable == "SNP_POS", ]$reference_file
  GC_SNP6 <- ref_files[ref_files$variable == "GC_SNP6", ]$reference_file

  lrrbaf <- utils::read.table(infile.logr.baf, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  SNPpos <- utils::read.table(SNP_POS_REF, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

  Tumor_LogR <- lrrbaf[rownames(SNPpos), 5, drop = FALSE]
  colnames(Tumor_LogR) <- samplename

  Tumor_BAF <- lrrbaf[rownames(SNPpos), 6, drop = FALSE]
  colnames(Tumor_BAF) <- samplename

  Normal_LogR <- lrrbaf[rownames(SNPpos), 3, drop = FALSE]
  colnames(Normal_LogR) <- samplename

  Normal_BAF <- lrrbaf[rownames(SNPpos), 4, drop = FALSE]
  colnames(Normal_BAF) <- samplename

  # replace 2's by NA
  Tumor_BAF[Tumor_BAF == 2] <- NA
  Normal_BAF[Normal_BAF == 2] <- NA

  # Tumor_LogR: correct difference between copy number only probes and other probes
  CNprobes <- substring(rownames(SNPpos), 1, 2) == "CN"

  Tumor_LogR[CNprobes, 1] <- Tumor_LogR[CNprobes, 1] - mean(Tumor_LogR[CNprobes, 1], na.rm = TRUE)
  Tumor_LogR[!CNprobes, 1] <- Tumor_LogR[!CNprobes, 1] - mean(Tumor_LogR[!CNprobes, 1], na.rm = TRUE)

  Normal_LogR[CNprobes, 1] <- Normal_LogR[CNprobes, 1] - mean(Normal_LogR[CNprobes, 1], na.rm = TRUE)
  Normal_LogR[!CNprobes, 1] <- Normal_LogR[!CNprobes, 1] - mean(Normal_LogR[!CNprobes, 1], na.rm = TRUE)

  # limit the number of digits:
  Tumor_LogR <- round(Tumor_LogR, 4)
  Normal_LogR <- round(Normal_LogR, 4)

  data.table::fwrite(cbind(SNPpos, Tumor_BAF), paste(outfile.tumor.BAF, "_noGCcorr.txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
  data.table::fwrite(cbind(SNPpos, Normal_BAF), paste(outfile.normal.BAF, "_noGCcorr.txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)

  # read into ASCAT and make GC corrected input:
  data.table::fwrite(cbind(SNPpos, Tumor_LogR), paste(outfile.tumor.LogR, "_noGCcorr.txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
  data.table::fwrite(cbind(SNPpos, Normal_LogR), paste(outfile.normal.LogR, "_noGCcorr.txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)

  # ======================================= above previous prepareGCcorrect, below runGCcorrect ==============================================

  # TODO: This must be a dapted to not hardcode the chromosome names
  gender <- utils::read.table(birdseed_report_file, sep = "\t", skip = 66, header = TRUE)
  sex <- as.vector(gender[, "computed_gender"])
  sex[sex == "female"] <- "XX"
  sex[sex == "male"] <- "XY"
  sex[sex == "unknown"] <- NA

  ascat.bc <- ASCAT::ascat.loadData(paste(outfile.tumor.LogR, "_noGCcorr.txt", sep = ""), paste(outfile.tumor.BAF, "_noGCcorr.txt", sep = ""), paste(outfile.normal.LogR, "_noGCcorr.txt", sep = ""), paste(outfile.normal.BAF, "_noGCcorr.txt", sep = ""), chrs = chr_names, gender = sex, genomeVersion = genomebuild)
  ASCAT::ascat.plotRawData(ascat.bc)
  ascat.bc <- ASCAT::ascat.correctLogR(ascat.bc, GC_SNP6)

  # Make sure the right column names are added here, because these are expected by fitcopynumber
  colnames(ascat.bc$SNPpos) <- c("Chromosome", "Position")

  # Determine SNPs with BAF between 0.3-0.7 from normal => these are supposed to be heterozygous
  is.het <- (ascat.bc$Germline_BAF >= 0.3 & ascat.bc$Germline_BAF <= 0.7)
  dat <- cbind(ascat.bc$SNPpos, round(ascat.bc$Germline_LogR, 4))
  dat <- dat[which(is.het), ]
  colnames(dat) <- c("Chromosome", "Position", samplename)
  data.table::fwrite(dat, file = outfile.normal.LogR, row.names = FALSE, quote = FALSE, sep = "\t")

  select <- !is.na(ascat.bc$Germline_BAF)
  dat <- cbind(ascat.bc$SNPpos, round(ascat.bc$Germline_BAF, 4))
  colnames(dat) <- c("Chromosome", "Position", samplename)
  data.table::fwrite(dat[which(select), ], file = outfile.normal.BAF, row.names = FALSE, quote = FALSE, sep = "\t")

  # Save the probe ids plus their BAF for only the germline heterozygous mutations
  select <- !is.na(ascat.bc$Tumor_BAF)
  dat <- cbind(row.names(ascat.bc$SNPpos), ascat.bc$Tumor_BAF)
  dat <- dat[which(select & is.het), ]
  data.table::fwrite(dat, file = outfile.probeBAF, row.names = FALSE, quote = FALSE, col_names = FALSE, sep = "\t")

  # Save tumour BAF and LogR directly. Include homozygous SNPs here.
  dat <- cbind(ascat.bc$SNPpos, round(ascat.bc$Tumor_BAF, 4))
  dat <- dat[which(select), ]
  colnames(dat) <- c("Chromosome", "Position", samplename)
  data.table::fwrite(dat, file = outfile.tumor.BAF, row.names = FALSE, quote = FALSE, sep = "\t")

  select <- !is.na(ascat.bc$Tumor_LogR)
  dat <- cbind(ascat.bc$SNPpos, round(ascat.bc$Tumor_LogR, 4))
  dat <- dat[which(select), ]
  colnames(dat) <- c("Chromosome", "Position", samplename)
  data.table::fwrite(dat, file = outfile.tumor.LogR, row.names = FALSE, quote = FALSE, sep = "\t")
}


#' Prepares data for impute
#'
#' The raw BAF and LogR data have been dumped into separate files. Now the data
#' needs to be prepared to go into Impute2, which is essentially morphing it into
#' the correct format. This function does that per chromosome and can therefore
#' be run in parallel for each chromosome.
#' @param infile_germlineBAF Germline BAF file generated by \code{cel2baf_logr}
#' @param infile_tumourBAF Tumour BAF file generated by \code{cel2baf_logr}
#' @param outFileStart Prefix of the filenames where the Impute2 input will be written. These will be extended with the chromosome
#' @param chrom Char with the chromosome for which an Impute2 file is produced
#' @param chr_names A vector of chromosome names that can be considered. This vector can just contain the chromosome for which the Impute2 file is produced, but can contain all chromosomes.
#' @param problem_loci_file A string that points to a file with problematic loci that should be removed from the data
#' @param snp6_reference_info_file String to the SNP6 reference info file that comes with Battenberg SNP6
#' @param imputeinfofile String to the impute 1000 genomes reference info file that comes with Battenberg
#' @param is_male Boolean that is True if the donor is male, False when female
#' @param heterozygous_filter BAF cutoff for calling homozygous SNPs
#' @author dw9 jd
#' @export
generate_impute_input_snp6 <- function(
  infile_germlineBAF,
  infile_tumourBAF,
  outFileStart,
  chrom,
  chr_names,
  problem_loci_file,
  snp6_reference_info_file,
  imputeinfofile,
  is_male,
  heterozygous_filter = "none"
) {
  ref_files <- parse_snp6_ref_file(snp6_reference_info_file)
  ANNO_FILE <- ref_files[ref_files$variable == "ANNO_FILE", "reference_file"]

  impute_info <- parse_imputeinfofile(imputeinfofile, is_male, chrom = chrom)

  known_SNPs <- data.table::rbindlist(
    lapply(impute_info$impute_legend, function(f) {
      data.table::fread(f, header = TRUE)
    })
  )

  allele_levels <- c("A", "C", "G", "T")
  data.table::set(
    known_SNPs,
    j = "position",
    value = as.integer(known_SNPs[["position"]])
  )
  data.table::set(
    known_SNPs,
    j = "a0",
    value = factor(known_SNPs[["a0"]], levels = allele_levels)
  )
  data.table::set(
    known_SNPs,
    j = "a1",
    value = factor(known_SNPs[["a1"]], levels = allele_levels)
  )

  if (!is.na(problem_loci_file) && problem_loci_file != "NA") {
    problemSNPs <- data.table::fread(problem_loci_file, header = TRUE)
    bad_pos <- problemSNPs[
      problemSNPs[["Chr"]] == chrom,
      problemSNPs[["Pos"]]
    ]
    known_SNPs <- known_SNPs[!(known_SNPs[["position"]] %in% bad_pos)]
  }

  knownSNP6data <- data.table::fread(ANNO_FILE, skip = "#", header = TRUE)
  knownSNP6data <- knownSNP6data[knownSNP6data[["Chromosome"]] == chrom]

  complement <- c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
  neg_strand <- knownSNP6data[["Strand"]] == "-"

  data.table::set(
    knownSNP6data,
    i = which(neg_strand),
    j = "Allele.A",
    value = complement[knownSNP6data[["Allele.A"]][neg_strand]]
  )
  data.table::set(
    knownSNP6data,
    i = which(neg_strand),
    j = "Allele.B",
    value = complement[knownSNP6data[["Allele.B"]][neg_strand]]
  )

  knownSNP6data <- knownSNP6data[!duplicated(knownSNP6data[["Physical.Position"]])]

  data.table::set(
    knownSNP6data,
    j = "Allele.A",
    value = factor(knownSNP6data[["Allele.A"]], levels = allele_levels)
  )
  data.table::set(
    knownSNP6data,
    j = "Allele.B",
    value = factor(knownSNP6data[["Allele.B"]], levels = allele_levels)
  )

  germline_snp_data <- data.table::fread(infile_germlineBAF, header = TRUE)
  chr_col <- names(germline_snp_data)[1]
  germline_snp_data <- germline_snp_data[germline_snp_data[[chr_col]] == chrom]

  tumour_snp_data <- data.table::fread(infile_tumourBAF, header = TRUE)
  chr_col <- names(tumour_snp_data)[1]
  tumour_snp_data <- tumour_snp_data[tumour_snp_data[[chr_col]] == chrom]

  data.table::setnames(germline_snp_data, c("Chr", "Pos", "nBAF"))
  data.table::setnames(tumour_snp_data, c("Chr", "Pos", "tBAF"))

  snp_data <- merge(germline_snp_data, tumour_snp_data, by = c("Chr", "Pos"))

  anno_subset <- data.table::data.table(
    Physical.Position = knownSNP6data[["Physical.Position"]],
    Allele.A = knownSNP6data[["Allele.A"]],
    Allele.B = knownSNP6data[["Allele.B"]]
  )

  matched.info <- merge(
    anno_subset,
    snp_data,
    by.x = "Physical.Position",
    by.y = "Pos"
  )

  combined.info <- merge(
    matched.info,
    known_SNPs,
    by.x = "Physical.Position",
    by.y = "position"
  )

  idx_match <- combined.info[["Allele.A"]] == combined.info[["a0"]] &
    combined.info[["Allele.B"]] == combined.info[["a1"]]

  idx_flip <- combined.info[["Allele.A"]] == combined.info[["a1"]] &
    combined.info[["Allele.B"]] == combined.info[["a0"]]

  combined.info1 <- combined.info[idx_match]
  combined.info2 <- combined.info[idx_flip]

  data.table::set(
    combined.info2,
    j = "nBAF",
    value = 1.0 - combined.info2[["nBAF"]]
  )
  data.table::set(
    combined.info2,
    j = "tBAF",
    value = 1.0 - combined.info2[["tBAF"]]
  )

  all.info <- data.table::rbindlist(list(combined.info1, combined.info2))
  all.info <- all.info[order(all.info[["Physical.Position"]])]

  is_het_vec <- all.info[["nBAF"]] >= 0.3 & all.info[["nBAF"]] <= 0.7

  utils::write.csv(
    all.info[is_het_vec, setdiff(names(all.info), "nBAF"), drop = FALSE],
    file = paste0(outFileStart, chrom, "_withAlleleFreq.csv"),
    quote = FALSE,
    row.names = FALSE
  )

  if (heterozygous_filter != "none") {
    minBaf <- min(heterozygous_filter, 1.0 - heterozygous_filter)
    maxBaf <- max(heterozygous_filter, 1.0 - heterozygous_filter)

    is_het <- all.info[["nBAF"]] >= 0.3 & all.info[["nBAF"]] <= 0.7
    is_hom_ref <- all.info[["nBAF"]] <= minBaf
    is_hom_alt <- all.info[["nBAF"]] >= maxBaf

    keep <- is_het | is_hom_ref | is_hom_alt
    subset <- all.info[keep]

    out.data <- data.table::data.table(
      snp.names = paste0("snp", seq_len(nrow(subset))),
      ID = subset[["id"]],
      Pos = subset[["Physical.Position"]],
      a0 = subset[["a0"]],
      a1 = subset[["a1"]],
      G1 = as.integer(is_hom_ref[keep]),
      G2 = as.integer(is_het[keep]),
      G3 = as.integer(is_hom_alt[keep])
    )
  } else {
    subset <- all.info[is_het_vec]

    out.data <- data.table::data.table(
      snp.names = paste0("snp", seq_len(nrow(subset))),
      ID = subset[["id"]],
      Pos = subset[["Physical.Position"]],
      a0 = subset[["a0"]],
      a1 = subset[["a1"]],
      G1 = 0L,
      G2 = 1L,
      G3 = 0L
    )
  }

  data.table::fwrite(
    out.data,
    file = paste0(outFileStart, chrom, ".txt"),
    col.names = FALSE,
    quote = FALSE,
    sep = " "
  )

  if (chrom == "chrX") {
    sample_g_data <- data.frame(
      ID_1 = c(0, "INDIVI1"),
      ID_2 = c(0, "INDIVI1"),
      missing = c(0, 0),
      sex = c("D", 2)
    )
    data.table::fwrite(
      sample_g_data,
      file = paste0(outFileStart, "sample_g.txt"),
      sep = " "
    )
  }
}
#' Prepare SNP6 data for haplotype construction
#'
#' This function performs part of the Battenberg SNP6 pipeline: Extract BAF and logR from the CEL files
#' and performing GC content correction.
#'
#' @param tumour_cel_file Full path to a CEL file containing the tumour raw data
#' @param normal_cel_file Full path to a CEL file containing the normal raw data
#' @param tumourname Identifier to be used for tumour output files
#' @param chrom_names A vector containing the names of chromosomes to be included
#' @param snp6_reference_info_file Full path to the SNP6 reference info file
#' @param apt_probeset_genotype_exe Full path to the apt.probeset.genotype executable (Default: expected in $PATH)
#' @param apt_probeset_summarize_exe Full path to the apt.probeset.summarize executable (Default: expected in $PATH)
#' @param norm_geno_clust_exe  Full path to the norm_geno_clust_exe executable (Default: expected in $PATH)
#' @param birdseed_report_file Name of the birdseed output file. This is a temp output file of one of the internally called functions of which the name cannot be defined. Don't change this parameter. (Default: birdseed.report.txt)
#' @author sd11
#' @export
prepare_snp6 <- function(
  tumour_cel_file, normal_cel_file,
  tumourname, chrom_names,
  snp6_reference_info_file,
  apt_probeset_genotype_exe = "apt-probeset-genotype",
  apt_probeset_summarize_exe = "apt-probeset-summarize",
  norm_geno_clust_exe = "normalize_affy_geno_cluster.pl",
  birdseed_report_file = "birdseed.report.txt",
  genomebuild = "hg19"
) {
  # Extract the LogR and BAF from both tumour and normal cel files.
  cel2baf_logr(
    normal_cel_file = normal_cel_file,
    tumour_cel_file = tumour_cel_file,
    output_file = paste(tumourname, "_lrr_baf.txt", sep = ""),
    snp6_reference_info_file = snp6_reference_info_file,
    apt_probeset_genotype_exe = apt_probeset_genotype_exe,
    apt_probeset_summarize_exe = apt_probeset_summarize_exe,
    norm_geno_clust_exe = norm_geno_clust_exe
  )

  gc_correct(
    samplename = tumourname,
    infile.logr.baf = paste(tumourname, "_lrr_baf.txt", sep = ""),
    outfile.tumor.LogR = paste(tumourname, "_mutantLogR.tab", sep = ""),
    outfile.tumor.BAF = paste(tumourname, "_mutantBAF.tab", sep = ""),
    outfile.normal.LogR = paste(tumourname, "_germlineLogR.tab", sep = ""),
    outfile.normal.BAF = paste(tumourname, "_germlineBAF.tab", sep = ""),
    outfile.probeBAF = paste(tumourname, "_probeBAF.txt", sep = ""),
    snp6_reference_info_file = snp6_reference_info_file,
    birdseed_report_file = birdseed_report_file,
    chr_names = chrom_names,
    genomebuild = genomebuild
  )
}
