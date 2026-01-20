#' Run the Battenberg pipeline
#' @param analysis The mode of Battenberg copy number analysis to be undertaken:
#' 'paired' for tumour-normal pair, 'cell_line' for Cell line tumour-only and
#' 'germline' for germline CNV of normal sample (Default: 'paired')
#' @param samplename Sample identifier (tumour or germline), this is used as a
#' prefix for the output files. If allele counts are supplied separately, they
#' are expected to have this identifier as prefix.
#' @param normalname Matched normal identifier, this is used as a prefix for the
#' output files. If allele counts are supplied separately, they are expected to
#' have this identifier as prefix.
#' @param sample_data_file A BAM or CEL file for the sample
#' @param normal_data_file A BAM or CEL file for the
#' normal-pair (paired analysis)
#' @param imputeinfofile Full path to a Battenberg impute info file with
#' pointers to Impute2 reference data
#' @param g1000prefix Full prefix path to 1000 Genomes SNP loci data, as part of
#' the Battenberg reference data
#' @param problemloci Full path to a problem loci file that contains SNP
#' loci that should be filtered out
#' @param gccorrectprefix Full prefix path to GC content files, as part of the
#' Battenberg reference data, not required for SNP6 data (Default: NULL)
#' @param repliccorrectprefix Full prefix path to replication timing files,
#' as part of the Battenberg reference data, not required
#' for SNP6 data (Default: NULL)
#' @param g1000allelesprefix Full prefix path to 1000 Genomes SNP alleles data,
#' as part of the Battenberg reference data, not required for SNP6 data
#' (Default: NA)
#' @param ismale A boolean set to TRUE if the donor is male, set to FALSE if
#' female, not required for SNP6 data (Default: NA)
#' @param data_type String that contains either wgs or snp6 depending on the
#' supplied input data (Default: wgs)
#' @param allele_counts_dir Directory containing the allele counts files (Required for WGS/CellLine/Germline).
#' @param impute_results_dir Directory containing the imputed haplotype results (Required for phasing).
#' @param nthreads The number of concurrent processes to use while running the
#' Battenberg pipeline (Default: 8)
#' @param platform_gamma Platform scaling factor,
#' suggestions are set to 1 for wgs and to 0.55 for snp6 (Default: 1)
#' @param phasing_gamma Gamma parameter used when correcting phasing mistakes
#' (Default: 1)
#' @param segmentation_gamma The gamma parameter
#' controls the size of the penalty
#' of starting a new segment during segmentation.
#' It is therefore the key parameter
#' for controlling the number of segments (Default: 10)
#' @param segmentation_gamma_multisample The gamma parameter
#' controls the size of the penalty of starting a new segment
#' during mutlisample segmentation. It is the
#' key parameter for controlling the number of segments (Default: 10)
#' @param segmentation_kmin Kmin represents the minimum number of
#' probes/SNPs that a segment should consist of (Default: 3)
#' @param phasing_kmin Kmin used when correcting for phasing mistakes
#' (Default: 3)
#' @param clonality_dist_metric Distance metric to use when
#' choosing purity/ploidy combinations (Default: 0)
#' @param ascat_dist_metric Distance metric to use when choosing purity/ploidy
#' combinations (Default: 1)
#' @param min_ploidy Minimum ploidy to be considered (Default: 1.6)
#' @param max_ploidy Maximum ploidy to be considered (Default: 4.8)
#' @param min_rho Minimum purity to be considered (Default: 0.1)
#' @param max_rho Maximum purity to be considered (Default: 1.0)
#' @param min_goodness Minimum goodness of fit required for a purity/ploidy
#' combination to be accepted as a solution (Default: 0.63)
#' @param uninformative_baf_threshold The threshold beyond which BAF becomes
#' uninformative (Default: 0.51)
#' @param min_normal_depth Minimum depth required in the matched normal
#' for a SNP to be considered as part of the wgs analysis (Default: 10)
#' @param min_base_qual Minimum base quality required for a read to
#' be counted when allele counting (Default: 20)
#' @param min_map_qual Minimum mapping quality required for a read to
#' be counted when allele counting (Default: 35)
#' @param max_allowed_state The maximum CN state allowed (Default 250)
#' @param cn_upper_limit Maximum number of copy number that can be called
#' (Default 1000)
#' @param calc_seg_baf_option Sets way to calculate BAF per segment: 1=mean,
#' 2=median, 3=ifelse median==0 | 1, mean, median (Default (paired): 3,
#' cell_line & germline: 1)
#' @param externalhaplotypefile Vcf containing externally
#' obtained haplotype blocks (Default: NA)
#' @param write_battenberg_phasing Write the Battenberg phasing results
#' as vcf to disk, e.g. for multisample cases (Default: TRUE)
#' @param multisample_maxlag Maximal number of upstream SNPs used in the
#' multisample haplotyping to inform the haplotype at another SNP (Default: 100)
#' @param multisample_relative_weight_balanced Relative weight to give to
#' haplotype info from a sample without allelic imbalance
#' in the region (Default: 0.25)
#' @param snp6_reference_info_file Reference info file for SNP6 data (Default: NA)
#' @param enhanced_grid_search Flag to determine if the grid search should be performed with a higher number of steps (Default: FALSE)
#' @param usebeagle Logical, if TRUE, expects Beagle output (VCF) in impute_results_dir and converts to IMPUTE format (Default: FALSE)
#' @param verbose_logging Print out more information during the run
#' (Default: FALSE)
#' @param skip_preprocessing Boolean, if TRUE skips the initial allele counting and GC correction (Default: FALSE)
#' @param preprocessed_data_dir Directory where existing preprocessed .tab files are located. If provided and skip_preprocessing is TRUE, files will be copied to local directory. (Default: NA)
#' @param logging_path Path to write log files to (Default: ".")
#'
#' @useDynLib Battenberg, .registration = TRUE
#' @importFrom data.table :=
#' @author sd11, jdemeul, Naser Ansari-Pour, Julio Cesar Cortes Rios
#' @export
battenberg <- function(
  analysis = "paired",
  samplename,
  normalname,
  sample_data_file,
  normal_data_file,
  imputeinfofile,
  g1000prefix,
  problemloci,
  allele_counts_dir,
  impute_results_dir,
  gccorrectprefix = NULL,
  repliccorrectprefix = NULL,
  g1000allelesprefix = NA,
  ismale = NA,
  data_type = "wgs",
  nthreads = 8,
  platform_gamma = 1,
  phasing_gamma = 1,
  segmentation_gamma = 10,
  segmentation_kmin = 3,
  phasing_kmin = 1,
  clonality_dist_metric = 0,
  ascat_dist_metric = 1,
  min_ploidy = 1.6,
  max_ploidy = 4.8,
  min_rho = 0.1,
  max_rho = 1.0,
  min_goodness = 0.63,
  uninformative_baf_threshold = 0.51,
  min_normal_depth = 10,
  min_base_qual = 20,
  min_map_qual = 35,
  max_allowed_state = 250,
  cn_upper_limit = 1000,
  calc_seg_baf_option = 3,
  externalhaplotypefile = NA,
  write_battenberg_phasing = TRUE,
  multisample_relative_weight_balanced = 0.25,
  multisample_maxlag = 90,
  segmentation_gamma_multisample = 5,
  snp6_reference_info_file = NA,
  apt_probeset_genotype_exe = "apt-probeset-genotype",
  apt_probeset_summarize_exe = "apt-probeset-summarize",
  norm_geno_clust_exe = "normalize_affy_geno_cluster.pl",
  birdseed_report_file = "birdseed.report.txt",
  heterozygous_filter = "none",
  prior_breakpoints_file = NULL,
  genomebuild = "hg38",
  chrom_coord_file = NULL,
  enhanced_grid_search = FALSE,
  verbose_logging = FALSE,
  usebeagle = FALSE,
  skip_preprocessing = FALSE,
  preprocessed_data_dir = NA,
  logging_path = "."
) {
  libs <- .libPaths()

  # Set global thread limits based on user configuration
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(nthreads)
  }
  Sys.setenv(OMP_NUM_THREADS = nthreads)
  Sys.setenv(MKL_NUM_THREADS = nthreads)
  Sys.setenv(OPENBLAS_NUM_THREADS = nthreads)

  log_setup(logging_path, verbose_logging)

  # Inform the user about the thread configuration
  log_info(strrep("-", 60))
  log_info("Battenberg Thread Configuration:")
  log_info("  - Total thread budget: {nthreads}")
  log_info("  - The pipeline will dynamically allocate these cores between")
  log_info("    sample-level and logic-level parallelism.")
  log_info(strrep("-", 60))

  log_info("Starting analysis for {samplename}")


  if (analysis == "cell_line") {
    calc_seg_baf_option <- 1
    phasing_gamma <- 1
    phasing_kmin <- 2
    segmentation_gamma <- 20
    segmentation_kmin <- 3
    # no matched normal required, but we  are
    # generating normal counts which have this name coded
    normalname <- paste0(samplename, "_normal")
    # other cell_line specific parameter values
    min_ploidy <- min_ploidy
    max_ploidy <- max_ploidy
    min_rho <- 0.99
    max_rho <- 1.01
  }
  if (analysis == "germline") {
    calc_seg_baf_option <- 1
    phasing_gamma <- 3
    phasing_kmin <- 1
    segmentation_gamma <- 3
    segmentation_kmin <- 3
    # no matched normal required,
    # but we  are generating normal counts which have this name coded
    normalname <- paste0(samplename, "_normal")
    min_ploidy <- 1.5
    max_ploidy <- 2.5
    min_rho <- 0.99
    max_rho <- 1.01
  }

  if (data_type == "wgs" && is.na(ismale)) {
    log_failure("Please provide a boolean denominator whether \\
    this sample represents a male donor")
  }

  if (data_type == "wgs" && is.na(g1000allelesprefix)) {
    log_failure("Please provide a path to 1000 Genomes allele reference files")
  }

  if (data_type == "wgs" && is.null(gccorrectprefix)) {
    log_failure("Please provide a path to GC content reference files")
  }

  if (data_type == "wgs" && !file.exists(problemloci)) {
    log_failure("Please provide a path to a problematic loci file")
  }

  if (!file.exists(imputeinfofile)) {
    log_failure("Please provide a path to an impute info file")
  }

  # check whether the impute_info.txt file contains correct paths
  # check whether the impute_info.txt file contains correct paths
  check_imputeinfofile(
    imputeinfofile = imputeinfofile,
    is_male = ismale,
    usebeagle = usebeagle
  )

  # check whether multisample case
  nsamples <- length(samplename)

  if (data_type == "wgs" || data_type == "WGS") {
    if (nsamples > 1) {
      log_info("Running Battenberg in multisample mode on {nsamples} samples: \\
                {paste(samplename, collapse = ', ')}")
    }
    chrom_names <- get_chrom_names(imputeinfofile, ismale, analysis = analysis)
  } else if (data_type == "snp6" || data_type == "SNP6") {
    if (nsamples > 1) {
      log_failure("Battenberg multisample mode has \\
       not been tested with SNP6 data")
    }
    chrom_names <- get_chrom_names(imputeinfofile, TRUE)
  }
  # Global parameter validation
  if (!missing(allele_counts_dir) && !is.na(allele_counts_dir) && !dir.exists(allele_counts_dir)) {
    log_failure("allele_counts_dir does not exist: {allele_counts_dir}")
  }
  if (!missing(impute_results_dir) && !is.na(impute_results_dir) && !dir.exists(impute_results_dir)) {
    log_failure("impute_results_dir does not exist: {impute_results_dir}")
  }

  log_info(chrom_names)
  for (sampleidx in 1:nsamples) {
    if (data_type == "wgs" || data_type == "WGS") {
      # Setup for parallel computing
      if (nthreads > 1 && !skip_preprocessing) {
        # In preprocessing, we run samples sequentially in a for loop.
        # So each sample can use the FULL nthreads budget for chromosome-level parallelism.
        clp <- parallel::makeCluster(nthreads, outfile = "")
        doParallel::registerDoParallel(clp)
      }

      if (!skip_preprocessing) {
        if (analysis == "paired") {
          if (is.null(normalname) || is.na(normalname)) {
            log_failure("No normal sample is specified for \\
                'paired analysis' - a normal paired BAM is required")
          }
          prepare_wgs(
            chrom_names = chrom_names,
            tumourbam = sample_data_file[sampleidx],
            normalbam = normal_data_file,
            tumourname = samplename[sampleidx],
            normalname = normalname,
            g1000allelesprefix = g1000allelesprefix,
            g1000prefix = g1000prefix,
            gccorrectprefix = gccorrectprefix,
            repliccorrectprefix = repliccorrectprefix,
            min_base_qual = min_base_qual,
            min_map_qual = min_map_qual,
            allele_counts_dir = allele_counts_dir,
            min_normal_depth = min_normal_depth,
            nthreads = nthreads,
            libs = libs
          )
        } else if (analysis == "cell_line") {
          prepare_wgs_cell_line(
            chrom_names = chrom_names,
            chrom_coord = chrom_coord_file,
            tumourbam = sample_data_file[sampleidx],
            tumourname = samplename[sampleidx],
            g1000lociprefix = g1000prefix,
            g1000allelesprefix = g1000allelesprefix,
            gamma_ivd = 1e5,
            kmin_ivd = 50,
            centromere_noise_seg_size = 1e6,
            centromere_dist = 5e5,
            min_het_dist = 1e5,
            gamma_logr = 100,
            length_adjacent = 5e4,
            gccorrectprefix = gccorrectprefix,
            repliccorrectprefix = repliccorrectprefix,
            min_base_qual = min_base_qual,
            min_map_qual = min_map_qual,
            allele_counts_dir = allele_counts_dir,
            min_normal_depth = min_normal_depth,
            libs = libs
          )
        } else if (analysis == "germline") {
          prepare_wgs_germline(
            chrom_names = chrom_names,
            chrom_coord = chrom_coord_file,
            germlinebam = sample_data_file[sampleidx],
            germlinename = samplename[sampleidx],
            g1000lociprefix = g1000prefix,
            g1000allelesprefix = g1000allelesprefix,
            gamma_ivd = 1e5,
            kmin_ivd = 50,
            centromere_noise_seg_size = 1e6,
            centromere_dist = 5e5,
            min_het_dist = 2e3,
            gamma_logr = 100,
            length_adjacent = 5e4,
            gccorrectprefix = gccorrectprefix,
            repliccorrectprefix = repliccorrectprefix,
            min_base_qual = min_base_qual,
            min_map_qual = min_map_qual,
            allele_counts_dir = allele_counts_dir,
            min_normal_depth = min_normal_depth,
            libs = libs
          )
        }
      } else {
        log_info("Skipping preprocessing (allele counting and GC correction) for sample '{samplename[sampleidx]}'")

        # If a preprocessed directory is provided, copy the files to current working directory
        if (!is.na(preprocessed_data_dir) && dir.exists(preprocessed_data_dir)) {
          log_info("Providing existing preprocessed files from {preprocessed_data_dir}")

          files_to_copy <- c(
            paste0(samplename[sampleidx], "_mutantBAF.tab"),
            paste0(samplename[sampleidx], "_normalBAF.tab"),
            paste0(samplename[sampleidx], "_mutantLogR.tab"),
            paste0(samplename[sampleidx], "_normalLogR.tab"),
            paste0(samplename[sampleidx], "_alleleCounts.tab"),
            paste0(samplename[sampleidx], "_mutantLogR_gcCorrected.tab"),
            paste0(samplename[sampleidx], "_GCwindowCorrelations.txt")
          )

          # Also copy allele frequency files if they exist there, as they are needed for haplotyping
          freq_files <- list.files(preprocessed_data_dir, pattern = paste0("^", samplename[sampleidx], "_alleleFrequencies_chr.*\\.txt$"))
          files_to_copy <- c(files_to_copy, freq_files)

          for (f in files_to_copy) {
            src <- file.path(preprocessed_data_dir, f)
            if (file.exists(src)) {
              log_info("Copying {f} to current directory")
              file.copy(src, ".", overwrite = TRUE)
            } else if (!grepl("gcCorrected|Correlations", f)) {
              # Some files might be optional or missing depending on analysis mode,
              # but essential ones should be warned about
              log_warning("Expected preprocessed file {f} not found in {preprocessed_data_dir}")
            }
          }
        }
      }

      # Kill the threads
      if (nthreads > 1 && !skip_preprocessing) {
        parallel::stopCluster(clp)
      }
    } else if (data_type == "snp6" || data_type == "SNP6") {
      prepare_snp6(
        tumour_cel_file = sample_data_file[sampleidx],
        normal_cel_file = normal_data_file,
        tumourname = samplename[sampleidx],
        chrom_names = chrom_names,
        snp6_reference_info_file = snp6_reference_info_file,
        apt_probeset_genotype_exe = apt_probeset_genotype_exe,
        apt_probeset_summarize_exe = apt_probeset_summarize_exe,
        norm_geno_clust_exe = norm_geno_clust_exe,
        birdseed_report_file = birdseed_report_file,
        genomebuild = genomebuild
      )
    } else {
      log_failure("Unknown data type provided, please provide wgs or snp6")
      q(save = "no", status = 1)
    }

    # Removed } else (end of if !skip_preprocessing) as skipping logic is now handled by presence of directories/files inside prepare functions or removed entirely.


    if (data_type == "snp6" || data_type == "SNP6") {
      # Infer what the gender is - WGS requires it to be specified
      gender <- infer_gender_birdseed(birdseed_report_file)
      ismale <- gender == "male"
    }


    if (TRUE) {
      # if external phasing data is provided (as a vcf), split into chromosomes for use in haplotype reconstruction
      if (!is.na(externalhaplotypefile) && file.exists(externalhaplotypefile)) {
        externalhaplotypeprefix <- paste0(normalname, "_external_haplotypes_chr")

        # if these files exist already, no need to split again
        if (any(!file.exists(paste0(externalhaplotypeprefix, seq_along(chrom_names), ".vcf")))) {
          log_info("Splitting external phasing data from '{externalhaplotypefile}'")
          split_input_haplotypes(
            chrom_names = chrom_names,
            externalhaplotypefile = externalhaplotypefile,
            outprefix = externalhaplotypeprefix
          )
        } else {
          log_info("No need to split, external haplotype files per chromosome found")
        }
      } else {
        externalhaplotypeprefix <- NA
      }

      # Setup for parallel computing
      # Setup for parallel computing
      if (nthreads > 1) {
        clp <- parallel::makeCluster(nthreads, outfile = "")
        doParallel::registerDoParallel(clp)
      }

      # Reconstruct haplotypes
      # mclapply(seq_along(chrom_names), function(chrom) {
      do_haplotyping <- function(i) {
        .libPaths(libs)
        chrom <- chrom_names[i]
        if (analysis == "germline") {
          log_info("germline chrom {chrom}")
          run_haplotyping_germline(
            chrom = chrom,
            germlinename = samplename[sampleidx],
            normalname = normalname,
            ismale = ismale,
            imputeinfofile = imputeinfofile,
            problemloci = problemloci,
            impute_results_dir = impute_results_dir,
            min_normal_depth = min_normal_depth,
            chrom_names = chrom_names,
            snp6_reference_info_file = NA,
            heterozygous_filter = NA,
            usebeagle = usebeagle
          )
        } else {
          .libPaths(libs)
          chrom <- chrom_names[i]
          log_info("chrom {chrom}")
          run_haplotyping(
            chrom = chrom,
            tumourname = samplename[sampleidx],
            normalname = normalname,
            ismale = ismale,
            imputeinfofile = imputeinfofile,
            problemloci = problemloci,
            impute_results_dir = impute_results_dir,
            min_normal_depth = min_normal_depth,
            chrom_names = chrom_names,
            snp6_reference_info_file = snp6_reference_info_file,
            heterozygous_filter = heterozygous_filter,
            externalhaplotypeprefix = externalhaplotypeprefix,
            usebeagle = usebeagle
          )
        }
      }
      run_parallel_or_serial(
        iterator = seq_along(chrom_names),
        func = do_haplotyping,
        libs = libs
      )

      # Kill the threads as from here its all single core
      # Kill the threads as from here its all single core
      if (nthreads > 1) {
        parallel::stopCluster(clp)
      }

      # Combine all the BAF output into a single file
      concatenate_baf_files(
        input_start = paste(samplename[sampleidx], "_chr", sep = ""),
        input_end = "_heterozygousMutBAFs_haplotyped.txt",
        output_file = paste(samplename[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
        chr_names = chrom_names
      )
    }

    # Segment the phased and haplotyped BAF data
    segment_baf_phased(
      samplename = samplename[sampleidx],
      inputfile = paste(samplename[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
      outputfile = paste(samplename[sampleidx], ".BAFsegmented.txt", sep = ""),
      prior_breakpoints_file = prior_breakpoints_file,
      gamma = segmentation_gamma,
      phasegamma = phasing_gamma,
      kmin = segmentation_kmin,
      phasekmin = phasing_kmin,
      calc_seg_baf_option = calc_seg_baf_option
    )

    if (nsamples > 1 || write_battenberg_phasing) {
      # Write the Battenberg phasing information to disk as a vcf
      write_battenberg_phasing(
        tumourname = samplename[sampleidx],
        SNPfiles = paste0(
          samplename[sampleidx], "_alleleFrequencies_chr",
          chrom_names, ".txt"
        ),
        imputedHaplotypeFiles = paste0(
          samplename[sampleidx],
          "_impute_output_chr", chrom_names,
          "_allHaplotypeInfo.txt"
        ),
        bafsegmented_file = paste0(samplename[sampleidx], ".BAFsegmented.txt"),
        outprefix = paste0(samplename[sampleidx], "_Battenberg_phased_chr"),
        chrom_names = chrom_names,
        include_homozygous = FALSE
      )
    }
  }

  # if this is a multisample run, combine the battenberg phasing outputs, incorporate it and resegment
  if (nsamples > 1) {
    log_info("Constructing multisample phasing")
    multisamplehaplotypeprefix <- paste0(normalname, "_multisample_haplotypes_chr")


    if (nthreads > 1) {
      clp <- parallel::makeCluster(nthreads, outfile = "")
      doParallel::registerDoParallel(clp)
    }

    run_parallel_or_serial(seq_along(chrom_names), function(i) {
      chrom <- chrom_names[i]
      log_info("multisample phasing chrom {chrom}")

      get_multisample_phasing(
        chrom = chrom,
        bbphasingprefixes = paste(samplename, "_Battenberg_phased_chr", sep = ""),
        maxlag = multisample_maxlag,
        relative_weight_balanced = multisample_relative_weight_balanced,
        outprefix = multisamplehaplotypeprefix
      )
    }, libs)

    # continue over all samples to incorporate the multisample phasing
    for (sampleidx in 1:nsamples) {
      # rename the original files without multisample phasing info
      MutBAFfiles <- paste0(samplename[sampleidx], "_chr", chrom_names, "_heterozygousMutBAFs_haplotyped.txt")
      heterozygousdatafiles <- paste0(samplename[sampleidx], "_chr", chrom_names, "_heterozygousData.png")
      raffiles <- paste0(samplename[sampleidx], "_RAFseg_chr", chrom_names, ".png")
      segfiles <- paste0(samplename[sampleidx], "_segment_chr", chrom_names, ".png")
      haplotypedandbafsegmentedfiles <- paste0(samplename[sampleidx], c("_heterozygousMutBAFs_haplotyped.txt", ".BAFsegmented.txt"))

      file.copy(
        from = MutBAFfiles,
        to = gsub(
          pattern = ".txt$", replacement = "_noMulti.txt",
          x = MutBAFfiles
        ), overwrite = TRUE
      )
      file.copy(
        from = heterozygousdatafiles,
        to = gsub(
          pattern = ".png$", replacement = "_noMulti.png",
          x = heterozygousdatafiles
        ), overwrite = TRUE
      )
      file.copy(
        from = raffiles,
        to = gsub(
          pattern = ".png$", replacement = "_noMulti.png",
          x = raffiles
        ), overwrite = TRUE
      )
      file.copy(
        from = segfiles,
        to = gsub(
          pattern = ".png$", replacement = "_noMulti.png",
          x = segfiles
        ), overwrite = TRUE
      )
      file.copy(
        from = haplotypedandbafsegmentedfiles,
        to = gsub(
          pattern = ".txt$", replacement = "_noMulti.txt",
          x = haplotypedandbafsegmentedfiles
        ), overwrite = TRUE
      )
      # done renaming, next sections will overwrite orignals

      run_parallel_or_serial(seq_along(chrom_names), function(i) {
        chrom <- chrom_names[i]
        log_info("sample in nsamples chrom {chrom}")

        # Reconstruct haplotypes from external file
        input_known_haplotypes(
          chrom = chrom,
          chrom_names = chrom_names,
          imputedHaplotypeFile = paste(samplename[sampleidx],
            "_impute_output_chr", chrom,
            "_allHaplotypeInfo.txt",
            sep = ""
          ),
          externalHaplotypeFile = paste(multisamplehaplotypeprefix, chrom,
            ".vcf",
            sep = ""
          ),
          oldfilesuffix = "_noMulti.txt"
        )

        # Get BAFs for the specific chromosome
        GetChromosomeBAFs(
          chrom = chrom,
          SNP_file = paste(samplename[sampleidx], "_alleleFrequencies_chr",
            chrom, ".txt",
            sep = ""
          ),
          haplotypeFile = paste(samplename[sampleidx], "_impute_output_chr",
            chrom, "_allHaplotypeInfo.txt",
            sep = ""
          ),
          samplename = samplename[sampleidx],
          outfile = paste(samplename[sampleidx], "_chr", chrom,
            "_heterozygousMutBAFs_haplotyped.txt",
            sep = ""
          ),
          chr_names = chrom_names,
          minCounts = min_normal_depth
        )

        # Plot the intermediate results
        plot_haplotype_data(
          haplotyped_baf_file = paste(samplename[sampleidx], "_chr", chrom,
            "_heterozygousMutBAFs_haplotyped.txt",
            sep = ""
          ),
          image_file_name = paste(samplename[sampleidx], "_chr", chrom,
            "_heterozygousData.png",
            sep = ""
          ),
          samplename = samplename[sampleidx],
          chrom = chrom
        )
      }, libs)
    }

    # Kill the threads as from here its single core
    # Kill the threads as from here its single core
    if (nthreads > 1) {
      parallel::stopCluster(clp)
    }

    for (sampleidx in 1:nsamples) {
      # Combine all the BAF output into a single file
      concatenate_baf_files(
        input_start = paste0(samplename[sampleidx], "_chr"),
        input_end = "_heterozygousMutBAFs_haplotyped.txt",
        output_file = paste0(samplename[sampleidx], "_heterozygousMutBAFs_haplotyped.txt"),
        chr_names = chrom_names
      )
    }
    # Segment the phased and haplotyped BAF data
    segment_baf_phased_multisample(
      samplename = samplename,
      inputfile = paste(samplename, "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
      outputfile = paste(samplename, ".BAFsegmented.txt", sep = ""),
      prior_breakpoints_file = prior_breakpoints_file,
      gamma = segmentation_gamma_multisample,
      calc_seg_baf_option = calc_seg_baf_option,
      GENOMEBUILD = genomebuild
    )
  }

  # Setup for parallel computing
  # Setup for parallel computing
  if (nthreads > 1) {
    # Dynamic Budgeting: Divide total nthreads by the number of samples being run in parallel.
    # If we have 40 cores and 2 samples, each sample gets 20 cores (inner_threads).
    # If we have more samples than cores, each sample gets 1 core.
    num_sample_workers <- min(nsamples, nthreads)
    clp <- parallel::makeCluster(num_sample_workers, outfile = "")
    doParallel::registerDoParallel(clp)
  }

  # Use the universal helper to process each sample
  run_parallel_or_serial(seq_len(nsamples), function(sampleidx) {
    # Scoping ensures this function sees 'samplename', 'libs', etc.
    log_info("Fitting final copy number and calling subclones for sample '{samplename[sampleidx]}'")

    # Determine file paths based on data type and analysis mode
    if (data_type == "wgs" || data_type == "WGS") {
      logr_file <- paste(samplename[sampleidx], "_mutantLogR_gcCorrected.tab", sep = "")
      if (analysis == "paired") {
        allelecounts_file <- paste(samplename[sampleidx], "_alleleCounts.tab", sep = "")
      } else {
        allelecounts_file <- NULL
      }
    }

    # Calculate safe inner threads to avoid thrashing
    # If NO parallel grid search, force sequential execution
    inner_threads <- max(1, floor(nthreads / min(nsamples, nthreads)))
    log_info(
      "Dynamic Threading: budget={nthreads}, workers={min(nsamples, nthreads)} -> inner_threads={inner_threads} (per sample)"
    )
    # Parallel workers will now report their index and error details if they fail
    fit_copy_number(
      samplename = samplename[sampleidx],
      outputfile_prefix = paste(samplename[sampleidx], "_", sep = ""),
      inputfile_baf_segmented = paste(samplename[sampleidx], ".BAFsegmented.txt", sep = ""),
      inputfile_baf = paste(samplename[sampleidx], "_mutantBAF.tab", sep = ""),
      inputfile_logr = logr_file,
      dist_choice = clonality_dist_metric,
      ascat_dist_choice = ascat_dist_metric,
      min_ploidy = min_ploidy,
      max_ploidy = max_ploidy,
      min_rho = min_rho,
      max_rho = max_rho,
      min_goodness = min_goodness,
      uninformative_baf_threshold = uninformative_baf_threshold,
      gamma_param = platform_gamma,
      use_preset_rho_psi = FALSE,
      preset_rho = NA,
      preset_psi = NA,
      read_depth = 30,
      analysis = analysis,
      nthreads = inner_threads,
      enhanced_grid_search = enhanced_grid_search
    )

    # Fit a second CN state (subclonal)
    log_info("call_subclones")
    call_subclones(
      sample_name = samplename[sampleidx],
      baf_segmented_file = paste(samplename[sampleidx], ".BAFsegmented.txt", sep = ""),
      logr_file = logr_file,
      rho_psi_file = paste(samplename[sampleidx], "_rho_and_psi.txt", sep = ""),
      output_file = paste(samplename[sampleidx], "_copynumber.txt", sep = ""),
      output_figures_prefix = paste(samplename[sampleidx], "_subclones_chr",
        sep = ""
      ),
      output_gw_figures_prefix = paste(samplename[sampleidx],
        "_BattenbergProfile",
        sep = ""
      ),
      masking_output_file = paste(samplename[sampleidx],
        "_segment_masking_details.txt",
        sep = ""
      ),
      prior_breakpoints_file = prior_breakpoints_file,
      chr_names = chrom_names,
      gamma = platform_gamma,
      segmentation_gamma = NA,
      siglevel = 0.05,
      maxdist = 0.01,
      max_allowed_state = max_allowed_state,
      nthreads = inner_threads,
      cn_upper_limit = cn_upper_limit,
      noperms = 1000,
      calc_seg_baf_option = calc_seg_baf_option,
      verbose_logging = verbose_logging
    )

    # Handle Male ChrX if applicable
    if (ismale && "X" %in% chrom_names) {
      log_info("callChrXsubclones")
      callChrXsubclones(
        tumourname = samplename[sampleidx],
        X_gamma = 1000,
        X_kmin = 100,
        genomebuild = genomebuild,
        AR = TRUE,
        prior_breakpoints_file = prior_breakpoints_file,
        chrom_names = chrom_names,
        data_type = data_type
      )
    }

    # Cleanup/Post-hoc visualisations
    log_info("make_posthoc_plots")
    make_posthoc_plots(
      samplename = samplename[sampleidx],
      logr_file = logr_file,
      bafsegmented_file = paste(samplename[sampleidx], ".BAFsegmented.txt", sep = ""),
      logrsegmented_file = paste(samplename[sampleidx], ".logRsegmented.txt", sep = ""),
      allelecounts_file = allelecounts_file
    )

    # Generate refit suggestions
    log_info("cnfit_to_refit_suggestions")
    cnfit_to_refit_suggestions(
      samplename = samplename[sampleidx],
      subclones_file = paste(samplename[sampleidx], "_copynumber_extended.txt", sep = ""),
      rho_psi_file = paste(samplename[sampleidx], "_rho_and_psi.txt", sep = ""),
      gamma_param = platform_gamma
    )
  }, libs)

  # Kill the threads as last part again is single core
  # Kill the threads as last part again is single core
  if (nthreads > 1) {
    parallel::stopCluster(clp)
  }

  if (nsamples > 1) {
    log_info("Assessing mirrored subclonal allelic imbalance (MSAI)")
    call_multisample_MSAI(
      rdsprefix = multisamplehaplotypeprefix,
      subclonesfiles = paste0(samplename, "_copynumber_extended.txt"),
      chrom_names = chrom_names,
      tumournames = samplename,
      plotting = TRUE
    )
  }
}
