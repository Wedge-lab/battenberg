#' Battenberg Command Line Interface
#' @description Parses command line arguments and executes the main battenberg function.
#' @export
battenberg_cli <- function() {
  options(error = function() {
    # Get the raw calls
    calls <- sys.calls()

    msg <- sprintf("Fatal Error: %s\n\n--- Call Stack ---", geterrmessage())
    for (i in seq_along(calls)) {
      msg <- paste(msg, sprintf("[%2d] %s", i, deparse(calls[[i]], width.cutoff = 500)[1]), sep = "\n")
    }
    log_failure("{msg}")
    quit(save = "no", status = 1)
  })
  options(show.error.messages = TRUE)
  options(keep.source = TRUE)
  options(width = 10000)

  option_list <- list(
    # Core Analysis & Sample Info
    optparse::make_option(c("-a", "--analysis"),
      type = "character", default = "paired",
      help = "Analysis type: paired, cell_line, germline"
    ),
    optparse::make_option(c("-t", "--samplename"),
      type = "character",
      help = "Tumour/Sample identifier"
    ),
    optparse::make_option(c("-n", "--normalname"),
      type = "character",
      help = "Matched normal identifier"
    ),
    optparse::make_option(c("--sample_data_file"),
      type = "character",
      help = "BAM/CEL for sample"
    ),
    optparse::make_option(c("--normal_data_file"),
      type = "character",
      help = "BAM/CEL for normal"
    ),
    optparse::make_option(c("--ismale"),
      type = "logical",
      default = NA,
      help = "TRUE/FALSE for donor sex"
    ),

    # Reference Paths
    optparse::make_option(c("--imputeinfofile"),
      type = "character",
      help = "Path to impute info file"
    ),
    optparse::make_option(c("--g1000prefix"),
      type = "character",
      help = "Prefix for 1000G SNP loci"
    ),
    optparse::make_option(c("--g1000allelesprefix"),
      type = "character",
      default = NA,
      help = "Prefix for 1000G alleles"
    ),
    optparse::make_option(c("--gccorrectprefix"),
      type = "character",
      default = NULL,
      help = "Prefix for GC correction"
    ),
    optparse::make_option(c("--repliccorrectprefix"),
      type = "character",
      default = NULL,
      help = "Prefix for replication timing"
    ),
    optparse::make_option(c("--problemloci"),
      type = "character",
      help = "Path to problem loci file"
    ),
    optparse::make_option(c("--genomebuild"),
      type = "character",
      default = "hg38",
      help = "hg19 or hg38"
    ),
    optparse::make_option(c("--chrom_coord_file"),
      type = "character",
      default = NULL
    ),
    optparse::make_option(c("--allele_counts_dir"),
      type = "character", default = NA,
      help = "Directory containing pre-calculated allele counts"
    ),
    optparse::make_option(c("--impute_results_dir"),
      type = "character", default = NA,
      help = "Directory containing pre-calculated imputation results"
    ),

    # Executables & Hardware
    optparse::make_option(c("--nthreads"),
      type = "integer", default = 8
    ),
    optparse::make_option(c("--data_type"),
      type = "character", default = "wgs"
    ),

    # Tuning Parameters (Gamma & Kmin)
    optparse::make_option(c("--platform_gamma"),
      type = "double", default = 1
    ),
    optparse::make_option(c("--phasing_gamma"),
      type = "double", default = 1
    ),
    optparse::make_option(c("--segmentation_gamma"),
      type = "double", default = 10
    ),
    optparse::make_option(c("--segmentation_gamma_multisample"),
      type = "double", default = 5
    ),
    optparse::make_option(c("--segmentation_kmin"),
      type = "integer", default = 3
    ),
    optparse::make_option(c("--phasing_kmin"),
      type = "integer", default = 1
    ),

    # Grid Search / ASCAT Params
    optparse::make_option(c("--clonality_dist_metric"),
      type = "integer", default = 0
    ),
    optparse::make_option(c("--ascat_dist_metric"),
      type = "integer", default = 1
    ),
    optparse::make_option(c("--min_ploidy"),
      type = "double", default = 1.6
    ),
    optparse::make_option(c("--max_ploidy"),
      type = "double", default = 4.8
    ),
    optparse::make_option(c("--min_rho"),
      type = "double", default = 0.1
    ),
    optparse::make_option(c("--max_rho"),
      type = "double", default = 1.0
    ),
    optparse::make_option(c("--min_goodness"),
      type = "double", default = 0.63
    ),
    optparse::make_option(c("--uninformative_baf_threshold"),
      type = "double", default = 0.51
    ),
    optparse::make_option(c("--enhanced_grid_search"),
      type = "logical", default = FALSE, action = "store_true"
    ),
    optparse::make_option(c("--skip_preprocessing"),
      type = "logical", default = FALSE, action = "store_true"
    ),
    optparse::make_option(c("--preprocessed_data_dir"),
      type = "character", default = NA
    ),

    # Quality Thresholds
    optparse::make_option(c("--min_normal_depth"),
      type = "integer", default = 10
    ),
    optparse::make_option(c("--min_base_qual"),
      type = "integer", default = 20
    ),
    optparse::make_option(c("--min_map_qual"),
      type = "integer", default = 35
    ),
    optparse::make_option(c("--max_allowed_state"),
      type = "integer", default = 250
    ),
    optparse::make_option(c("--cn_upper_limit"),
      type = "integer", default = 1000
    ),
    optparse::make_option(c("--calc_seg_baf_option"),
      type = "integer", default = 3
    ),

    # Beagle Specifics
    optparse::make_option(c("--usebeagle"),
      type = "logical", default = FALSE, action = "store_true"
    ),
    optparse::make_option(c("--prior_breakpoints_file"),
      type = "character", default = NULL
    ),
    optparse::make_option(c("--externalhaplotypefile"),
      type = "character", default = NA
    ),
    optparse::make_option(c("--write_battenberg_phasing"),
      type = "logical", default = TRUE
    ),

    # Multisample & SNP6 Legacy/Special
    optparse::make_option(c("--multisample_maxlag"),
      type = "integer",
      default = 90
    ),
    optparse::make_option(c("--multisample_relative_weight_balanced"),
      type = "double", default = 0.25
    ),
    optparse::make_option(c("--snp6_reference_info_file"),
      type = "character", default = NA
    ),
    optparse::make_option(c("--apt_probeset_genotype_exe"),
      type = "character", default = "apt-probeset-genotype"
    ),
    optparse::make_option(c("--apt_probeset_summarize_exe"),
      type = "character", default = "apt-probeset-summarize"
    ),
    optparse::make_option(c("--norm_geno_clust_exe"),
      type = "character", default = "normalize_affy_geno_cluster.pl"
    ),
    optparse::make_option(c("--birdseed_report_file"),
      type = "character", default = "birdseed.report.txt"
    ),
    optparse::make_option(c("--heterozygous_filter"),
      type = "character", default = "none"
    ),

    # Logging & Debug
    optparse::make_option(c("--verbose_logging"),
      type = "logical",
      default = FALSE, action = "store_true"
    ),
    optparse::make_option(c("--logging_path"),
      type = "character", default = "."
    )
  )

  # Parse arguments
  parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(parser)

  # Remove the 'help' flag which optparse adds automatically
  opt$help <- NULL

  log_info(strrep("=", 60))
  log_info("BATTENBERG CLI: EXECUTION PARAMETERS")
  log_info(strrep("=", 60))

  # Sort names so they are easy to find in the log
  opt_names <- sort(names(opt))
  for (name in opt_names) {
    # Cleanly format each argument and its value
    val <- opt[[name]]
    log_info(sprintf("%-40s : %s", name, paste(val, collapse = ", ")))
  }
  log_info(strrep("=", 60))

  # Execute main function
  do.call(battenberg, opt)
}
