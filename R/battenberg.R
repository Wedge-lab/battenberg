
#' Run the Battenberg pipeline
#' 
#' @param tumourname Tumour identifier, this is used as a prefix for the output files. If allele counts are supplied separately, they are expected to have this identifier as prefix.
#' @param normalname Matched normal identifier, this is used as a prefix for the output files. If allele counts are supplied separately, they are expected to have this identifier as prefix.
#' @param tumour_data_file A BAM or CEL file for the tumour
#' @param normal_data_file A BAM or CEL file for the normal
#' @param imputeinfofile Full path to a Battenberg impute info file with pointers to Impute2 reference data
#' @param g1000prefix Full prefix path to 1000 Genomes SNP loci data, as part of the Battenberg reference data
#' @param problemloci Full path to a problem loci file that contains SNP loci that should be filtered out
#' @param gccorrectprefix Full prefix path to GC content files, as part of the Battenberg reference data, not required for SNP6 data (Default: NULL)
#' @param repliccorrectprefix Full prefix path to replication timing files, as part of the Battenberg reference data, not required for SNP6 data (Default: NULL)
#' @param g1000allelesprefix Full prefix path to 1000 Genomes SNP alleles data, as part of the Battenberg reference data, not required for SNP6 data (Default: NA)
#' @param ismale A boolean set to TRUE if the donor is male, set to FALSE if female, not required for SNP6 data (Default: NA)
#' @param data_type String that contains either wgs or snp6 depending on the supplied input data (Default: wgs)
#' @param impute_exe Pointer to the Impute2 executable (Default: impute2, i.e. expected in $PATH)
#' @param allelecounter_exe Pointer to the alleleCounter executable (Default: alleleCounter, i.e. expected in $PATH)
#' @param nthreads The number of concurrent processes to use while running the Battenberg pipeline (Default: 8)
#' @param platform_gamma Platform scaling factor, suggestions are set to 1 for wgs and to 0.55 for snp6 (Default: 1)
#' @param phasing_gamma Gamma parameter used when correcting phasing mistakes (Default: 1)
#' @param segmentation_gamma The gamma parameter controls the size of the penalty of starting a new segment during segmentation. It is therefore the key parameter for controlling the number of segments (Default: 10)
#' @param segmentation_kmin Kmin represents the minimum number of probes/SNPs that a segment should consist of (Default: 3)
#' @param phasing_kmin Kmin used when correcting for phasing mistakes (Default: 3)
#' @param clonality_dist_metric  Distance metric to use when choosing purity/ploidy combinations (Default: 0)
#' @param ascat_dist_metric Distance metric to use when choosing purity/ploidy combinations (Default: 1)
#' @param min_ploidy Minimum ploidy to be considered (Default: 1.6)
#' @param max_ploidy Maximum ploidy to be considered (Default: 4.8)
#' @param min_rho Minimum purity to be considered (Default: 0.1)
#' @param min_goodness Minimum goodness of fit required for a purity/ploidy combination to be accepted as a solution (Default: 0.63)
#' @param uninformative_BAF_threshold The threshold beyond which BAF becomes uninformative (Default: 0.51)
#' @param min_normal_depth Minimum depth required in the matched normal for a SNP to be considered as part of the wgs analysis (Default: 10)
#' @param min_base_qual Minimum base quality required for a read to be counted when allele counting (Default: 20)
#' @param min_map_qual Minimum mapping quality required for a read to be counted when allele counting (Default: 35)
#' @param calc_seg_baf_option Sets way to calculate BAF per segment: 1=mean, 2=median, 3=ifelse median==0 | 1, mean, median (Default: 3)
#' @param skip_allele_counting Provide TRUE when allele counting can be skipped (i.e. its already done) (Default: FALSE)
#' @param skip_preprocessing Provide TRUE when preprocessing is already complete (Default: FALSE)
#' @param skip_phasing  Provide TRUE when phasing is already complete (Default: FALSE)
#' @param snp6_reference_info_file Reference files for the SNP6 pipeline only (Default: NA)
#' @param apt.probeset.genotype.exe Helper tool for extracting data from CEL files, SNP6 pipeline only (Default: apt-probeset-genotype)
#' @param apt.probeset.summarize.exe  Helper tool for extracting data from CEL files, SNP6 pipeline only (Default: apt-probeset-summarize)
#' @param norm.geno.clust.exe  Helper tool for extracting data from CEL files, SNP6 pipeline only (Default: normalize_affy_geno_cluster.pl)
#' @param birdseed_report_file Sex inference output file, SNP6 pipeline only (Default: birdseed.report.txt)
#' @param heterozygousFilter Legacy option to set a heterozygous SNP filter, SNP6 pipeline only (Default: "none")
#' @param prior_breakpoints_file A two column file with prior breakpoints to be used during segmentation (Default: NULL)
#' @author sd11
#' @export
battenberg = function(tumourname, normalname, tumour_data_file, normal_data_file, imputeinfofile, g1000prefix, problemloci, gccorrectprefix=NULL,
                      repliccorrectprefix=NULL, g1000allelesprefix=NA, ismale=NA, data_type="wgs", impute_exe="impute2", allelecounter_exe="alleleCounter", nthreads=8, platform_gamma=1, phasing_gamma=1,
                      segmentation_gamma=10, segmentation_kmin=3, phasing_kmin=1, clonality_dist_metric=0, ascat_dist_metric=1, min_ploidy=1.6,
                      max_ploidy=4.8, min_rho=0.1, min_goodness=0.63, uninformative_BAF_threshold=0.51, min_normal_depth=10, min_base_qual=20, 
                      min_map_qual=35, calc_seg_baf_option=3, skip_allele_counting=F, skip_preprocessing=F, skip_phasing=F,
                      snp6_reference_info_file=NA, apt.probeset.genotype.exe="apt-probeset-genotype", apt.probeset.summarize.exe="apt-probeset-summarize", 
                      norm.geno.clust.exe="normalize_affy_geno_cluster.pl", birdseed_report_file="birdseed.report.txt", heterozygousFilter="none",
                      prior_breakpoints_file=NULL) {
  
  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")
  
  if (data_type=="wgs" & is.na(ismale)) {
    stop("Please provide a boolean denominator whether this sample represents a male donor")
  }
  
  if (data_type=="wgs" & is.na(g1000allelesprefix)) {
    stop("Please provide a path to 1000 Genomes allele reference files")
  }
  
  if (data_type=="wgs" & is.null(gccorrectprefix)) {
    stop("Please provide a path to GC content reference files")
  }

  if (!file.exists(problemloci)) {
       stop("Please provide a path to a problematic loci file")
  }

  if (!file.exists(imputeinfofile)) {
	  stop("Please provide a path to an impute info file")
  }

  # check whether the impute_info.txt file contains correct paths
  check.imputeinfofile(imputeinfofile, ismale)

  if (data_type=="wgs" | data_type=="WGS") {
    chrom_names = get.chrom.names(imputeinfofile, ismale)
    logr_file = paste(tumourname, "_mutantLogR_gcCorrected.tab", sep="")
    allelecounts_file = paste(tumourname, "_alleleCounts.tab", sep="")
  } else if (data_type=="snp6" | data_type=="SNP6") {
    chrom_names = get.chrom.names(imputeinfofile, TRUE)
    logr_file = paste(tumourname, "_mutantLogR.tab", sep="")
    allelecounts_file = NULL
  }
 
  if (!skip_preprocessing) {
    if (data_type=="wgs" | data_type=="WGS") {
      # Setup for parallel computing
      clp = parallel::makeCluster(nthreads)
      doParallel::registerDoParallel(clp)
            
      prepare_wgs(chrom_names=chrom_names, 
                  tumourbam=tumour_data_file, 
                  normalbam=normal_data_file, 
                  tumourname=tumourname, 
                  normalname=normalname, 
                  g1000allelesprefix=g1000allelesprefix, 
                  g1000prefix=g1000prefix, 
                  gccorrectprefix=gccorrectprefix,
                  repliccorrectprefix=repliccorrectprefix,
                  min_base_qual=min_base_qual, 
                  min_map_qual=min_map_qual, 
                  allelecounter_exe=allelecounter_exe, 
                  min_normal_depth=min_normal_depth, 
                  nthreads=nthreads,
                  skip_allele_counting=skip_allele_counting)
      
      # Kill the threads
      parallel::stopCluster(clp)
      
    } else if (data_type=="snp6" | data_type=="SNP6") {
      
      prepare_snp6(tumour_cel_file=tumour_data_file, 
                   normal_cel_file=normal_data_file, 
                   tumourname=tumourname, 
                   chrom_names=chrom_names, 
                   snp6_reference_info_file=snp6_reference_info_file, 
                   apt.probeset.genotype.exe=apt.probeset.genotype.exe, 
                   apt.probeset.summarize.exe=apt.probeset.summarize.exe,
                   norm.geno.clust.exe=norm.geno.clust.exe, 
                   birdseed_report_file=birdseed_report_file)
      
    } else {
      print("Unknown data type provided, please provide wgs or snp6")
      q(save="no", status=1)
    }
  }
  
  if (data_type=="snp6" | data_type=="SNP6") {
    # Infer what the gender is - WGS requires it to be specified
    gender = infer_gender_birdseed(birdseed_report_file)
    ismale = gender == "male"
  }
 
  if (!skip_phasing) {
    # Setup for parallel computing
    clp = parallel::makeCluster(nthreads)
    doParallel::registerDoParallel(clp)
    
    # Reconstruct haplotypes 
    # mclapply(1:length(chrom_names), function(chrom) {
    foreach::foreach (chrom=1:length(chrom_names)) %dopar% {
      print(chrom)
      
      run_haplotyping(chrom=chrom, 
                      tumourname=tumourname, 
                      normalname=normalname, 
                      ismale=ismale, 
                      imputeinfofile=imputeinfofile, 
                      problemloci=problemloci, 
                      impute_exe=impute_exe, 
                      min_normal_depth=min_normal_depth,
  		                chrom_names=chrom_names,
  		                snp6_reference_info_file=snp6_reference_info_file,
  		                heterozygousFilter=heterozygousFilter)
    }#, mc.cores=nthreads)
    
    # Kill the threads as from here its all single core
    parallel::stopCluster(clp)
    
    # Combine all the BAF output into a single file
    combine.baf.files(inputfile.prefix=paste(tumourname, "_chr", sep=""), 
                      inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt", 
                      outputfile=paste(tumourname, "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                      no.chrs=length(chrom_names))
  }
  
  # Segment the phased and haplotyped BAF data
  segment.baf.phased(samplename=tumourname,
                     inputfile=paste(tumourname, "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
                     outputfile=paste(tumourname, ".BAFsegmented.txt", sep=""),
                     prior_breakpoints_file=prior_breakpoints_file,
                     gamma=segmentation_gamma,
                     phasegamma=phasing_gamma,
                     kmin=segmentation_kmin,
                     phasekmin=phasing_kmin,
                     calc_seg_baf_option=calc_seg_baf_option)
  
  # Fit a clonal copy number profile
  fit.copy.number(samplename=tumourname,
                  outputfile.prefix=paste(tumourname, "_", sep=""),
                  inputfile.baf.segmented=paste(tumourname, ".BAFsegmented.txt", sep=""), 
                  inputfile.baf=paste(tumourname,"_mutantBAF.tab", sep=""), 
                  inputfile.logr=logr_file, 
                  dist_choice=clonality_dist_metric, 
                  ascat_dist_choice=ascat_dist_metric, 
                  min.ploidy=min_ploidy, 
                  max.ploidy=max_ploidy, 
                  min.rho=min_rho, 
                  min.goodness=min_goodness, 
                  uninformative_BAF_threshold=uninformative_BAF_threshold, 
                  gamma_param=platform_gamma, 
                  use_preset_rho_psi=F, 
                  preset_rho=NA, 
                  preset_psi=NA, 
                  read_depth=30)
  
  # Go over all segments, determine which segements are a mixture of two states and fit a second CN state
  callSubclones(sample.name=tumourname, 
                baf.segmented.file=paste(tumourname, ".BAFsegmented.txt", sep=""), 
                logr.file=logr_file, 
                rho.psi.file=paste(tumourname, "_rho_and_psi.txt",sep=""), 
                output.file=paste(tumourname,"_subclones.txt", sep=""), 
                output.figures.prefix=paste(tumourname,"_subclones_chr", sep=""), 
                output.gw.figures.prefix=paste(tumourname,"_BattenbergProfile", sep=""),
                masking_output_file=paste(tumourname, "_segment_masking_details.txt", sep=""),
                prior_breakpoints_file=prior_breakpoints_file,
                chr_names=chrom_names, 
                gamma=platform_gamma, 
                segmentation.gamma=NA, 
                siglevel=0.05, 
                maxdist=0.01, 
                noperms=1000,
                calc_seg_baf_option=calc_seg_baf_option)
  
  # Make some post-hoc plots
  make_posthoc_plots(samplename=tumourname, 
                     logr_file=logr_file, 
                     subclones_file=paste(tumourname, "_subclones.txt", sep=""), 
                     rho_psi_file=paste(tumourname, "_rho_and_psi.txt", sep=""), 
                     bafsegmented_file=paste(tumourname, ".BAFsegmented.txt", sep=""), 
                     logrsegmented_file=paste(tumourname, ".logRsegmented.txt", sep=""), 
                     allelecounts_file=allelecounts_file)
  
  # Save refit suggestions for a future rerun
  cnfit_to_refit_suggestions(samplename=tumourname, 
                             subclones_file=paste(tumourname, "_subclones.txt", sep=""),
                             rho_psi_file=paste(tumourname, "_rho_and_psi.txt", sep=""),
                             gamma_param=platform_gamma)
  
}

