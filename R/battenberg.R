
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
#' @param segmentation_gamma_multisample The gamma parameter controls the size of the penalty of starting a new segment during mutlisample segmentation. It is the key parameter for controlling the number of segments (Default: 10)
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
#' @param usebeagle Should use beagle5 instead of impute2 Default: FALSE
#' @param beaglejar Full path to Beagle java jar file Default: NA
#' @param beagleref.template Full path template to Beagle reference files where the chromosome is replaced by 'CHROMNAME' Default: NA
#' @param beagleplink.template Full path template to Beagle plink files where the chromosome is replaced by 'CHROMNAME' Default: NA
#' @param beaglemaxmem Integer Beagle max heap size in Gb  Default: 10
#' @param beaglenthreads Integer number of threads used by beagle5 Default:1
#' @param beaglewindow Integer size of the genomic window for beagle5 (cM) Default:40
#' @param beagleoverlap Integer size of the overlap between windows beagle5 Default:4
#' @param javajre Path to the Java JRE executable, only required for haplotype reconstruction with Beagle (default java, i.e. in $PATH)
#' @param snp6_reference_info_file Reference files for the SNP6 pipeline only (Default: NA)
#' @param apt.probeset.genotype.exe Helper tool for extracting data from CEL files, SNP6 pipeline only (Default: apt-probeset-genotype)
#' @param apt.probeset.summarize.exe  Helper tool for extracting data from CEL files, SNP6 pipeline only (Default: apt-probeset-summarize)
#' @param norm.geno.clust.exe  Helper tool for extracting data from CEL files, SNP6 pipeline only (Default: normalize_affy_geno_cluster.pl)
#' @param birdseed_report_file Sex inference output file, SNP6 pipeline only (Default: birdseed.report.txt)
#' @param heterozygousFilter Legacy option to set a heterozygous SNP filter, SNP6 pipeline only (Default: "none")
#' @param prior_breakpoints_file A two column file with prior breakpoints to be used during segmentation (Default: NULL)
#' @param GENOMEBUILD Genome build upon which the 1000G SNP coordinates were obtained (Default: hg19; options: "hg19" or "hg38")  
#' @param externalhaplotypefile Vcf containing externally obtained haplotype blocks (Default: NA)
#' @param write_battenberg_phasing Write the Battenberg phasing results as vcf to disk, e.g. for multisample cases (Default: TRUE)
#' @param multisample_maxlag Maximal number of upstream SNPs used in the multisample haplotyping to inform the haplotype at another SNP (Default: 100)
#' @param multisample_relative_weight_balanced Relative weight to give to haplotype info from a sample without allelic imbalance in the region (Default: 0.25)
#' @author sd11, jdemeul, Naser Ansari-Pour
#' @export
battenberg = function(tumourname, normalname, tumour_data_file, normal_data_file, imputeinfofile, g1000prefix, problemloci, gccorrectprefix=NULL,
                      repliccorrectprefix=NULL, g1000allelesprefix=NA, ismale=NA, data_type="wgs", impute_exe="impute2", allelecounter_exe="alleleCounter", nthreads=8, platform_gamma=1, phasing_gamma=1,
                      segmentation_gamma=10, segmentation_kmin=3, phasing_kmin=1, clonality_dist_metric=0, ascat_dist_metric=1, min_ploidy=1.6,
                      max_ploidy=4.8, min_rho=0.1, min_goodness=0.63, uninformative_BAF_threshold=0.51, min_normal_depth=10, min_base_qual=20,
                      min_map_qual=35, calc_seg_baf_option=3, skip_allele_counting=F, skip_preprocessing=F, skip_phasing=F, externalhaplotypefile = NA,
                      usebeagle=FALSE,
                      beaglejar=NA,
                      beagleref.template=NA,
                      beagleplink.template=NA,
                      beaglemaxmem=10,
                      beaglenthreads=1,
                      beaglewindow=40,
                      beagleoverlap=4,
                      javajre="java",
                      write_battenberg_phasing = T, multisample_relative_weight_balanced = 0.25, multisample_maxlag = 100, segmentation_gamma_multisample = 5,
                      snp6_reference_info_file=NA, apt.probeset.genotype.exe="apt-probeset-genotype", apt.probeset.summarize.exe="apt-probeset-summarize",
                      norm.geno.clust.exe="normalize_affy_geno_cluster.pl", birdseed_report_file="birdseed.report.txt", heterozygousFilter="none",
                      prior_breakpoints_file=NULL, GENOMEBUILD="hg19") {
  
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
  check.imputeinfofile(imputeinfofile = imputeinfofile, is.male = ismale, usebeagle = usebeagle)
  
  # check whether multisample case
  nsamples <- length(tumourname)
  if (nsamples > 1) {
    if (length(skip_allele_counting) < nsamples) {
      skip_allele_counting = rep(skip_allele_counting[1], nsamples)
    }
    if (length(skip_preprocessing) < nsamples) {
      skip_preprocessing = rep(skip_preprocessing[1], nsamples)
    }
    if (length(skip_phasing) < nsamples) {
      skip_phasing = rep(skip_phasing[1], nsamples)
    }
  }
  
  
  if (data_type=="wgs" | data_type=="WGS") {
    if (nsamples > 1) {
      print(paste0("Running Battenberg in multisample mode on ", nsamples, " samples: ", paste0(tumourname, collapse = ", ")))
    }
    chrom_names = get.chrom.names(imputeinfofile, ismale)
  } else if (data_type=="snp6" | data_type=="SNP6") {
    if (nsamples > 1) {
      stop(paste0("Battenberg multisample mode has not been tested with SNP6 data"))
    }
    chrom_names = get.chrom.names(imputeinfofile, TRUE)
    logr_file = paste(tumourname, "_mutantLogR.tab", sep="")
    allelecounts_file = NULL
  }
  
  for (sampleidx in 1:nsamples) {
    
    if (!skip_preprocessing[sampleidx]) {
      if (data_type=="wgs" | data_type=="WGS") {
        # Setup for parallel computing
        clp = parallel::makeCluster(nthreads)
        doParallel::registerDoParallel(clp)
        
        prepare_wgs(chrom_names=chrom_names,
                    tumourbam=tumour_data_file[sampleidx],
                    normalbam=normal_data_file,
                    tumourname=tumourname[sampleidx],
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
                    skip_allele_counting=skip_allele_counting[sampleidx],
                    skip_allele_counting_normal = (sampleidx > 1))
        
        # Kill the threads
        parallel::stopCluster(clp)
        
      } else if (data_type=="snp6" | data_type=="SNP6") {
        
        prepare_snp6(tumour_cel_file=tumour_data_file[sampleidx],
                     normal_cel_file=normal_data_file,
                     tumourname=tumourname[sampleidx],
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
    
    
    if (!skip_phasing[sampleidx]) {
      
      # if external phasing data is provided (as a vcf), split into chromosomes for use in haplotype reconstruction
      if (!is.na(externalhaplotypefile) && file.exists(externalhaplotypefile)) {
        externalhaplotypeprefix <- paste0(normalname, "_external_haplotypes_chr")
        
        # if these files exist already, no need to split again
        if (any(!file.exists(paste0(externalhaplotypeprefix, 1:length(chrom_names), ".vcf")))) {
          
          print(paste0("Splitting external phasing data from ", externalhaplotypefile))
          split_input_haplotypes(chrom_names = chrom_names,
                                 externalhaplotypefile = externalhaplotypefile,
                                 outprefix = externalhaplotypeprefix)
        } else {
          print("No need to split, external haplotype files per chromosome found")
        }
      } else {
        externalhaplotypeprefix <- NA
      }
      
      # Setup for parallel computing
      clp = parallel::makeCluster(nthreads)
      doParallel::registerDoParallel(clp)
      
      # Reconstruct haplotypes
      # mclapply(1:length(chrom_names), function(chrom) {
      foreach::foreach (chrom=1:length(chrom_names)) %dopar% {
        print(chrom)
        
        run_haplotyping(chrom=chrom,
                        tumourname=tumourname[sampleidx],
                        normalname=normalname,
                        ismale=ismale,
                        imputeinfofile=imputeinfofile,
                        problemloci=problemloci,
                        impute_exe=impute_exe,
                        min_normal_depth=min_normal_depth,
                        chrom_names=chrom_names,
                        snp6_reference_info_file=snp6_reference_info_file,
                        heterozygousFilter=heterozygousFilter,
                        usebeagle=usebeagle,
                        beaglejar=beaglejar,
                        beagleref=gsub("CHROMNAME",if(chrom==23) "X" else chrom, beagleref.template),
                        beagleplink=gsub("CHROMNAME",if(chrom==23) "X" else chrom, beagleplink.template),
                        beaglemaxmem=beaglemaxmem,
                        beaglenthreads=beaglenthreads,
                        beaglewindow=beaglewindow,
                        beagleoverlap=beagleoverlap,
                        externalhaplotypeprefix=externalhaplotypeprefix,
                        use_previous_imputation=(sampleidx > 1))
      }#, mc.cores=nthreads)
      
      # Kill the threads as from here its all single core
      parallel::stopCluster(clp)
      
      # Combine all the BAF output into a single file
      combine.baf.files(inputfile.prefix=paste(tumourname[sampleidx], "_chr", sep=""),
                        inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt",
                        outputfile=paste(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                        no.chrs=length(chrom_names))
    }
    
    # Segment the phased and haplotyped BAF data
    segment.baf.phased(samplename=tumourname[sampleidx],
                       inputfile=paste(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
                       outputfile=paste(tumourname[sampleidx], ".BAFsegmented.txt", sep=""),
                       prior_breakpoints_file=prior_breakpoints_file,
                       gamma=segmentation_gamma,
                       phasegamma=phasing_gamma,
                       kmin=segmentation_kmin,
                       phasekmin=phasing_kmin,
                       calc_seg_baf_option=calc_seg_baf_option)
    
    if (nsamples > 1 | write_battenberg_phasing) {
      # Write the Battenberg phasing information to disk as a vcf
      write_battenberg_phasing(tumourname = tumourname[sampleidx],
                               SNPfiles = paste0(tumourname[sampleidx], "_alleleFrequencies_chr", 1:length(chrom_names), ".txt"),
                               imputedHaplotypeFiles = paste0(tumourname[sampleidx], "_impute_output_chr", 1:length(chrom_names), "_allHaplotypeInfo.txt"),
                               bafsegmented_file = paste0(tumourname[sampleidx], ".BAFsegmented.txt"),
                               outprefix = paste0(tumourname[sampleidx], "_Battenberg_phased_chr"),
                               chrom_names = chrom_names,
                               include_homozygous = F)
    }
    
  }
  
  
  # if this is a multisample run, combine the battenberg phasing outputs, incorporate it and resegment
  if (nsamples > 1) {
    print("Constructing multisample phasing")
    multisamplehaplotypeprefix <- paste0(normalname, "_multisample_haplotypes_chr")

    
    # Setup for parallel computing
    clp = parallel::makeCluster(nthreads)
    doParallel::registerDoParallel(clp)
    
    # Reconstruct haplotypes
    # mclapply(1:length(chrom_names), function(chrom) {
    foreach::foreach (chrom=1:length(chrom_names)) %dopar% {
      print(chrom)
      
      get_multisample_phasing(chrom = chrom,
                              bbphasingprefixes = paste0(tumourname, "_Battenberg_phased_chr"),
                              maxlag = multisample_maxlag,
                              relative_weight_balanced = multisample_relative_weight_balanced,
                              outprefix = multisamplehaplotypeprefix)
      
    }

    
    # continue over all samples to incorporate the multisample phasing
    for (sampleidx in 1:nsamples) {
      
      # rename the original files without multisample phasing info
      MutBAFfiles <- paste0(tumourname[sampleidx], "_chr", 1:length(chrom_names), "_heterozygousMutBAFs_haplotyped.txt")
      heterozygousdatafiles <- paste0(tumourname[sampleidx], "_chr", 1:length(chrom_names), "_heterozygousData.png")
      raffiles <- paste0(tumourname[sampleidx], "_RAFseg_chr", chrom_names, ".png")
      segfiles <- paste0(tumourname[sampleidx], "_segment_chr", chrom_names, ".png")
      haplotypedandbafsegmentedfiles <- paste0(tumourname[sampleidx], c("_heterozygousMutBAFs_haplotyped.txt", ".BAFsegmented.txt"))
      
      file.copy(from = MutBAFfiles, to = gsub(pattern = ".txt$", replacement = "_noMulti.txt", x = MutBAFfiles), overwrite = T)
      file.copy(from = heterozygousdatafiles, to = gsub(pattern = ".png$", replacement = "_noMulti.png", x = heterozygousdatafiles), overwrite = T)
      file.copy(from = raffiles, to = gsub(pattern = ".png$", replacement = "_noMulti.png", x = raffiles), overwrite = T)
      file.copy(from = segfiles, to = gsub(pattern = ".png$", replacement = "_noMulti.png", x = segfiles), overwrite = T)
      file.copy(from = haplotypedandbafsegmentedfiles, to = gsub(pattern = ".txt$", replacement = "_noMulti.txt", x = haplotypedandbafsegmentedfiles), overwrite = T)
      # done renaming, next sections will overwrite orignals
      
      
      foreach::foreach (chrom=1:length(chrom_names)) %dopar% {
        print(chrom)
        
        input_known_haplotypes(chrom = chrom,
                               chrom_names = chrom_names,
                               imputedHaplotypeFile = paste0(tumourname[sampleidx], "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"),
                               externalHaplotypeFile = paste0(multisamplehaplotypeprefix, chrom, ".vcf"),
                               oldfilesuffix = "_noMulti.txt")
        
        GetChromosomeBAFs(chrom=chrom,
                          SNP_file=paste0(tumourname[sampleidx], "_alleleFrequencies_chr", chrom, ".txt"),
                          haplotypeFile=paste0(tumourname[sampleidx], "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"),
                          samplename=tumourname[sampleidx],
                          outfile=paste0(tumourname[sampleidx], "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"),
                          chr_names=chrom_names,
                          minCounts=min_normal_depth)
        
        # Plot what we have until this point
        plot.haplotype.data(haplotyped.baf.file=paste0(tumourname[sampleidx], "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"),
                            imageFileName=paste0(tumourname[sampleidx],"_chr",chrom,"_heterozygousData.png"),
                            samplename=tumourname[sampleidx],
                            chrom=chrom,
                            chr_names=chrom_names)
      }
      
    }
    
    # Kill the threads as from here its single core
    parallel::stopCluster(clp)
    
    for (sampleidx in 1:nsamples) {
      
      # Combine all the BAF output into a single file
      combine.baf.files(inputfile.prefix=paste0(tumourname[sampleidx], "_chr"),
                        inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt",
                        outputfile=paste0(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt"), 
                        no.chrs=length(chrom_names))
      
    }
    
    # Segment the phased and haplotyped BAF data
    segment.baf.phased.multisample(samplename=tumourname,
                                   inputfile=paste(tumourname, "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
                                   outputfile=paste(tumourname, ".BAFsegmented.txt", sep=""),
                                   prior_breakpoints_file=prior_breakpoints_file,
                                   gamma=segmentation_gamma_multisample,
                                   calc_seg_baf_option=calc_seg_baf_option)
    
  }
  
  
  # Setup for parallel computing
  clp = parallel::makeCluster(min(nthreads, nsamples))
  doParallel::registerDoParallel(clp)
  
  # for (sampleidx in 1:nsamples) {
  foreach::foreach (sampleidx=1:nsamples) %dopar% {
    print(paste0("Fitting final copy number and calling subclones for sample ", tumourname[sampleidx]))
    
    if (data_type=="wgs" | data_type=="WGS") {
      logr_file = paste(tumourname[sampleidx], "_mutantLogR_gcCorrected.tab", sep="")
      allelecounts_file = paste(tumourname[sampleidx], "_alleleCounts.tab", sep="")
    }
    
    # Fit a clonal copy number profile
    fit.copy.number(samplename=tumourname[sampleidx],
                    outputfile.prefix=paste(tumourname[sampleidx], "_", sep=""),
                    inputfile.baf.segmented=paste(tumourname[sampleidx], ".BAFsegmented.txt", sep=""),
                    inputfile.baf=paste(tumourname[sampleidx],"_mutantBAF.tab", sep=""),
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
    callSubclones(sample.name=tumourname[sampleidx],
                  baf.segmented.file=paste(tumourname[sampleidx], ".BAFsegmented.txt", sep=""),
                  logr.file=logr_file,
                  rho.psi.file=paste(tumourname[sampleidx], "_rho_and_psi.txt",sep=""),
                  output.file=paste(tumourname[sampleidx],"_subclones.txt", sep=""),
                  output.figures.prefix=paste(tumourname[sampleidx],"_subclones_chr", sep=""),
                  output.gw.figures.prefix=paste(tumourname[sampleidx],"_BattenbergProfile", sep=""),
                  masking_output_file=paste(tumourname[sampleidx], "_segment_masking_details.txt", sep=""),
                  prior_breakpoints_file=prior_breakpoints_file,
                  chr_names=chrom_names, 
                  gamma=platform_gamma, 
                  segmentation.gamma=NA, 
                  siglevel=0.05, 
                  maxdist=0.01, 
                  noperms=1000,
                  calc_seg_baf_option=calc_seg_baf_option)
    
    # If patient is male, get copy number status of ChrX based only on logR segmentation (due to hemizygosity of SNPs)
    if (ismale){
      callChrXsubclones(TUMOURNAME=tumourname[sampleidx],
                        X_GAMMA=1000,
                        X_KMIN=100,
                        GENOMEBUILD=GENOMEBUILD,
                        AR=TRUE)
    }
    
    # Make some post-hoc plots
    make_posthoc_plots(samplename=tumourname[sampleidx],
                       logr_file=logr_file,
                       subclones_file=paste(tumourname[sampleidx], "_subclones.txt", sep=""),
                       rho_psi_file=paste(tumourname[sampleidx], "_rho_and_psi.txt", sep=""),
                       bafsegmented_file=paste(tumourname[sampleidx], ".BAFsegmented.txt", sep=""),
                       logrsegmented_file=paste(tumourname[sampleidx], ".logRsegmented.txt", sep=""),
                       allelecounts_file=allelecounts_file)
    
    # Save refit suggestions for a future rerun
    cnfit_to_refit_suggestions(samplename=tumourname[sampleidx],
                               subclones_file=paste(tumourname[sampleidx], "_subclones.txt", sep=""),
                               rho_psi_file=paste(tumourname[sampleidx], "_rho_and_psi.txt", sep=""),
                               gamma_param=platform_gamma)
    
  }
  
  # Kill the threads as last part again is single core
  parallel::stopCluster(clp)
  
  if (nsamples > 1) {
    print("Assessing mirrored subclonal allelic imbalance (MSAI)")
    call_multisample_MSAI(rdsprefix = multisamplehaplotypeprefix,
                          subclonesfiles = paste0(tumourname, "_subclones.txt"),
                          chrom_names = chrom_names,
                          tumournames = tumourname,
                          plotting = T)
  }
}

