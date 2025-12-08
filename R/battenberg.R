
#' Run the Battenberg pipeline
#'
#' @param analysis The mode of Battenberg copy number analysis to be undertaken: 'paired' for tumour-normal pair, 'cell_line' for Cell line tumour-only and 'germline' for germline CNV of normal sample (Default: 'paired')
#' @param samplename Sample identifier (tumour or germline), this is used as a prefix for the output files. If allele counts are supplied separately, they are expected to have this identifier as prefix.
#' @param normalname Matched normal identifier, this is used as a prefix for the output files. If allele counts are supplied separately, they are expected to have this identifier as prefix.
#' @param sample_data_file A BAM or CEL file for the sample
#' @param normal_data_file A BAM or CEL file for the normal-pair (paired analysis)
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
#' @param max_rho Maximum purity to be considered (Default: 1.0)
#' @param min_goodness Minimum goodness of fit required for a purity/ploidy combination to be accepted as a solution (Default: 0.63)
#' @param uninformative_BAF_threshold The threshold beyond which BAF becomes uninformative (Default: 0.51)
#' @param min_normal_depth Minimum depth required in the matched normal for a SNP to be considered as part of the wgs analysis (Default: 10)
#' @param min_base_qual Minimum base quality required for a read to be counted when allele counting (Default: 20)
#' @param min_map_qual Minimum mapping quality required for a read to be counted when allele counting (Default: 35)
#' @param max_allowed_state The maximum CN state allowed (Default 250)
#' @param cn_upper_limit Maximum number of copy number that can be called (Default 1000)
#' @param calc_seg_baf_option Sets way to calculate BAF per segment: 1=mean, 2=median, 3=ifelse median==0 | 1, mean, median (Default (paired): 3, cell_line & germline: 1)
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
#' @param genomebuild Genome build upon which the 1000G SNP coordinates were obtained (Default: hg19; options: "hg19" or "hg38")  
#' @param externalhaplotypefile Vcf containing externally obtained haplotype blocks (Default: NA)
#' @param write_battenberg_phasing Write the Battenberg phasing results as vcf to disk, e.g. for multisample cases (Default: TRUE)
#' @param multisample_maxlag Maximal number of upstream SNPs used in the multisample haplotyping to inform the haplotype at another SNP (Default: 100)
#' @param multisample_relative_weight_balanced Relative weight to give to haplotype info from a sample without allelic imbalance in the region (Default: 0.25)
#' @param enhanced_grid_search Should use multi-start, parallelized and multi-approach grid search (Default: FALSE)
#' @author sd11, jdemeul, Naser Ansari-Pour, Julio Cesar Cortes Rios
#' @export
battenberg = function(analysis="paired",
                      samplename,
                      normalname,
                      sample_data_file,
                      normal_data_file,
                      imputeinfofile,
                      g1000prefix,
                      problemloci,
                      gccorrectprefix=NULL,
                      repliccorrectprefix=NULL,
                      g1000allelesprefix=NA,
                      ismale=NA,
                      data_type="wgs",
                      impute_exe="impute2",
                      allelecounter_exe="alleleCounter",
                      nthreads=8,
                      platform_gamma=1,
                      phasing_gamma=1,
                      segmentation_gamma=10,
                      segmentation_kmin=3,
                      phasing_kmin=1,
                      clonality_dist_metric=0,
                      ascat_dist_metric=1,
                      min_ploidy=1.6,
                      max_ploidy=4.8,
                      min_rho=0.1,
                      max_rho=1.0,
                      min_goodness=0.63,
                      uninformative_BAF_threshold=0.51,
                      min_normal_depth=10,
                      min_base_qual=20,
                      min_map_qual=35,
                      max_allowed_state=250,
                      cn_upper_limit=1000,
                      calc_seg_baf_option=3,
                      skip_allele_counting=F,
                      skip_preprocessing=F,
                      skip_phasing=F,
                      externalhaplotypefile = NA,
                      usebeagle=FALSE,
                      beaglejar=NA,
                      beagleref.template=NA,
                      beagleplink.template=NA,
                      beaglemaxmem=10,
                      beaglenthreads=1,
                      beaglewindow=40,
                      beagleoverlap=4,
                      javajre="java",
                      write_battenberg_phasing = T,
                      multisample_relative_weight_balanced = 0.25,
                      multisample_maxlag = 90,
                      segmentation_gamma_multisample = 5,
                      snp6_reference_info_file=NA,
                      apt.probeset.genotype.exe="apt-probeset-genotype",
                      apt.probeset.summarize.exe="apt-probeset-summarize",
                      norm.geno.clust.exe="normalize_affy_geno_cluster.pl",
                      birdseed_report_file="birdseed.report.txt",
                      heterozygousFilter="none",
                      prior_breakpoints_file=NULL,
                      genomebuild="hg19",
                      chrom_coord_file=NULL,
		      enhanced_grid_search = F) {

  requireNamespace("foreach")
  requireNamespace("doParallel")
  requireNamespace("parallel")
  libs <- .libPaths()
 
  if (analysis == "cell_line"){
    calc_seg_baf_option=1
    phasing_gamma=1
    phasing_kmin=2
    segmentation_gamma=20
    segmentation_kmin=3
    # no matched normal required, but we  are generating normal counts which have this name coded
    normalname = paste0(samplename, "_normal")
    # other cell_line specific parameter values
    min_ploidy=min_ploidy
    max_ploidy=max_ploidy
    min_rho=0.99
    max_rho=1.01
  }
  if (analysis == "germline"){
    calc_seg_baf_option=1
    phasing_gamma=3
    phasing_kmin=1
    segmentation_gamma=3
    segmentation_kmin=3
    # no matched normal required, but we  are generating normal counts which have this name coded
    normalname = paste0(samplename, "_normal")
    min_ploidy=1.5
    max_ploidy=2.5
    min_rho=0.99
    max_rho=1.01
  }
  
  if (data_type=="wgs" & is.na(ismale)) {
    stop("Please provide a boolean denominator whether this sample represents a male donor")
  }
  
  if (data_type=="wgs" & is.na(g1000allelesprefix)) {
    stop("Please provide a path to 1000 Genomes allele reference files")
  }
  
  if (data_type=="wgs" & is.null(gccorrectprefix)) {
    stop("Please provide a path to GC content reference files")
  }
  
  if (data_type=="wgs" && !file.exists(problemloci)) {
    stop("Please provide a path to a problematic loci file")
  }
  
  if (!file.exists(imputeinfofile)) {
    stop("Please provide a path to an impute info file")
  }
  
  # check whether the impute_info.txt file contains correct paths
  check.imputeinfofile(imputeinfofile = imputeinfofile, is.male = ismale, usebeagle = usebeagle)
  
  # check whether multisample case
  nsamples <- length(samplename)
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
      print(paste0("Running Battenberg in multisample mode on ", nsamples, " samples: ", paste0(samplename, collapse = ", ")))
    }
    chrom_names = get.chrom.names(imputeinfofile, ismale, analysis=analysis)
  } else if (data_type=="snp6" | data_type=="SNP6") {
    if (nsamples > 1) {
      stop(paste0("Battenberg multisample mode has not been tested with SNP6 data"))
    }
    chrom_names = get.chrom.names(imputeinfofile, TRUE)
    logr_file = paste(samplename, "_mutantLogR.tab", sep="")
    allelecounts_file = NULL
  }
  print(chrom_names) 
  for (sampleidx in 1:nsamples) {
    if (!skip_preprocessing[sampleidx]) {
      if (data_type=="wgs" | data_type=="WGS") {
        # Setup for parallel computing
        clp = parallel::makeCluster(nthreads,outfile="")
        doParallel::registerDoParallel(clp)
        
        if (analysis == "paired"){
          
          if (is.null(normalname)|is.na(normalname)){
            stop("No normal sample is specified for 'paired analysis' - a normal paired BAM is required")
            }
          prepare_wgs(chrom_names=chrom_names,
                      tumourbam=sample_data_file[sampleidx],
                      normalbam=normal_data_file,
                      tumourname=samplename[sampleidx],
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
          
        } else if (analysis == "cell_line") {
          prepare_wgs_cell_line(chrom_names=chrom_names,
                                chrom_coord=chrom_coord_file,
                                tumourbam=sample_data_file,
                                tumourname=samplename,
                                g1000lociprefix=g1000prefix,
                                g1000allelesprefix=g1000allelesprefix, 
                                gamma_ivd=1e5,
                                kmin_ivd=50,
                                centromere_noise_seg_size=1e6,
                                centromere_dist=5e5,
                                min_het_dist=1e5, 
                                gamma_logr=100,
                                length_adjacent=5e4,
                                gccorrectprefix=gccorrectprefix, 
                                repliccorrectprefix=repliccorrectprefix,
                                min_base_qual=min_base_qual,
                                min_map_qual=min_map_qual, 
                                allelecounter_exe=allelecounter_exe,
                                min_normal_depth=min_normal_depth,
                                skip_allele_counting=skip_allele_counting[sampleidx])
        } else if (analysis == "germline"){
          
          prepare_wgs_germline(chrom_names=chrom_names,
                               chrom_coord=chrom_coord_file,
                               germlinebam=sample_data_file,
                               germlinename=samplename,
                               g1000lociprefix=g1000prefix,
                               g1000allelesprefix=g1000allelesprefix,
                               gamma_ivd=1e5,
                               kmin_ivd=50,
                               centromere_noise_seg_size=1e6,
                               centromere_dist=5e5,
                               min_het_dist=2e3,
                               gamma_logr=100,
                               length_adjacent=5e4,
                               gccorrectprefix=gccorrectprefix,
                               repliccorrectprefix=repliccorrectprefix,
                               min_base_qual=min_base_qual,
                               min_map_qual=min_map_qual,
                               allelecounter_exe=allelecounter_exe,
                               min_normal_depth=min_normal_depth,
                               skip_allele_counting=skip_allele_counting[sampleidx])
        }
        
        
        # Kill the threads
        parallel::stopCluster(clp)
        
      } else if (data_type=="snp6" | data_type=="SNP6") {
        
        prepare_snp6(tumour_cel_file=sample_data_file[sampleidx],
                     normal_cel_file=normal_data_file,
                     tumourname=samplename[sampleidx],
                     chrom_names=chrom_names,
                     snp6_reference_info_file=snp6_reference_info_file,
                     apt.probeset.genotype.exe=apt.probeset.genotype.exe,
                     apt.probeset.summarize.exe=apt.probeset.summarize.exe,
                     norm.geno.clust.exe=norm.geno.clust.exe,
                     birdseed_report_file=birdseed_report_file,
                     genomebuild=genomebuild)
        
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
      clp = parallel::makeCluster(nthreads,outfile="")
      doParallel::registerDoParallel(clp)
      
      # Reconstruct haplotypes
      # mclapply(1:length(chrom_names), function(chrom) {
      if (analysis=="germline"){
        foreach::foreach (i=1:length(chrom_names)) %dopar% {
          .libPaths(libs)
          chrom = chrom_names[i]
          print(chrom)
          
          run_haplotyping_germline(chrom=chrom,
                                   germlinename=samplename,
                                   normalname=normalname,
                                   ismale=ismale,
                                   imputeinfofile=imputeinfofile,
                                   problemloci=problemloci,
                                   impute_exe=impute_exe,
                                   min_normal_depth=min_normal_depth,
                                   chrom_names=chrom_names, 
                                   externalhaplotypeprefix = NA,
                                   use_previous_imputation=F,
                                   snp6_reference_info_file=NA, 
                                   heterozygousFilter=NA,
                                   usebeagle=usebeagle,
                                   beaglejar=beaglejar,
                                   beagleref=gsub("CHROMNAME",chrom,beagleref.template),
                                   beagleplink=gsub("CHROMNAME",chrom,beagleplink.template),
                                   beaglemaxmem=beaglemaxmem,
                                   beaglenthreads=beaglenthreads,
                                   beaglewindow=beaglewindow,
                                   beagleoverlap=beagleoverlap)      
        }
      } else {
        foreach::foreach (i=1:length(chrom_names)) %dopar% {
          .libPaths(libs)
          chrom = chrom_names[i]
          print(chrom)      
          run_haplotyping(chrom=chrom,
                          tumourname=samplename[sampleidx],
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
                          beagleref=gsub("CHROMNAME", chrom, beagleref.template),
                          beagleplink=gsub("CHROMNAME", chrom, beagleplink.template),
                          beaglemaxmem=beaglemaxmem,
                          beaglenthreads=beaglenthreads,
                          beaglewindow=beaglewindow,
                          beagleoverlap=beagleoverlap,
                          externalhaplotypeprefix=externalhaplotypeprefix,
                          use_previous_imputation=(sampleidx > 1))
        }
      }
      
      # Kill the threads as from here its all single core
      parallel::stopCluster(clp)
      
      # Combine all the BAF output into a single file
      combine.baf.files(inputfile.prefix=paste(samplename[sampleidx], "_chr", sep=""),
                        inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt",
                        outputfile=paste(samplename[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                        chr_names=chrom_names)
    }
    
    # Segment the phased and haplotyped BAF data
    segment.baf.phased(samplename=samplename[sampleidx],
                       inputfile=paste(samplename[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
                       outputfile=paste(samplename[sampleidx], ".BAFsegmented.txt", sep=""),
                       prior_breakpoints_file=prior_breakpoints_file,
                       gamma=segmentation_gamma,
                       phasegamma=phasing_gamma,
                       kmin=segmentation_kmin,
                       phasekmin=phasing_kmin,
                       calc_seg_baf_option=calc_seg_baf_option)
    
    if (nsamples > 1 | write_battenberg_phasing) {
      # Write the Battenberg phasing information to disk as a vcf
      write_battenberg_phasing(tumourname = samplename[sampleidx],
                               SNPfiles = paste0(samplename[sampleidx], "_alleleFrequencies_chr", chrom_names, ".txt"),
                               imputedHaplotypeFiles = paste0(samplename[sampleidx], "_impute_output_chr", chrom_names, "_allHaplotypeInfo.txt"),
                               bafsegmented_file = paste0(samplename[sampleidx], ".BAFsegmented.txt"),
                               outprefix = paste0(samplename[sampleidx], "_Battenberg_phased_chr"),
                               chrom_names = chrom_names,
                               include_homozygous = F)
    }
    
  }
  
  # if this is a multisample run, combine the battenberg phasing outputs, incorporate it and resegment
  if (nsamples > 1) {
    print("Constructing multisample phasing")
    multisamplehaplotypeprefix <- paste0(normalname, "_multisample_haplotypes_chr")
    
    
    # Setup for parallel computing
    clp = parallel::makeCluster(nthreads,outfile="")
    doParallel::registerDoParallel(clp)
    
    # Reconstruct haplotypes
    .libPaths()
    foreach::foreach (i=1:length(chrom_names)) %dopar% {
      .libPaths(libs)
    .libPaths()
      chrom = chrom_names[i]
      print(chrom)
      
      get_multisample_phasing(chrom = chrom,
                              bbphasingprefixes = paste0(samplename, "_Battenberg_phased_chr"),
                              maxlag = multisample_maxlag,
                              relative_weight_balanced = multisample_relative_weight_balanced,
                              outprefix = multisamplehaplotypeprefix)
    }
    
    # continue over all samples to incorporate the multisample phasing
    for (sampleidx in 1:nsamples) {
      
      # rename the original files without multisample phasing info
      MutBAFfiles <- paste0(samplename[sampleidx], "_chr", chrom_names, "_heterozygousMutBAFs_haplotyped.txt")
      heterozygousdatafiles <- paste0(samplename[sampleidx], "_chr", chrom_names, "_heterozygousData.png")
      raffiles <- paste0(samplename[sampleidx], "_RAFseg_chr", chrom_names, ".png")
      segfiles <- paste0(samplename[sampleidx], "_segment_chr", chrom_names, ".png")
      haplotypedandbafsegmentedfiles <- paste0(samplename[sampleidx], c("_heterozygousMutBAFs_haplotyped.txt", ".BAFsegmented.txt"))
      
      file.copy(from = MutBAFfiles, to = gsub(pattern = ".txt$", replacement = "_noMulti.txt", x = MutBAFfiles), overwrite = T)
      file.copy(from = heterozygousdatafiles, to = gsub(pattern = ".png$", replacement = "_noMulti.png", x = heterozygousdatafiles), overwrite = T)
      file.copy(from = raffiles, to = gsub(pattern = ".png$", replacement = "_noMulti.png", x = raffiles), overwrite = T)
      file.copy(from = segfiles, to = gsub(pattern = ".png$", replacement = "_noMulti.png", x = segfiles), overwrite = T)
      file.copy(from = haplotypedandbafsegmentedfiles, to = gsub(pattern = ".txt$", replacement = "_noMulti.txt", x = haplotypedandbafsegmentedfiles), overwrite = T)
      # done renaming, next sections will overwrite orignals
      
      
      foreach::foreach (i=1:length(chrom_names)) %dopar% {
        .libPaths(libs)
        chrom = chrom_names[i]
        print(chrom)
        
        input_known_haplotypes(chrom = chrom,
                               chrom_names = chrom_names,
                               imputedHaplotypeFile = paste0(samplename[sampleidx], "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"),
                               externalHaplotypeFile = paste0(multisamplehaplotypeprefix, chrom, ".vcf"),
                               oldfilesuffix = "_noMulti.txt")
        
        GetChromosomeBAFs(chrom=chrom,
                          SNP_file=paste0(samplename[sampleidx], "_alleleFrequencies_chr", chrom, ".txt"),
                          haplotypeFile=paste0(samplename[sampleidx], "_impute_output_chr", chrom, "_allHaplotypeInfo.txt"),
                          samplename=samplename[sampleidx],
                          outfile=paste0(samplename[sampleidx], "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"),
                          chr_names=chrom_names,
                          minCounts=min_normal_depth)

	# Plot what we have until this point
        plot.haplotype.data(haplotyped.baf.file=paste0(samplename[sampleidx], "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt"),
                            imageFileName=paste0(samplename[sampleidx],"_chr",chrom,"_heterozygousData.png"),
                            samplename=samplename[sampleidx],
                            chrom=chrom,
                            chr_names=chrom_names)
      }
      
    }
    
    # Kill the threads as from here its single core
    parallel::stopCluster(clp)
    
    for (sampleidx in 1:nsamples) {
      
      # Combine all the BAF output into a single file
      combine.baf.files(inputfile.prefix=paste0(samplename[sampleidx], "_chr"),
                        inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt",
                        outputfile=paste0(samplename[sampleidx], "_heterozygousMutBAFs_haplotyped.txt"), 
                        chr_names=chrom_names)
      
    }
    # Segment the phased and haplotyped BAF data
    segment.baf.phased.multisample(samplename=samplename,
                                   inputfile=paste(samplename, "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
                                   outputfile=paste(samplename, ".BAFsegmented.txt", sep=""),
                                   prior_breakpoints_file=prior_breakpoints_file,
                                   gamma=segmentation_gamma_multisample,
                                   calc_seg_baf_option=calc_seg_baf_option,
                                   GENOMEBUILD=genomebuild)
    
  }
  
  # Setup for parallel computing
  clp = parallel::makeCluster(min(nthreads, nsamples),outfile="")
  doParallel::registerDoParallel(clp)
  # for (sampleidx in 1:nsamples) {
  foreach::foreach (sampleidx=1:nsamples) %dopar% {
    .libPaths(libs)
    print(paste0("Fitting final copy number and calling subclones for sample ", samplename[sampleidx]))
    
    if (data_type=="wgs" | data_type=="WGS") {
      logr_file = paste(samplename[sampleidx], "_mutantLogR_gcCorrected.tab", sep="")
      if (analysis=="paired") {
        allelecounts_file = paste(samplename[sampleidx], "_alleleCounts.tab", sep="")
      } else {
      # Not produced by a number of analysis and is required for some plots. Setting to NULL  makes the pipeline not attempt to create these plots
        allelecounts_file = NULL
      }
    }
    
    # Fit a clonal copy number profile
    fit.copy.number(samplename=samplename[sampleidx],
                    outputfile.prefix=paste(samplename[sampleidx], "_", sep=""),
                    inputfile.baf.segmented=paste(samplename[sampleidx], ".BAFsegmented.txt", sep=""),
                    inputfile.baf=paste(samplename[sampleidx],"_mutantBAF.tab", sep=""),
                    inputfile.logr=logr_file,
                    dist_choice=clonality_dist_metric,
                    ascat_dist_choice=ascat_dist_metric,
                    min.ploidy=min_ploidy,
                    max.ploidy=max_ploidy,
                    min.rho=min_rho,
                    max.rho=max_rho,
                    min.goodness=min_goodness,
                    uninformative_BAF_threshold=uninformative_BAF_threshold,
                    gamma_param=platform_gamma,
                    use_preset_rho_psi=F,
                    preset_rho=NA,
                    preset_psi=NA,
                    read_depth=30,
                    analysis=analysis,
                    nthreads=nthreads,
                    enhanced_grid_search=enhanced_grid_search)
    
    # Go over all segments, determine which segements are a mixture of two states and fit a second CN state
    print("callSubclones")
    callSubclones(sample.name=samplename[sampleidx],
                  baf.segmented.file=paste(samplename[sampleidx], ".BAFsegmented.txt", sep=""),
                  logr.file=logr_file,
                  rho.psi.file=paste(samplename[sampleidx], "_rho_and_psi.txt",sep=""),
                  output.file=paste(samplename[sampleidx],"_copynumber.txt", sep=""),
                  output.figures.prefix=paste(samplename[sampleidx],"_subclones_chr", sep=""),
                  output.gw.figures.prefix=paste(samplename[sampleidx],"_BattenbergProfile", sep=""),
                  masking_output_file=paste(samplename[sampleidx], "_segment_masking_details.txt", sep=""),
                  prior_breakpoints_file=prior_breakpoints_file,
                  chr_names=chrom_names, 
                  gamma=platform_gamma, 
                  segmentation.gamma=NA, 
                  siglevel=0.05, 
                  maxdist=0.01, 
                  max_allowed_state=max_allowed_state, 
                  cn_upper_limit=cn_upper_limit, 
                  noperms=1000,
                  calc_seg_baf_option=calc_seg_baf_option)
    
    # If patient is male, get copy number status of ChrX based only on logR segmentation (due to hemizygosity of SNPs)
    # Only do this when X chromosome is included
    if (ismale & "X" %in% chrom_names){
      print("callChrXsubclones")
      callChrXsubclones(tumourname=samplename[sampleidx],
                        X_gamma=1000,
                        X_kmin=100,
                        genomebuild=genomebuild,
                        AR=TRUE,
                        prior_breakpoints_file=prior_breakpoints_file,
			chrom_names=chrom_names,
                        data_type=data_type)
    }
    
    # Make some post-hoc plots
    print("make_posthoc_plots")
    make_posthoc_plots(samplename=samplename[sampleidx],
                       logr_file=logr_file,
                       bafsegmented_file=paste(samplename[sampleidx], ".BAFsegmented.txt", sep=""),
                       logrsegmented_file=paste(samplename[sampleidx], ".logRsegmented.txt", sep=""),
                       allelecounts_file=allelecounts_file)
    
    # Save refit suggestions for a future rerun
    print("cnfit_to_refit_suggestions")
    cnfit_to_refit_suggestions(samplename=samplename[sampleidx],
                               subclones_file=paste(samplename[sampleidx], "_copynumber_extended.txt", sep=""),
                               rho_psi_file=paste(samplename[sampleidx], "_rho_and_psi.txt", sep=""),
                               gamma_param=platform_gamma)
  }
  
  # Kill the threads as last part again is single core
  parallel::stopCluster(clp)
  
  if (nsamples > 1) {
    print("Assessing mirrored subclonal allelic imbalance (MSAI)")
    call_multisample_MSAI(rdsprefix = multisamplehaplotypeprefix,
                          subclonesfiles = paste0(samplename, "_copynumber_extended.txt"),
                          chrom_names = chrom_names,
                          tumournames = samplename,
                          plotting = T)
  }
}
