
#' Run the Battenberg pipeline
#' 
#' @param tumourname
#' @param normalname
#' @param tumour_data_file A BAM or CEL file
#' @param normal_data_file A BAM or CEL file
#' @param ismale
#' @param imputeinfofile
#' @param g1000prefix
#' @param g1000allelesprefix
#' @param gccorrectprefix
#' @param problemloci
#' @param data_type Default: wgs
#' @param impute_exe Default: impute2
#' @param allelecounter_exe Default: alleleCounter
#' @param nthreads Default: 8
#' @param platform_gamma Default: 1
#' @param phasing_gamma Default: 1
#' @param segmentation_gamma Default: 10
#' @param segmentation_kmin Default: 3
#' @param phasing_kmin Default: 1
#' @param clonality_dist_metric Default: 0
#' @param ascat_dist_metric Default: 1
#' @param min_ploidy Default: 1.6
#' @param max_ploidy Default: 4.8
#' @param min_rho Default: 0.1
#' @param min_goodness Default: 0.63
#' @param uninformative_BAF_threshold Default: 0.51
#' @param min_normal_depth Default: 10
#' @param min_base_qual Default: 20
#' @param min_map_qual Default: 35
#' @param calc_seg_baf_option Default: 1
#' @param skip_allele_counting Default: FALSE
#' @param skip_preprocessing Default: FALSE
#' @param snp6_reference_info_file SNP6 pipeline only Default: NA
#' @param apt.probeset.genotype.exe SNP6 pipeline only Default: apt-probeset-genotype
#' @param apt.probeset.summarize.exe SNP6 pipeline only Default: apt-probeset-summarize
#' @param norm.geno.clust.exe SNP6 pipeline only Default: normalize_affy_geno_cluster.pl
#' @param birdseed_report_file SNP6 pipeline only Default: birdseed.report.txt
#' @author sd11
#' @export
battenberg = function(tumourname, normalname, tumour_data_file, normal_data_file, ismale, imputeinfofile, g1000prefix, g1000allelesprefix, gccorrectprefix, problemloci, 
                      data_type="wgs", impute_exe="impute2", allelecounter_exe="alleleCounter", nthreads=8, platform_gamma=1, phasing_gamma=1,
                      segmentation_gamma=10, segmentation_kmin=3, phasing_kmin=1, clonality_dist_metric=0, ascat_dist_metric=1, min_ploidy=1.6,
                      max_ploidy=4.8, min_rho=0.1, min_goodness=0.63, uninformative_BAF_threshold=0.51, min_normal_depth=10, min_base_qual=20, 
                      min_map_qual=35, calc_seg_baf_option=1, skip_allele_counting=F, skip_preprocessing=F,
                      snp6_reference_info_file=NA, apt.probeset.genotype.exe="apt-probeset-genotype", apt.probeset.summarize.exe="apt-probeset-summarize", 
                      norm.geno.clust.exe="normalize_affy_geno_cluster.pl", birdseed_report_file="birdseed.report.txt") {

  # Parallelism parameters
  # NTHREADS = 6
  
  # General static
  # IMPUTEINFOFILE = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_impute_v3/impute_info.txt"
  # G1000PREFIX = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr"
  # G1000PREFIX_AC = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr"
  # GCCORRECTPREFIX = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_"
  # IMPUTE_EXE = "impute2"
  
  
  # PLATFORM_GAMMA = 1
  # PHASING_GAMMA = 1
  # SEGMENTATION_GAMMA = 10
  # CLONALITY_DIST_METRIC = 0
  # ASCAT_DIST_METRIC = 1
  # MIN_PLOIDY = 1.6
  # MAX_PLOIDY = 4.8
  # MIN_RHO = 0.1
  # MIN_GOODNESS_OF_FIT = 0.63
  # BALANCED_THRESHOLD = 0.51
  # MIN_NORMAL_DEPTH = 10
  # MIN_BASE_QUAL = 20
  # MIN_MAP_QUAL = 35
  # CALC_SEG_BAF_OPTION = 1
  
  # WGS specific static
  # ALLELECOUNTER = "alleleCounter"
  # PROBLEMLOCI = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_probloci/probloci_270415.txt.gz"
  
  chrom_names = get.chrom.names(imputeinfofile, ismale)
  
  # Setup for parallel computing
  clp = makeCluster(nthreads)
  registerDoParallel(clp)
  
  if (!skip_preprocessing) {
    if (data_type=="wgs" | data_type=="WGS") {
      
      prepare_wgs(chrom_names=chrom_names, 
                  tumourbam=tumour_data_file, 
                  normalbam=normal_data_file, 
                  tumourname=tumourname, 
                  normalname=normalname, 
                  g1000allelesprefix=g1000allelesprefix, 
                  g1000prefix=g1000prefix, 
                  gccorrectprefix=gccorrectprefix, 
                  min_base_qual=min_base_qual, 
                  min_map_qual=min_map_qual, 
                  allelecounter_exe=allelecounter_exe, 
                  min_normal_depth=min_normal_depth, 
                  nthreads=nthreads,
                  skip_allele_counting=skip_allele_counting)
      
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
      # # Extract the LogR and BAF from both tumour and normal cel files.
      # cel2baf.logr(normal_cel_file=NORMALCEL, 
      #              tumour_cel_file=TUMOURCEL, 
      #              output_file=paste(tumourname, "_lrr_baf.txt", sep=""), 
      #              snp6_reference_info_file=SNP6_REF_INFO_FILE, 
      #              apt.probeset.genotype.exe=APT_PROBESET_GENOTYPE_EXE, 
      #              apt.probeset.summarize.exe=APT_PROBESET_SUMMARIZE_EXE, 
      #              norm.geno.clust.exe=NORM_GENO_CLUST_EXE)
      # 
      # gc.correct(samplename=tumourname, 
      #            infile.logr.baf=paste(tumourname, "_lrr_baf.txt", sep=""), 
      #            outfile.tumor.LogR=paste(tumourname, "_mutantLogR.tab", sep=""), 
      #            outfile.tumor.BAF=paste(tumourname, "_mutantBAF.tab", sep=""), 
      #            outfile.normal.LogR=paste(tumourname, "_germlineLogR.tab", sep=""), 
      #            outfile.normal.BAF=paste(tumourname, "_germlineBAF.tab", sep=""), 
      #            outfile.probeBAF=paste(tumourname, "_probeBAF.txt", sep=""),
      #            snp6_reference_info_file=SNP6_REF_INFO_FILE,
      #            birdseed_report_file=BIRDSEED_REPORT_FILE,
      #            chr_names=chrom_names)
      
      
      
      
    } else {
      print("Unknown data type provided, please provide wgs or snp6")
      q(save="no", status=1)
    }
  }
 
  # Reconstruct haplotypes 
  mclapply(1:length(chrom_names), function(chrom) {
    print(chrom)
    
    run_haplotyping(chrom=chrom, 
                    tumourname=tumourname, 
                    normalname=normalname, 
                    ismale=ismale, 
                    imputeinfofile=imputeinfofile, 
                    problemloci=problemloci, 
                    impute_exe=impute_exe, 
                    min_normal_depth=min_normal_depth,
		    chrom_names=chrom_names)
  }, mc.cores=nthreads)
  
  # Kill the threads as from here its all single core
  stopCluster(clp)
  
  # Combine all the BAF output into a single file
  combine.baf.files(inputfile.prefix=paste(tumourname, "_chr", sep=""), 
                    inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt", 
                    outputfile=paste(tumourname, "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                    no.chrs=length(chrom_names))
  
  # Segment the phased and haplotyped BAF data
  segment.baf.phased(samplename=tumourname,
                     inputfile=paste(tumourname, "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
                     outputfile=paste(tumourname, ".BAFsegmented.txt", sep=""),
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
                  inputfile.logr=paste(tumourname,"_mutantLogR_gcCorrected.tab", sep=""), 
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
                logr.file=paste(tumourname,"_mutantLogR_gcCorrected.tab", sep=""), 
                rho.psi.file=paste(tumourname, "_rho_and_psi.txt",sep=""), 
                output.file=paste(tumourname,"_subclones.txt", sep=""), 
                output.figures.prefix=paste(tumourname,"_subclones_chr", sep=""), 
                output.gw.figures.prefix=paste(tumourname,"_BattenbergProfile", sep=""),
                masking_output_file=paste(tumourname, "_segment_masking_details.txt", sep=""),
                sv_breakpoints_file="NA",
                chr_names=chrom_names, 
                gamma=platform_gamma, 
                segmentation.gamma=NA, 
                siglevel=0.05, 
                maxdist=0.01, 
                noperms=1000,
                calc_seg_baf_option=CALC_SEG_BAF_OPTION)
  
  # Make some post-hoc plots
  make_posthoc_plots(samplename=tumourname, 
                     logr_file=paste(tumourname, "_mutantLogR_gcCorrected.tab", sep=""), 
                     subclones_file=paste(tumourname, "_subclones.txt", sep=""), 
                     rho_psi_file=paste(tumourname, "_rho_and_psi.txt", sep=""), 
                     bafsegmented_file=paste(tumourname, ".BAFsegmented.txt", sep=""), 
                     logrsegmented_file=paste(tumourname, ".logRsegmented.txt", sep=""), 
                     allelecounts_file=paste(tumourname, "_alleleCounts.tab", sep=""))
  
  # Save refit suggestions for a future rerun
  cnfit_to_refit_suggestions(samplename=tumourname, 
                             subclones_file=paste(tumourname, "_subclones.txt", sep=""),
                             rho_psi_file=paste(tumourname, "_rho_and_psi.txt", sep=""),
                             gamma_param=PLATFORM_GAMMA)
  
  
  
  
    
}

