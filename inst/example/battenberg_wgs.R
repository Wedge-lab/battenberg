library(Battenberg)
library(optparse)

option_list = list(
  make_option(c("-a", "--analysis_type"), type="character", default="paired", help="Type of analysis to run: paired (tumour+normal), cell_line (only tumour), germline (only normal)", metavar="character"),
  make_option(c("-t", "--samplename"), type="character", default=NULL, help="Samplename of the tumour", metavar="character"),
  make_option(c("-n", "--normalname"), type="character", default=NULL, help="Samplename of the normal", metavar="character"),
  make_option(c("--tb"), type="character", default=NULL, help="Sample BAM file", metavar="character"),
  make_option(c("--nb"), type="character", default=NULL, help="Normal BAM file", metavar="character"),
  make_option(c("--beagle_jar"), type="character", default=NULL, help="Full path to beagle jar", metavar="character"),
  make_option(c("--beagle_ref_template"), type="character", default=NULL, help="Full path to beagle reference template", metavar="character"),
  make_option(c("--beagle_plink_template"), type="character", default=NULL, help="Full path to beagle plink maps template", metavar="character"),
  make_option(c("--sex"), type="character", default=NULL, help="Sex of the sample", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Directory where output will be written", metavar="character"),
  make_option(c("--skip_allelecount"), type="logical", default=FALSE, action="store_true", help="Provide when alleles don't have to be counted. This expects allelecount files on disk", metavar="character"),
  make_option(c("--skip_preprocessing"), type="logical", default=FALSE, action="store_true", help="Provide when pre-processing has previously completed. This expects the files on disk", metavar="character"),
  make_option(c("--skip_phasing"), type="logical", default=FALSE, action="store_true", help="Provide when phasing has previously completed. This expects the files on disk", metavar="character"),
  make_option(c("--cpu"), type="numeric", default=8, help="The number of CPU cores to be used by the pipeline (Default: 8)", metavar="character"),
  make_option(c("--bp"), type="character", default=NULL, help="Optional two column file (chromosome and position) specifying prior breakpoints to be used during segmentation", metavar="character"),
  make_option(c("-g", "--ref_genome_build"), type="character", default="hg19", help="Reference genome build to which the reads have been aligned. Options are hg19 and hg38", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

analysis = opt$analysis_type
SAMPLENAME = opt$samplename
NORMALNAME = opt$normalname
NORMALBAM = opt$nb
SAMPLEBAM = opt$tb
BEAGLEJAR = opt$beagle_jar
BEAGLEREF.template = opt$beagle_ref_template
BEAGLEPLINK.template = opt$beagle_plink_template
IS.MALE = opt$sex=="male" | opt$sex=="Male"
RUN_DIR = opt$output
SKIP_ALLELECOUNTING = opt$skip_allelecount
SKIP_PREPROCESSING = opt$skip_preprocessing
SKIP_PHASING = opt$skip_phasing
NTHREADS = opt$cpu
PRIOR_BREAKPOINTS_FILE = opt$bp
GENOMEBUILD = opt$ref_genome_build

supported_analysis = c("paired", "cell_line", "germline")
if (!analysis %in% supported_analysis) {
        stop(paste0("Requested analysis type ", analysis, " is not available. Please provide either of ", paste(supported_analysis, collapse=" ")))
}

supported_genome_builds = c("hg19", "hg38")
if (!GENOMEBUILD %in% supported_genome_builds) {
        stop(paste0("Provided genome build ", GENOMEBUILD, " is not supported. Please provide either of ", paste(supported_genome_builds, collapse=" ")))
}

###############################################################################
# 2022-01-07
# A pure R Battenberg v2.2.10 WGS pipeline implementation.
###############################################################################

JAVAJRE = "java"
ALLELECOUNTER = "alleleCounter"
IMPUTE_EXE = "impute2"
USEBEAGLE = T

# General static
if (GENOMEBUILD=="hg19") {
	impute_basedir = "/hps/research/gerstung/sdentro/reference/human/battenberg/"
	IMPUTEINFOFILE = file.path(impute_basedir, "battenberg_impute_v3/impute_info.txt")
	G1000PREFIX_AC = file.path(impute_basedir, "battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr")
	G1000PREFIX = file.path(impute_basedir, "battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr")
	GCCORRECTPREFIX = file.path(impute_basedir, "battenberg_wgs_gc_correction_1000g_v3_noNA/1000_genomes_GC_corr_chr_")
	REPLICCORRECTPREFIX = file.path(impute_basedir, "battenberg_wgs_replic_correction_1000g_v3/1000_genomes_replication_timing_chr_")
	
	# WGS specific static
	PROBLEMLOCI = "/hps/research/gerstung/sdentro/reference/human/battenberg/battenberg_probloci/probloci_270415.txt.gz"
	GENOME_VERSION = "b37"
	GENOMEBUILD = "hg19"
	BEAGLE_BASEDIR = "/hps/research/gerstung/sdentro/reference/human/battenberg/battenberg_beagle"
	BEAGLEJAR = file.path(BEAGLE_BASEDIR, "beagle.24Aug19.3e8.jar")
	BEAGLEREF.template = file.path(BEAGLE_BASEDIR, GENOME_VERSION, "chrCHROMNAME.1kg.phase3.v5a.b37.bref3")
	BEAGLEPLINK.template = file.path(BEAGLE_BASEDIR, GENOME_VERSION, "plink.chrCHROMNAME.GRCh37.map")

	CHROM_COORD_FILE = "/homes/sdentro/repo/battenberg/gcCorrect_chromosome_coordinates_hg19.txt"

} else if (GENOMEBUILD=="hg38") {
	
	BEAGLE_BASEDIR = "/hps/research/gerstung/sdentro/reference/human/battenberg_hg38"
	GENOMEBUILD = "hg38"
	IMPUTEINFOFILE = file.path(BEAGLE_BASEDIR, "imputation/impute_info.txt")
	G1000PREFIX_AC = file.path(BEAGLE_BASEDIR, "1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_allele_index_chr")
	G1000PREFIX = file.path(BEAGLE_BASEDIR, "1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chr")
	GCCORRECTPREFIX = file.path(BEAGLE_BASEDIR, "GC_correction_hg38/1000G_GC_chr")
	REPLICCORRECTPREFIX = file.path(BEAGLE_BASEDIR, "RT_correction_hg38/1000G_RT_chr")
	PROBLEMLOCI = file.path(BEAGLE_BASEDIR, "probloci/probloci.txt.gz")
	
	BEAGLEREF.template = file.path(BEAGLE_BASEDIR, "beagle5/chrCHROMNAME.1kg.phase3.v5a_GRCh38nounref.vcf.gz")
	BEAGLEPLINK.template = file.path(BEAGLE_BASEDIR, "beagle5/plink.chrCHROMNAME.GRCh38.map")
	BEAGLEJAR = file.path(BEAGLE_BASEDIR, "beagle.18May20.d20.jar")

	CHROM_COORD_FILE = "/homes/sdentro/repo/battenberg/chromosome_coordinates_hg38.txt"
} 

PLATFORM_GAMMA = 1
PHASING_GAMMA = 1
SEGMENTATION_GAMMA = 10
SEGMENTATIIN_KMIN = 3
PHASING_KMIN = 1
CLONALITY_DIST_METRIC = 0
ASCAT_DIST_METRIC = 1
MIN_PLOIDY = 1.6
MAX_PLOIDY = 4.8
MIN_RHO = 0.1
MIN_GOODNESS_OF_FIT = 0.63
BALANCED_THRESHOLD = 0.51
MIN_NORMAL_DEPTH = 10
MIN_BASE_QUAL = 20
MIN_MAP_QUAL = 35
#CALC_SEG_BAF_OPTION = 1
CALC_SEG_BAF_OPTION = 3
USEBEAGLE=TRUE
BEAGLE_MAX_MEM=15
BEAGLENTHREADS=1
BEAGLEWINDOW=40
BEAGLEOVERLAP=4

# Enable cairo device (needed to prevent 'X11 not available' errors)
options(bitmapType='cairo')

# Change to work directory and load the chromosome information
setwd(RUN_DIR)

battenberg(analysis=analysis,
	   samplename=SAMPLENAME,
           normalname=NORMALNAME,
           sample_data_file=SAMPLEBAM,
           normal_data_file=NORMALBAM, 
           ismale=IS.MALE, 
           imputeinfofile=IMPUTEINFOFILE, 
           g1000prefix=G1000PREFIX,
           g1000allelesprefix=G1000PREFIX_AC, 
           gccorrectprefix=GCCORRECTPREFIX, 
           repliccorrectprefix=REPLICCORRECTPREFIX, 
           problemloci=PROBLEMLOCI, 
           data_type="wgs",
           impute_exe=IMPUTE_EXE,
           allelecounter_exe=ALLELECOUNTER,
	   usebeagle=USEBEAGLE, ##set to TRUE to use beagle
           beaglejar=BEAGLEJAR, ##path
           beagleref=BEAGLEREF.template, ##pathtemplate
           beagleplink=BEAGLEPLINK.template, ##pathtemplate
           beaglemaxmem=BEAGLE_MAX_MEM,
           beaglenthreads=BEAGLENTHREADS,
           beaglewindow=BEAGLEWINDOW,
           beagleoverlap=BEAGLEOVERLAP,
           javajre=JAVAJRE,
           nthreads=NTHREADS,
           platform_gamma=PLATFORM_GAMMA,
           phasing_gamma=PHASING_GAMMA,
           segmentation_gamma=SEGMENTATION_GAMMA,
           segmentation_kmin=SEGMENTATIIN_KMIN,
           phasing_kmin=PHASING_KMIN,
           clonality_dist_metric=CLONALITY_DIST_METRIC,
           ascat_dist_metric=ASCAT_DIST_METRIC,
           min_ploidy=MIN_PLOIDY,
           max_ploidy=MAX_PLOIDY,
           min_rho=MIN_RHO,
           min_goodness=MIN_GOODNESS_OF_FIT,
           uninformative_BAF_threshold=BALANCED_THRESHOLD,
           min_normal_depth=MIN_NORMAL_DEPTH,
           min_base_qual=MIN_BASE_QUAL,
           min_map_qual=MIN_MAP_QUAL,
           calc_seg_baf_option=CALC_SEG_BAF_OPTION,
           skip_allele_counting=SKIP_ALLELECOUNTING,
           skip_preprocessing=SKIP_PREPROCESSING,
           skip_phasing=SKIP_PHASING,
           prior_breakpoints_file=PRIOR_BREAKPOINTS_FILE,
           genomebuild=GENOMEBUILD,
           chrom_coord_file=CHROM_COORD_FILE)
