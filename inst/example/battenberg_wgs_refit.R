library(Battenberg)
#library(dpclust3p)
library(optparse)

option_list = list(
  make_option(c("-a", "--analysis_type"), type="character", default="paired", help="Type of analysis to run: paired (tumour+normal), cell_line (only tumour), germline (only normal)", metavar="character"),
  make_option(c("-t", "--tumourname"), type="character", default=NULL, help="Samplename of the tumour", metavar="character"),
  make_option(c("--sex"), type="character", default=NULL, help="Sex of the sample", metavar="character"),
  make_option(c("-i", "--inputdir"), type="character", default=NULL, help="Directory where input will be read from", metavar="character"),
  make_option(c("--bp"), type="character", default=NULL, help="Optional two column file (chromosome and position) specifying prior breakpoints to be used during segmentation", metavar="character"),
  make_option(c("-g", "--ref_genome_build"), type="character", default="hg19", help="Reference genome build to which the reads have been aligned. Options are hg19 and hg38", metavar="character"),
  make_option(c("-r", "--rho"), type="numeric", default=NULL, help="Rho value to be used for refitting", metavar="character"),
  make_option(c("-p", "--psi"), type="numeric", default=NULL, help="Psi value to be used for refitting", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

tumourname = opt$tumourname
ismale = opt$sex=="male" | opt$sex=="Male"
preset_rho = opt$rho
preset_psi = opt$psi
prior_breakpoints_file = opt$bp
analysis = opt$analysis_type
inputdir = opt$inputdir
genomebuild = opt$ref_genome_build

logr_file = file.path(inputdir, paste0(tumourname, "_mutantLogR_gcCorrected.tab"))
baf_file = file.path(inputdir, paste0(tumourname,"_mutantBAF.tab"))
bafsegmented_file = file.path(inputdir, paste0(tumourname,".BAFsegmented.txt"))
subclones_file = file.path(inputdir, paste0(tumourname,"_subclones.txt"))


clonality_dist_metric = 0
ascat_dist_metric = 2
min_ploidy = 1.6
max_ploidy = 4.8
min_rho = 0.1
min_goodness = 0.63
uninformative_BAF_threshold = 0.51
platform_gamma = 1
calc_seg_baf_option = 1

if (genomebuild=="hg19") {
        impute_basedir = "/hps/research/gerstung/sdentro/reference/human/battenberg/"
        imputeinfofile = file.path(impute_basedir, "battenberg_impute_v3/impute_info.txt")

} else if (genomebuild=="hg38") {

        beagle_basedir = "/lustre/scratch117/casm/team219/sd11/reference/human/battenberg_hg38"
        imputeinfofile = file.path(beagle_basedir, "imputation/impute_info.txt")
}

#outdir = paste0(tumourname, "_rho", preset_rho, "_psi", preset_psi)
#dir.create(outdir, showWarnings=F)
#setwd(outdir)


chrom_names = get.chrom.names(imputeinfofile, ismale, analysis=analysis)

Battenberg::fit.copy.number(samplename=tumourname,
                    outputfile.prefix=paste(tumourname, "_", sep=""),
                    inputfile.baf.segmented=bafsegmented_file,
                    inputfile.baf=baf_file,
                    inputfile.logr=logr_file,
                    dist_choice=clonality_dist_metric,
                    ascat_dist_choice=ascat_dist_metric,
                    min.ploidy=min_ploidy,
                    max.ploidy=max_ploidy,
                    min.rho=min_rho,
                    min.goodness=min_goodness,
                    uninformative_BAF_threshold=uninformative_BAF_threshold,
                    gamma_param=platform_gamma,
                    use_preset_rho_psi=T,
                    preset_rho=preset_rho,
                    preset_psi=preset_psi,
                    read_depth=30,
                    analysis=analysis)
    
    # Go over all segments, determine which segements are a mixture of two states and fit a second CN state
Battenberg::callSubclones(sample.name=tumourname,
                  baf.segmented.file=bafsegmented_file,
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
    
    # If patient is male, get copy number status of ChrX based only on logR segmentation (due to hemizygosity of SNPs)
    # Only do this when X chromosome is included
    if (ismale & "X" %in% chrom_names){
      Battenberg::callChrXsubclones(TUMOURNAME=tumourname[sampleidx],
                        X_GAMMA=1000,
                        X_KMIN=100,
                        GENOMEBUILD=GENOMEBUILD,
                        AR=TRUE)
    }
