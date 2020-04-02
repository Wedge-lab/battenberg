library(Battenberg)
library(optparse)

option_list = list(
  make_option(c("-t", "--tumourname"), type="character", default=NULL, help="Samplename of the tumour", metavar="character"),
  make_option(c("-n", "--normalname"), type="character", default=NULL, help="Samplename of the normal", metavar="character"),
  make_option(c("--tb"), type="character", default=NULL, help="Tumour BAM file", metavar="character"),
  make_option(c("--nb"), type="character", default=NULL, help="Normal BAM file", metavar="character"),
  make_option(c("--sex"), type="character", default=NULL, help="Sex of the sample", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Directory where output will be written", metavar="character"),
  make_option(c("--cpu"), type="numeric", default=8, help="The number of CPU cores to be used by the pipeline (Default: 8)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

tumourname = opt$tumourname
normalname = opt$normalname
normalbam = opt$nb
tumourbam = opt$tb
ismale = opt$sex=="male" | opt$sex=="Male"
run_dir = opt$output
nthreads = opt$cpu

###############################################################################
# 2019-04-29
# A pure R Battenberg v2.2.9 allele counting pipeline implementation.
# sd11 [at] sanger.ac.uk
###############################################################################

# Reference files
imputeinfofile = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_impute_v3/impute_info.txt"
g1000allelesprefix = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr"

# allele counter parameters
min_base_qual = 20
min_map_qual = 35
allelecounter_exe = "alleleCounter"

setwd(run_dir)

# get all required chromosomes
chrom_names = get.chrom.names(imputeinfofile, ismale)

# Parallel computing setup
clp = parallel::makeCluster(nthreads)
doParallel::registerDoParallel(clp)

# run allele counter
foreach::foreach(i=1:length(chrom_names)) %dopar% {
  getAlleleCounts(bam.file=tumourbam,
                  output.file=paste(tumourname,"_alleleFrequencies_chr", i, ".txt", sep=""),
                  g1000.loci=paste(g1000allelesprefix, i, ".txt", sep=""),
                  min.base.qual=min_base_qual,
                  min.map.qual=min_map_qual,
                  allelecounter.exe=allelecounter_exe)
  
  getAlleleCounts(bam.file=normalbam,
                  output.file=paste(normalname,"_alleleFrequencies_chr", i, ".txt",  sep=""),
                  g1000.loci=paste(g1000allelesprefix, i, ".txt", sep=""),
                  min.base.qual=min_base_qual,
                  min.map.qual=min_map_qual,
                  allelecounter.exe=allelecounter_exe)
}

# Kill the threads
parallel::stopCluster(clp)
