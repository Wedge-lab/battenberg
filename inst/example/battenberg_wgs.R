args = commandArgs(TRUE)
TUMOURNAME = toString(args[1])
NORMALNAME = toString(args[2])
NORMALBAM = toString(args[3])
TUMOURBAM = toString(args[4])
IS.MALE = as.logical(args[5])
RUN_DIR = toString(args[6])

library(Battenberg)
library(doParallel)

###############################################################################
# 2017-02-07
# A pure R Battenberg v2.2.3 WGS pipeline implementation.
# sd11@sanger.ac.uk
###############################################################################


# Sample specific
#TUMOURNAME = "PD7422a"
#NORMALNAME = "PD7422b"
#IS.MALE = F
#TUMOURBAM = "/lustre/scratch110/sanger/sd11/epitax/bam/PD7422a.bam"
#NORMALBAM = "/lustre/scratch110/sanger/sd11/epitax/bam/PD7422b.bam"
#RUN_DIR = getwd()

# Parallelism parameters
NTHREADS = 6

# General static
IMPUTEINFOFILE = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_impute_v3/impute_info.txt"
G1000PREFIX = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr"
G1000PREFIX_AC = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_1000genomesloci2012_v3/1000genomesloci2012_chr"
GCCORRECTPREFIX = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_wgs_gc_correction_1000g_v3/1000_genomes_GC_corr_chr_"
IMPUTE_EXE = "impute2"


PLATFORM_GAMMA = 1
PHASING_GAMMA = 1
SEGMENTATION_GAMMA = 10
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
CALC_SEG_BAF_OPTION = 1

# WGS specific static
ALLELECOUNTER = "alleleCounter"
PROBLEMLOCI = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_probloci/probloci_270415.txt.gz"

chrom_names = get.chrom.names(IMPUTEINFOFILE, IS.MALE)

# Setup for parallel computing
clp = makeCluster(NTHREADS)
registerDoParallel(clp)

# Obtain allele counts for 1000 Genomes locations for both tumour and normal
#for (i in 1:length(chrom_names)) {
foreach(i=1:length(chrom_names), .export=c("getAlleleCounts")) %dopar% {
  getAlleleCounts(bam.file=TUMOURBAM,
                  output.file=paste(TUMOURNAME,"_alleleFrequencies_chr", i, ".txt", sep=""),
                  g1000.loci=paste(G1000PREFIX_AC, i, ".txt", sep=""),
                  min.base.qual=MIN_BASE_QUAL,
                  min.map.qual=MIN_MAP_QUAL,
                  allelecounter.exe=ALLELECOUNTER)
  
  getAlleleCounts(bam.file=NORMALBAM,
                  output.file=paste(NORMALNAME,"_alleleFrequencies_chr", i, ".txt",  sep=""),
                  g1000.loci=paste(G1000PREFIX_AC, i, ".txt", sep=""),
                  min.base.qual=MIN_BASE_QUAL,
                  min.map.qual=MIN_MAP_QUAL,
                  allelecounter.exe=ALLELECOUNTER)

}
# Obtain BAF and LogR from the raw allele counts
getBAFsAndLogRs(tumourAlleleCountsFile.prefix=paste(TUMOURNAME,"_alleleFrequencies_chr", sep=""), 
                normalAlleleCountsFile.prefix=paste(NORMALNAME,"_alleleFrequencies_chr", sep=""),
                figuresFile.prefix=paste(TUMOURNAME, "_", sep=''),
                BAFnormalFile=paste(TUMOURNAME,"_normalBAF.tab", sep=""), 
                BAFmutantFile=paste(TUMOURNAME,"_mutantBAF.tab", sep=""), 
                logRnormalFile=paste(TUMOURNAME,"_normalLogR.tab", sep=""), 
                logRmutantFile=paste(TUMOURNAME,"_mutantLogR.tab", sep=""), 
                combinedAlleleCountsFile=paste(TUMOURNAME,"_alleleCounts.tab", sep=""),
                chr_names=chrom_names, 
                g1000file.prefix=G1000PREFIX, 
                minCounts=MIN_NORMAL_DEPTH, 
                samplename=TUMOURNAME)
# Perform GC correction
gc.correct.wgs(Tumour_LogR_file=paste(TUMOURNAME,"_mutantLogR.tab", sep=""), 
               outfile=paste(TUMOURNAME,"_mutantLogR_gcCorrected.tab", sep=""),
               correlations_outfile=paste(TUMOURNAME, "_GCwindowCorrelations.txt", sep=""),
               gc_content_file_prefix=GCCORRECTPREFIX, 
               chrom_names=chrom_names)


# These steps are independent and can be run in parallel. This could be done through multi-threading or splitting these up into separate jobs on a cluster.
foreach(chrom=1:length(chrom_names), .export=c("generate.impute.input.wgs","run.impute","combine.impute.output","GetChromosomeBAFs","plot.haplotype.data")) %dopar% {
  print(chrom)
  # Transform the allele counts into something that impute understands
  generate.impute.input.wgs(chrom=chrom, 
                            tumour.allele.counts.file=paste(TUMOURNAME,"_alleleFrequencies_chr", chrom, ".txt", sep=""), 
                            normal.allele.counts.file=paste(NORMALNAME,"_alleleFrequencies_chr", chrom, ".txt", sep=""), 
                            output.file=paste(TUMOURNAME, "_impute_input_chr", chrom, ".txt", sep=""), 
                            imputeinfofile=IMPUTEINFOFILE,
                            is.male=IS.MALE,
                            problemLociFile=PROBLEMLOCI, 
                            useLociFile=NA)

  # Run impute on the files
  run.impute(inputfile=paste(TUMOURNAME, "_impute_input_chr", chrom, ".txt", sep=""), 
             outputfile.prefix=paste(TUMOURNAME, "_impute_output_chr", chrom, ".txt", sep=""), 
             is.male=IS.MALE, 
             imputeinfofile=IMPUTEINFOFILE, 
             impute.exe=IMPUTE_EXE, 
             region.size=5000000, 
             chrom=chrom)
  
  # As impute runs in windows across a chromosome we need to assemble the output
  combine.impute.output(inputfile.prefix=paste(TUMOURNAME, "_impute_output_chr", chrom, ".txt", sep=""), 
                        outputfile=paste(TUMOURNAME, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep=""), 
                        is.male=IS.MALE, 
                        imputeinfofile=IMPUTEINFOFILE, 
                        region.size=5000000, 
                        chrom=chrom)

  # Transform the impute output into haplotyped BAFs
  GetChromosomeBAFs(chrom=chrom, 
                    SNP_file=paste(TUMOURNAME, "_alleleFrequencies_chr", chrom, ".txt", sep=""), 
                    haplotypeFile=paste(TUMOURNAME, "_impute_output_chr", chrom, "_allHaplotypeInfo.txt", sep=""), 
                    samplename=TUMOURNAME, 
                    outfile=paste(TUMOURNAME, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                    chr_names=chrom_names, 
                    minCounts=MIN_NORMAL_DEPTH)

  # Plot what we have until this point
  plot.haplotype.data(haplotyped.baf.file=paste(TUMOURNAME, "_chr", chrom, "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                      imageFileName=paste(TUMOURNAME,"_chr",chrom,"_heterozygousData.png",sep=""), 
                      samplename=TUMOURNAME, 
                      chrom=chrom, 
                      chr_names=chrom_names)

  # Cleanup temp Impute output
  unlink(paste(TUMOURNAME, "_impute_output_chr", chrom, "*K.txt*", sep=""))
}

# Kill the threads as from here its all single core
stopCluster(clp)

# Combine all the BAF output into a single file
combine.baf.files(inputfile.prefix=paste(TUMOURNAME, "_chr", sep=""), 
                  inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt", 
                  outputfile=paste(TUMOURNAME, "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                  no.chrs=length(chrom_names))

# Segment the phased and haplotyped BAF data
segment.baf.phased(samplename=TUMOURNAME,
                     inputfile=paste(TUMOURNAME, "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
                     outputfile=paste(TUMOURNAME, ".BAFsegmented.txt", sep=""),
                     gamma=SEGMENTATION_GAMMA,
                     phasegamma=PHASING_GAMMA,
                     kmin=3,
                     phasekmin=1,
                     calc_seg_baf_option=CALC_SEG_BAF_OPTION)

# Fit a clonal copy number profile
fit.copy.number(samplename=TUMOURNAME,
                outputfile.prefix=paste(TUMOURNAME, "_", sep=""),
                inputfile.baf.segmented=paste(TUMOURNAME, ".BAFsegmented.txt", sep=""), 
                inputfile.baf=paste(TUMOURNAME,"_mutantBAF.tab", sep=""), 
                inputfile.logr=paste(TUMOURNAME,"_mutantLogR_gcCorrected.tab", sep=""), 
                dist_choice=CLONALITY_DIST_METRIC, 
                ascat_dist_choice=ASCAT_DIST_METRIC, 
                min.ploidy=MIN_PLOIDY, 
                max.ploidy=MAX_PLOIDY, 
                min.rho=MIN_RHO, 
                min.goodness=MIN_GOODNESS_OF_FIT, 
                uninformative_BAF_threshold=BALANCED_THRESHOLD, 
                gamma_param=PLATFORM_GAMMA, 
                use_preset_rho_psi=F, 
                preset_rho=NA, 
                preset_psi=NA, 
                read_depth=30)

# Go over all segments, determine which segements are a mixture of two states and fit a second CN state
callSubclones(sample.name=TUMOURNAME, 
              baf.segmented.file=paste(TUMOURNAME, ".BAFsegmented.txt", sep=""), 
              logr.file=paste(TUMOURNAME,"_mutantLogR_gcCorrected.tab", sep=""), 
              rho.psi.file=paste(TUMOURNAME, "_rho_and_psi.txt",sep=""), 
              output.file=paste(TUMOURNAME,"_subclones.txt", sep=""), 
              output.figures.prefix=paste(TUMOURNAME,"_subclones_chr", sep=""), 
              output.gw.figures.prefix=paste(TUMOURNAME,"_BattenbergProfile", sep=""),
              masking_output_file=paste(TUMOURNAME, "_segment_masking_details.txt", sep=""),
              sv_breakpoints_file="NA",
              chr_names=chrom_names, 
              gamma=PLATFORM_GAMMA, 
              segmentation.gamma=NA, 
              siglevel=0.05, 
              maxdist=0.01, 
              noperms=1000,
              calc_seg_baf_option=CALC_SEG_BAF_OPTION)

# Make some post-hoc plots
make_posthoc_plots(samplename=TUMOURNAME, 
                   logr_file=paste(TUMOURNAME, "_mutantLogR_gcCorrected.tab", sep=""), 
                   subclones_file=paste(TUMOURNAME, "_subclones.txt", sep=""), 
                   rho_psi_file=paste(TUMOURNAME, "_rho_and_psi.txt", sep=""), 
                   bafsegmented_file=paste(TUMOURNAME, ".BAFsegmented.txt", sep=""), 
                   logrsegmented_file=paste(TUMOURNAME, ".logRsegmented.txt", sep=""), 
                   allelecounts_file=paste(TUMOURNAME, "_alleleCounts.tab", sep=""))

# Save refit suggestions for a future rerun
cnfit_to_refit_suggestions(samplename=TUMOURNAME, 
                           subclones_file=paste(TUMOURNAME, "_subclones.txt", sep=""),
                           rho_psi_file=paste(TUMOURNAME, "_rho_and_psi.txt", sep=""),
                           gamma_param=PLATFORM_GAMMA)

