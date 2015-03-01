# Parameters:
TUMOURNAME = "PD7422a"
NORMALNAME = "PD7422b"
IS.MALE = F
TUMOURBAM = "/lustre/scratch110/sanger/sd11/epitax/bam/PD7422a.bam"
NORMALBAM = "/lustre/scratch110/sanger/sd11/epitax/bam/PD7422b.bam"

IMPUTEINFOFILE = "/nfs/users/nfs_s/sd11/repo/battenberg/impute_info.txt"
G1000PREFIX = "/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_1000genomesloci2012/1000genomesAlleles2012_chr"
PROBLEMLOCI = "/nfs/users/nfs_s/sd11/repo/battenberg/probloci.txt"
ALLELECOUNTER = "alleleCounter"
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

chrom_names = get.chrom.names(IMPUTEINFOFILE, IS.MALE)

# Obtain allele counts for 1000 Genomes locations for both tumour and normal
for (i in 1:length(chrom_names)) {
  getAlleleCounts(bam.file=TUMOURBAM,
                  output.file=paste(TUMOURNAME,"_alleleFrequencies_chr", sep=""),
                  g1000.loci=paste(G1000PREFIX, i, ".txt", sep=""),
                  min.base.qual=MIN_BASE_QUAL,
                  min.map.qual=MIN_MAP_QUAL,
                  allelecounter.exe=ALLELECOUNTER)
  
  getAlleleCounts(bam.file=NORMALBAM,
                  output.file=paste(NORMALNAME,"_alleleFrequencies_chr", sep=""),
                  g1000.loci=paste(G1000PREFIX, i, ".txt", sep=""),
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

# These steps are independent and can be run in parallel. This could be done through multi-threading or splitting these up into separate jobs on a cluster.
for (chrom in 1:length(chrom_names)) {
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
             impute.exe="impute2", 
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
                      imageFileName=paste(TUMOURNAME,"_chr",chr_name,"_heterozygousData.png",sep=""), 
                      samplename=TUMOURNAME, 
                      chrom=chrom, 
                      chr_names=chrom_names)
}

# Combine all the BAF output into a single file
# TODO: Clean up the parameters here
combine.baf.files(inputfile.prefix=paste(TUMOURNAME, "_chr", sep=""), 
                  inputfile.postfix="_heterozygousMutBAFs_haplotyped.txt", 
                  outputfile=paste(TUMOURNAME, "_heterozygousMutBAFs_haplotyped.txt", sep=""),
                  no.chrs=length(chrom_names),
                  gamma=10,
                  phasegamma=3,
                  kmin=3,
                  phasekmin=3)

# Segment the phased and haplotyped BAF data
segment.baf.phased(inputfile=paste(TUMOURNAME, "_heterozygousMutBAFs_haplotyped.txt", sep=""), 
                   outputfile=paste(TUMOURNAME, ".BAFsegmented.txt", sep=""))

# Fit a clonal copy number profile
# TODO: Check parameters, clean up internally defined output names
fit.copy.number(samplename=TUMOURNAME,
                outputfile.prefix=paste(TUMOURNAME, "_", sep=""),
                inputfile.baf.segmented=paste(TUMOURNAME, ".BAFsegmented.txt", sep=""), 
                inputfile.baf=paste(TUMOURNAME,"_mutantBAF.tab", sep=""), 
                inputfile.logr=paste(TUMOURNAME,"_mutantLogR.tab", sep=""), 
                dist_choice=CLONALITY_DIST_METRIC, 
                ascat_dist_choice=ASCAT_DIST_METRIC, 
                min.ploidy=MIN_PLOIDY, 
                max.ploidy=MAX_PLOIDY, 
                min.rho=MIN_RHO, 
                min.goodness=MIN_GOODNESS_OF_FIT, 
                uninformative_BAF_threshold=BALANCED_THRESHOLD, 
                gamma_param=1, 
                use_preset_rho_psi=F, 
                preset_rho=NA, 
                preset_psi=NA, 
                read_depth=30)

# Go over all segments, determine which segements are a mixture of two states and fit a second CN state
callSubclones(sample.name=TUMOURNAME, 
              baf.segmented.file=paste(TUMOURNAME, ".BAFsegmented.txt", sep=""), 
              logr.file=paste(TUMOURNAME,"_mutantLogR.tab", sep=""), 
              rho.psi.file=paste(TUMOURNAME, "_rho_and_psi.txt",sep=""), 
              output.file=paste(TUMOURNAME,"_subclones.txt", sep=""), 
              output.figures.prefix=paste(TUMOURNAME,"_subclones_chr", sep=""), 
              chr_names=chrom_names, 
              gamma=1, 
              segmentation.gamma=NA, 
              siglevel=0.05, 
              maxdist=0.01, 
              noperms=1000)