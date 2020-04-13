# Battenberg

This repository contains code for the whole genome sequencing subclonal copy number caller Battenberg, as described in [Nik-Zainal, Van Loo, Wedge, et al. (2012), Cell](https://www.ncbi.nlm.nih.gov/pubmed/22608083).

## Installation instructions

The instructions below will install the latest stable Battenberg version.

#### Prerequisites

Installing from Github requires devtools and Battenberg requires readr, splines, RColorBrewer and ASCAT, while the pipeline requires parallel and doParallel. From the command line run:

```
R -q -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("devtools", "splines", "readr", "doParallel", "ggplot2", "RColorBrewer", "gridExtra", "gtools", "parallel"));'
R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")'
```

#### Installation from Github

To install Battenberg, run the following from the command line:

```
R -q -e 'devtools::install_github("Wedge-Oxford/battenberg")'
```

#### Required reference files

Battenberg requires reference files, for now for GRCh37 only, that can be downloaded from here: https://ora.ox.ac.uk/objects/uuid:2c1fec09-a504-49ab-9ce9-3f17bac531bc

The bundle contains the following files:

  * battenberg_1000genomesloci2012_v3.tar.gz
  * battenberg_impute_1000G_v3.tar.gz
  * probloci_270415.txt.gz
  * battenberg_wgs_gc_correction_1000g_v3.tar.gz
  * battenberg_wgs_replic_correction_1000g_v3.tar.gz
  * battenberg_snp6_exe.tgz (SNP6 only)
  * battenberg_snp6_ref.tgz (SNP6 only)
  
#### Pipeline

Go into ```inst/example``` for example WGS and SNP6 R-only pipelines.

## Description of the output

### Key output files

* `[samplename]_subclones.txt` contains the copy number data (see table below)
* `[samplename]_rho_and_psi.txt` contains the purity estimate (make sure to use the FRAC_genome, rho field in the second row, first column)
* `[samplename]_BattenbergProfile*png` shows the profile (the two variants show subclonal copy number in a different way)
* `[samplename]_subclones_chr*.png` show detailed figures of the copy number calls per chromosome
* `[samplename]_distance.png` This shows the purity and ploidy solution space and can be used to pick alternative solutions

The copy number profile saved in the `[samplename]_subclones.txt` is a tab delimited file in text format. Within this file there is a line for each segment in the tumour genome.
Each segment will have either one or two copy number states:

* If there is one state that line represents the clonal copy number (i.e. all tumour cells have this state)
* If there are two states that line represents subclonal copy number (i.e. there are two populations of cells, each with a different state)

A copy number state consists of a major and a minor allele and their frequencies, which together add give the total copy number for that segment and an estimate fraction of tumour cells that carry each allele.

The following columns are available in the Battenberg output:

| Column | Description |
| ------------- | ------------- |
| chr | The chromosome of the segment |
| startpos | Start position on the chromosome |
| endpos | End position on the chromosome |
| BAF | The B-allele frequency of the segment |
| pval | P-value that is obtained when testing whether this segment should be represented by one or two states. A low p-value will result in the fitting of a second copy number state |
| LogR | The log ratio of normalised tumour coverage versus its matched normal sequencing sample |
| ntot | An internal total copy number value used to determine the priority of solutions. NOTE: This is not the total copy number of this segment! |
| nMaj1_A | The major allele copy number of state 1 from solution A |
| nMin1_A | The minor allele copy number of state 1 from solution A |
| frac1_A | Fraction of tumour cells carrying state 1 in solution A |
| nMaj2_A | The major allele copy number of state 2 from solution A. This value can be NA |
| nMin2_A | The minor allele copy number of state 2 from solution A. This value can be NA |
| frac2_A | Fraction of tumour cells carrying state 2 in solution A. This value can be NA |
| SDfrac_A | Standard deviation on the BAF of SNPs in this segment, can be used as a measure of uncertainty |
| SDfrac_A_BS | Bootstrapped standard deviation |
| frac1_A_0.025 | Associated 95% confidence interval of the bootstrap measure of uncertainty |

Followed by possible equivalent solutions B to F with the same columns as defined above for solution A (due to the way a profile is fit Battenberg can generate a series of equivalent solutions that are reported separately in the output).

### Plots for QC

It also produces a number plots that show the raw data and are useful for QC (and their raw data files denoted by *.tab)

* `[samplename].tumour.png` and `[samplename].germline.png` show the raw BAF and logR
* `[samplename]_coverage.png` contains coverage divided by the mean coverage of both tumour and normal
* `[samplename]_alleleratio.png` shows BAF*logR, a rough approximation of what the data looks like shortly before copy number calling

### Intermediate figures

Finally, a range of plots show intermediate steps and can occasionally be useful

* `[samplename]_chr*_heterozygousData.png` shows reconstructed haplotype blocks in the characteristic Battenberg cake pattern
* `[samplename]_RAFseg_chr*.png` and `[samplename]_segment_chr*.png` contains segmentation data for step 1 and step 2 respectively
* `[samplename]_nonroundedprofile.png` shows the copy number profile without rounding to integers
* `[samplename]_copynumberprofile.png` shows the copy number profile with (including subclonal copy number) rounding to integers

## Advice for including structural variant breakpoints

Battenberg can take prior breakpoints, from structural variants (SVs) for example, as input. SV breakpoints are typically much more precise and a pair of SVs can be closer together then what typically can be obtained from a BAF or coverage track. It is therefore adventageous to include prior breakpoints in a Battenberg run. However, including too many (as in 100s) incorrect breakpoints can have adverse effects by allowing many small segments to be affected by noise where there isn't any signal and increasing the runtime of the pipeline. It is therefore advised to `filter prior breakpoints from SVs such that the genome is slightly oversegmented.` Finally, some SV types, such as inversions, do not constitute a change in copy number and therefore also add breakpoints that should not be considered. It is therefore also advised to `filter breakpoints from SVs that do not cause a change in copynumber, such as inversions`.

## Docker - experimental

Battenberg can be run inside a Docker container. Please follow the instructions below.

#### Installation

```
git clone git@github.com:Wedge-Oxford/battenberg.git
cd battenberg
docker build -t battenberg:2.2.9 .
```

#### Reference data

First, download the Battenberg reference data from the URL provided further in this README. Then in the ```impute_info.txt``` file, replace the paths to the reference files with ```/opt/battenberg_reference```. I.e. the path to the first legend file should become:

```
/opt/battenberg_reference/1000genomes_2012_v3_impute/ALL_1000G_phase1integrated_v3_chr1_impute.legend
```

#### Run interactively

These commands run the Battenberg pipeline within a Docker container in interactive mode. This command assumes the data is available locally in `$PWD/data/pcawg/HCC1143_ds` and the reference files have been placed in `$PWD/battenberg_reference`

```
docker run -it -v `pwd`/data/pcawg/HCC1143_ds:/mnt/battenberg/ -v `pwd`/battenberg_reference:/opt/battenberg_reference battenberg:2.2.8
```

Within the Docker terminal run the pipeline, in this case on the ICGC PCAWG testing data available [here](https://s3-eu-west-1.amazonaws.com/wtsi-pancancer/testdata/HCC1143_ds.tar).

```
R CMD BATCH '--no-restore-data --no-save --args -t HCC1143 -n HCC1143_BL --nb /mnt/battenberg/HCC1143_BL.bam --tb /mnt/battenberg/HCC1143.bam --sex female -o /mnt/battenberg/' /usr/local/bin/battenberg_wgs.R /mnt/battenberg/battenberg.Rout
```
  
### Building a release

In RStudio: In the Build tab, click Check Package

Then open the ```NAMESPACE``` file and edit:

```
S3method(plot,haplotype.data)
```  

to:

```
export(plot.haplotype.data)
```



##### hg38 for Beagle5

Modified original code to derive the input vcf for Beagle5 and hg38:

```
#!/bin/bash
#
# READ_ME file (08 Dec 2015)
# 
# 1000 Genomes Project Phase 3 data release (version 5a) in VCF format for use with Beagle version 4.x
#    Data Source: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/*
#
# NOTES:
# 
# 1) Markers with <5 copies of the reference allele or <5 copies of the non-reference alleles have been excluded.
#
# 2) Structural variants have been excluded.
#
# 3) All non-unique identifiers in the ID column are removed from the ID column
#
# 4) Additional marker filtering may be performed using the gtstats.jar and filterlines.jar utilities
#
# 5) Sample information is in files: 
#      integrated_call_samples.20130502.ALL.ped
#      integrated_call_samples_v3.20130502.ALL.panel
#
# 6) Male haploid chromosome X genotypes are encoded as diploid homozygous genotypes.
#
############################################################################
#  The following shell script was used to create the files in this folder  #
############################################################################
#

## required if loading modules
module load Java
module load HTSlib

## wget for GRCh37
## wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/*

# wget for GRCh38 (liftover from hg38)
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/*
# see article: https://wellcomeopenresearch.org/articles/4-50
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr*
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/utilities/filterlines.jar
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/utilities/gtstats.jar
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/utilities/remove.ids.jar
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/utilities/simplify-vcf.jar

## Downloadable from Beagle5 web page
gts="java -ea -jar gtstats.jar"
fl="java -ea -jar filterlines.jar"
rmids="java -ea -jar remove.ids.jar"
simplify="java -ea -jar simplify-vcf.jar"
min_minor="5"

## Running directory 
src="./"
#mkdir ${src}
#cd ${src}
#cd -

## Go through autosomes and prepare vcf
for chr in $(seq 4 5); do
echo "chr${chr}"
#input="${src}ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz"
input="${src}ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
#vcf_removefield="${input}_removedfield.vcg.gz"
vcf="chr${chr}.1kg.phase3.v5a_GRCh38.vcf.gz"
excl="chr${chr}.1kg.phase3.v5a_GRCh38.excl"

#zcat ${input} | awk '/^[^#]/ {gsub(";GRC.*","",$9);print}' > ${vcf_removefield}
zcat ${input} | grep -v 'SVTYPE' | ${gts} | ${fl} \# -13 ${min_minor} 1000000 | cut -f2 > ${excl}
zcat ${input} | grep -v 'SVTYPE' | grep -v '^#' | cut -f3 | tr ';' '\n' | sort | uniq -d > chr${chr}.dup.id

# BEGIN: add 4 duplicate markers to exclusion list
#if [ ${chr} == "8" ]; then echo "88316919"; fi >> ${excl}
#if [ ${chr} == "12" ]; then echo ""; fi >> ${excl}
#if [ ${chr} == "14" ]; then echo "21181798"; fi  >> ${excl}
#if [ ${chr} == "17" ]; then echo "1241338"; fi  >> ${excl}
# END:  add 4 duplicate markers to exclusion list

zcat ${input} | grep -v 'SVTYPE' | ${fl} \# \-2 ${excl} | ${simplify} | ${rmids} chr${chr}.dup.id | bgzip -c > ${vcf}
tabix ${vcf}
done

## Same for chromosome X                                                                                                  
chr="X"
#in="${src}ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz"
in="${src}ALL.chr${chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
vcf="chr${chr}.1kg.phase3.v5a_GRCh38.vcf.gz"
excl="chr${chr}.1kg.phase3.v5a_GRCh38.excl"

### sed command recodes male chrom X genotypes as homozygote diploid genotypes
### sed command makes temporary change to floating point numbers and diploid genotypes to permit use of word-boundary '\>'
cat_in=" zcat ${in} | grep -v 'SVTYPE' | sed \
-e 's/\t0\t\.\t/\t111\tPASS\t/' \
-e 's/0\tPASS/222\tPASS/' \
-e 's/\([0-9]\)\./\1WXYZ/g' \
-e 's/\([0-9]\)|\([0-9]\)/\1X\2/g' \
-e 's/\t\([0-9]\)\>/\t\1|\1/g' \
-e 's/\([0-9]\)WXYZ/\1./g' \
-e 's/\([0-9]\)X\([0-9]\)/\1|\2/g';"

echo ${chr}
  eval ${cat_in} | grep -v '^#' | cut -f3 | tr ';' '\n' | sort | uniq -d > chr${chr}.dup.id
  eval ${cat_in} | ${gts} | ${fl} '\#' -13 ${min_minor} 1000000 | cut -f2 > ${excl}

  # BEGIN: add duplicate markers to exclusion list
  #if [ ${chr} == "X" ]; then echo "5457254"; fi >> ${excl}
  #if [ ${chr} == "X" ]; then echo "32344545"; fi >> ${excl}
  #if [ ${chr} == "X" ]; then echo "68984437"; fi >> ${excl}
  # END:  add duplicate markers to exclusion list

eval ${cat_in} | ${fl} \# \-2 ${excl} | ${simplify} | ${rmids} chr${chr}.dup.id | bgzip -c > ${vcf}
tabix ${vcf}
```

Run R code to generate loci, allele and gc_content files:

```
##########################################################################
## set working directory to where the vcf files are located
setwd("./")
##########################################################################
ffs <- dir(full=T,pattern="1kg.phase3.v5a_GRCh38.vcf.gz$")
MC.CORES <- 10 ## set number of cores to use
##########################################################################


##########################################################################
library(parallel)
##########################################################################
## Further remove unreferenced alleles and duplicates
##########################################################################
mclapply(ffs[length(ffs):1],function(x)
{
    cat(paste(x,"read file"))
    out  <- gsub(".vcf.gz","nounref.vcf.gz",x)
    cmd <- paste0("zcat ",
                  x,
                  " | grep -v '",
                  "\\",
                  ".",
                  "|","' ",
                  " | grep -v '",
                  "|",
                  "\\",
                  ".",
                  "' "
                 ,"|"," awk '!seen[$2]++' |  gzip > ",
                  out)
    system(cmd)
},mc.cores=MC.CORES)
##########################################################################

ffs <- dir(full=T,pattern="1kg.phase3.v5a_GRCh38nounref.vcf.gz$")
## Generate Loci files
##########################################################################
mclapply(ffs[length(ffs):1],function(x)
{
    cat(paste(x,"read file"))
    out  <- gsub(".vcf.gz","_loci.txt",x)
    cmd <- paste0("zcat ",
                  x," | grep -v '^#' | awk -v OFS='\t' '{print $1, $2}' > ",
                  out)
    system(cmd)
},mc.cores=MC.CORES)
##########################################################################
## with chr string (for BAMs that contain "chr")
mclapply(ffs[length(ffs):1],function(x)
{
    ##cat(paste(x,"read file"))
    out  <- gsub(".vcf.gz","_loci_chrstring.txt",x)
    cmd <- paste0("zcat ",
                  x," | grep -v '^#' | awk -v OFS='\t' '{print ",
                  '"chr"'
                 ,"$1, $2}' > ",
                  out)
    system(noquote(cmd))
},mc.cores=MC.CORES)
##########################################################################
## Generate Allele Files
##########################################################################
mclapply(ffs[length(ffs):1],function(x)
{
    cat(paste(x,"read file"))
    out  <- gsub(".vcf.gz","_allele_letter.txt",x)
    cmd <- paste0("zcat ",
                  x," | grep -v '^#' | awk -v OFS='\t' '{print $2, $4, $5}' > ",
                  out)
    system(cmd)
},mc.cores=MC.CORES)
##########################################################################
## Convert Alleles into Indices for BB input
indices <- c("A"=1,"C"=2,"G"=3,"T"=4)
##########################################################################
mclapply(ffs[length(ffs):1],function(x)
{
    cat(".")
    inp <- gsub(".vcf.gz","_allele_letter.txt",x)
    out <- gsub(".vcf.gz","_allele_index.txt",x)
    df <- as.data.frame(data.table::fread(inp))
    ref <- indices[df[,2]]
    alt <- indices[df[,3]]
    ndf <- data.frame(position=df[,1],
                      a0=ref,
                      a1=alt)
    write.table(ndf,file=out,
                row.names=F,col.names=T,sep="\t",quote=F)
},mc.cores=5)
##########################################################################


##########################################################################
## Symlink loci to change the names allowing for a "prefix" before chr in BB
##########################################################################
allfs <- dir(full=T)
allfs_loci <- allfs[grepl("loci.txt$",allfs)]
tnull <- lapply(allfs_loci,function(x)
{
    cmd <- paste0("ln -s ",x," ",gsub("chr(.*?)\\.(.*).txt","\\2_chr\\1.txt",x))
    system(cmd)
    if(grepl("chrX",x))
    {
        cmd <- paste0("ln -s ",x," ",gsub("chr(.*?)\\.(.*).txt","\\2_chr23.txt",x))
        system(cmd)
    }
})
##########################################################################
allfs <- dir(full=T)
allfs_loci <- allfs[grepl("loci_chrstring.txt$",allfs)]
tnull <- lapply(allfs_loci,function(x)
{
    cmd <- paste0("ln -s ",x," ",gsub("chr(.*?)\\.(.*).txt","\\2_chr\\1.txt",x))
    system(cmd)
    if(grepl("chrX",x))
    {
        cmd <- paste0("ln -s ",x," ",gsub("chr(.*?)\\.(.*).txt","\\2_chr23.txt",x))
        system(cmd)
    }
})
##########################################################################
## Symlink alleles: same as for loci, symlink to change name for the
## use of prefixes
##########################################################################
allfs <- dir(full=T)
allfs_index <- allfs[grepl("allele_index",allfs)]
tnull <- lapply(allfs_index,function(x)
{
    cmd <- paste0("ln -s ",x," ",gsub("chr(.*?)\\.(.*).txt","\\2_chr\\1.txt",x))
    system(cmd)
    if(grepl("chrX",x))
    {
        cmd <- paste0("ln -s ",x," ",gsub("chr(.*?)\\.(.*).txt","\\2_chr23.txt",x))
        system(cmd)
    }
})
       
##########################################################################
##  Derive GC content files
##########################################################################
library(Rsamtools)
library(data.table)
library(Biostrings)
##########################################################################
gcTrack <- function(chr,
                    pos,
                    dna,
                    window=5000)
{
    gc <- rowSums(letterFrequencyInSlidingView(dna[[chr]],
                                               window,
                                               c("G","C")))/window
    gc[pos]
}
getRefGenome <- function (fasta = FASTA, CHRS = paste0("", c(1:22, "X", "Y", 
    "MT"))) 
{
    dna <- Biostrings::readDNAStringSet(fasta, format = "fasta")
    dna <- lapply(1:length(CHRS), function(x) dna[[x]])
    names(dna) <- CHRS
    return(dna)
}
##########################################################################
FASTA <- "genome.fa" ## Link to genome reference fasta file 
CHRS <- paste0("", c(1:22, "X"))
BEAGLELOCI.template <- "chrCHROMNAME.1kg.phase3.v5a_GRCh38nounref_loci.txt"
##########################################################################
REFGENOME <- getRefGenome(fasta = FASTA, CHRS = CHRS) ## Loads genome
in memory
##########################################################################
OUTDIR <- "1000genomes_2012_v3_gcContent_hg38"
system(paste0("mkdir ",OUTDIR))
##########################################################################

##########################################################################
windows <- c(25,
             50,
             100,
             200,
             500,
             1000,
             2000,
             5000,
             10000,
             20000,
             50000,
             100000,
             200000,
             500000,
             1000000,
             2000000,
             5000000,
             10000000)
names(windows) <- formatC(windows,format="f",digits=0)
names(windows) <- gsub("000000$","Mb",names(windows))
names(windows) <- gsub("000$","kb",names(windows))
names(windows) <- sapply(names(windows),function(x) if(grepl("[0-9]$",x)) paste0(x,"bp") else x)
##########################################################################

writeGC <- function(gccontent,chr,outdir)
{
    write.table(gccontent,
                file=gzfile(paste0(outdir,"/1000_genomes_GC_corr_chr_",chr,".txt.gz")),
                col.names=T,
                row.names=T,quote=F,sep="\t")
}

mclapply(CHRS,function(chr)
{
	cat(chr)
	snppos <- as.data.frame(data.table::fread(gsub("CHROMNAME",chr,BEAGLELOCI.template)))
    gccontent <- sapply(windows,function(window) gcTrack(chr=chr,
                                                         pos=snppos[,2],
                                                         dna=REFGENOME,
                                                         window=window))*100
    gccontent <- cbind(rep(chr,nrow(gccontent)),snppos[,2],gccontent)
    colnames(gccontent)[1:2] <- c("chr","Position")
    rownames(gccontent) <- snppos[,2]
    writeGC(gccontent,chr,OUTDIR)
},mc.cores=10)
   
```

###### Example run

To run using Beagle5, simply parametrise the same way you would run
under impute2. It should be back compatible, so you can run impute2
by setting usebeagle=FALSE. And it uses the same input files needed for the
pipeline, i.e. 1000G loci/alleles + ref panel + prob loci + imputeinfo file etc.

The map plink files for Beagle can be downloaded from:
http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/


```
BEAGLEJAR <- "$PATHTOBEAGLEFILES/beagle.24Aug19.3e8.jar"
BEAGLEREF.template <- "$PATHTOBEAGLEFILES/chrCHROMNAME.1kg.phase3.v5a.b37.bref3"
BEAGLEPLINK.template <- "$PATHTOBEAGLEFILES/plink.chrCHROMNAME.GRCh37.map"

timed <- system.time(battenberg(tumourname=TUMOURNAME,
                                normalname=NORMALNAME,
                                tumour_data_file=TUMOURBAM,
                                normal_data_file=NORMALBAM,
                                imputeinfofile=IMPUTEINFOFILE,
                                g1000prefix=G1000PREFIX,
                                problemloci=PROBLEMLOCI,
                                gccorrectprefix=GCCORRECTPREFIX,
                                repliccorrectprefix=REPLICCORRECTPREFIX,
                                g1000allelesprefix=G1000PREFIX_AC,
                                ismale=IS_MALE,
                                data_type="wgs",
                                impute_exe="impute2",
                                allelecounter_exe="alleleCounter",
                                nthreads=NTHREADS,
                                platform_gamma=1,
                                phasing_gamma=1,
                                segmentation_gamma=10,
                                segmentation_kmin=3,
                                phasing_kmin=1,
                                clonality_dist_metric=0,
                                ascat_dist_metric=1,
                                min_ploidy=1.6,
                                max_ploidy=4.8, min_rho=0.1,
                                min_goodness=0.63,
                                uninformative_BAF_threshold=0.51,
                                min_normal_depth=10,
                                min_base_qual=20,
                                min_map_qual=35,
                                calc_seg_baf_option=1,
                                skip_allele_counting=F,
                                skip_preprocessing=F,
                                skip_phasing=F,
                                usebeagle=USEBEAGLE, ##set to TRUE to use beagle
                                beaglejar=BEAGLEJAR, ##path
                                beagleref=BEAGLEREF.template, ##pathtemplate
                                beagleplink=BEAGLEPLINK.template, ##pathtemplate
                                beaglemaxmem=15, 
                                beaglenthreads=1,
                                beaglewindow=40,
                                beagleoverlap=4,
                                snp6_reference_info_file=NA,
                                apt.probeset.genotype.exe="apt-probeset-genotype",
                                apt.probeset.summarize.exe="apt-probeset-summarize",
                                norm.geno.clust.exe="normalize_affy_geno_cluster.pl",
                                birdseed_report_file="birdseed.report.txt",
                                heterozygousFilter="none",
                                prior_breakpoints_file=NULL))
```

