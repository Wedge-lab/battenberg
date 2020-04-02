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
