# Battenberg
-----

## Advised installation and running instructions

Please visit the [cgpBattenberg page](https://github.com/cancerit/cgpBattenberg), download the code from there and run the setup.sh script. The battenberg.pl script can then be used to run the pipeline.


## Advanced installation instructions

The instructions below will install the latest stable Battenberg version. Please take this approach only when you'd like to do something not covered by cgpBattenberg.

### Standalone

#### Prerequisites

Installing from Github requires devtools and Battenberg requires readr, RColorBrewer and ASCAT. The pipeline requires doParallel. From the command line run:

```
R -q -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("devtools", "readr", "doParallel", "ggplot2", "RColorBrewer", "gridExtra", "gtools"));'
R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")'
```

#### Installation from Github

To install Battenberg, run the following from the command line:

```
R -q -e 'devtools::install_github("Wedge-Oxford/battenberg")'
```

#### Required reference files

Battenberg requires a number of reference files that should be downloaded.

  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_1000genomesloci2012_v3.tar.gz
  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_impute_1000G_v3.tar.gz
  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/probloci_270415.txt.gz
  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_wgs_gc_correction_1000g_v3.tar.gz
  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_snp6_exe.tgz (SNP6 only)
  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_snp6_ref.tgz (SNP6 only)
  
#### Pipeline

Go into inst/example for example WGS and SNP6 R-only pipelines.
  
### Docker - experimental

Battenberg can be run inside a Docker container. Please follow the instructions below.

#### Installation

```
wget https://raw.githubusercontent.com/Wedge-Oxford/battenberg/dev/Dockerfile
docker build -t battenberg:2.2.8 .
```

#### Download reference data

To do

#### Run interactively

These commands run the Battenberg pipeline within a Docker container in interactive mode. This command assumes the data is available locally in `$PWD/data/pcawg/HCC1143_ds` and the reference files have been placed in `$PWD/battenberg_reference`

```
docker run -it -v `pwd`/data/pcawg/HCC1143_ds:/mnt/battenberg/ -v `pwd`/battenberg_reference:/opt/battenberg_reference battenberg:2.2.8
```

Within the Docker terminal run the pipeline, in this case on the ICGC PCAWG testing data available [here](https://s3-eu-west-1.amazonaws.com/wtsi-pancancer/testdata/HCC1143_ds.tar).

```
R CMD BATCH '--no-restore-data --no-save --args HCC1143 HCC1143_BL /mnt/battenberg/HCC1143_BL.bam /mnt/battenberg/HCC1143.bam FALSE /mnt/battenberg/' /usr/local/bin/battenberg_wgs.R /mnt/battenberg/battenberg.Rout
```
  
### Building a release

In RStudio: In the Build tab, click Check Package

Then open the NAMESPACE file and edit:

```
S3method(plot,haplotype.data)
```  

to:

```
export(plot.haplotype.data)
```
