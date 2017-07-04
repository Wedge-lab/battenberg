# Battenberg
-----

## Advised installation and running instructions

Please visit the [cgpBattenberg page](https://github.com/cancerit/cgpBattenberg), download the code from there and run the setup.sh script. The battenberg.pl script can then be used to run the pipeline


## Advanced installation instructions

At the moment Battenberg can only be installed directly from Github. The instructions below will install the latest development version, which might not always work out of the box.

### Prerequisites

Installing from Github requires devtools and Battenberg requires readr, RColorBrewer and ASCAT. The pipeline requires doParallel. From the command line run:

  > R -q -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("devtools", "readr", "doParallel", "ggplot2", "RColorBrewer"));'
  > R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")'

### Installation from Github

To install Battenberg, run the following from the command line:

  > R -q -e 'devtools::install_github("sdentro/battenberg")'

### Required reference files

Battenberg requires a number of reference files that should be downloaded.

  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_impute_1000G_v1.tar.gz
  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_snp6_exe.tgz (SNP6 only)
  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_snp6_ref.tgz (SNP6 only)
  * GC content file to be added
  * ignore loci file to be added
  
### Pipeline

Go into inst/example for example WGS and SNP6 R-only pipelines.
  
  