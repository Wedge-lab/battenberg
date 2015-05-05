# Battenberg
-----

## Installation instructions

At the moment Battenberg can only be installed directly from Github. The instructions below will install the latest development version, which might not always work out of the box.

### Prerequisites

Installing from Github requires devtools and Battenberg requires doParallel. From the command line run:

  > R -q -e 'source("http://bioconductor.org/biocLite.R"); biocLite("devtools"); biocLite("doParallel")'

### Installation from Github

To install Battenberg, run the following from the command line:

  > R -q -e 'devtools::install_github("sdentro/battenberg", auth_token="797b8d19e041bdc6404d0f605e427537e44db9b6")'

### Required reference files

Battenberg requires a number of reference files that should be downloaded.

  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_snp6_exe.tgz
  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_snp6_ref.tgz
  * ftp://ftp.sanger.ac.uk/pub/teams/113/Battenberg/battenberg_impute_1000G_v1.tar.gz