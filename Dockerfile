FROM rocker/r-ver:4.1.2

RUN apt-get update

# for R install
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev zlib1g-dev curl

# for Rhtslib
RUN apt-get install -y libbz2-dev liblzma-dev 

RUN R -e "install.packages(c('devtools', 'testthat', 'roxygen2', 'BiocManager'))" 

# Battenberg prerequisites
RUN R -q -e 'BiocManager::install(c("devtools", "splines", "readr", "doParallel", "ggplot2", "RColorBrewer", "gridExtra", "gtools", "parallel", "igordot/copynumber", "VariantAnnotation"))'
RUN R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")'

# Battenberg
RUN R -q -e 'devtools::install_github("Wedge-Oxford/battenberg")'
