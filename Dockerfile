FROM r-base

# Add dependencies
RUN apt-get update && apt-get install -y libxml2 libxml2-dev libcurl4-gnutls-dev r-cran-rgl

RUN R -q -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("devtools","RColorBrewer","ggplot2","gridExtra","readr"))'
RUN R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")'
RUN R -q -e 'devtools::install_github("Wedge-Oxford/battenberg")'

RUN wget https://raw.githubusercontent.com/Wedge-Oxford/battenberg/master/inst/example/battenberg_wgs.R
RUN wget https://raw.githubusercontent.com/Wedge-Oxford/battenberg/master/inst/example/battenberg_snp6.R

ADD battenberg_wgs.R /usr/local/bin/battenberg_wgs.R
ADD battenberg_snp6.R /usr/local/bin/battenberg_snp6.R
