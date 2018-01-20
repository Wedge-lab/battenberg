#FROM r-base
FROM ubuntu:16.04

# Add dependencies
RUN apt-get update && apt-get install -y libxml2 libxml2-dev libcurl4-gnutls-dev r-cran-rgl git libssl-dev curl

RUN mkdir /tmp/downloads

RUN curl -sSL -o tmp.tar.gz --retry 10 https://github.com/samtools/htslib/archive/1.2.1.tar.gz && \
    mkdir /tmp/downloads/htslib && \
    tar -C /tmp/downloads/htslib --strip-components 1 -zxf tmp.tar.gz && \
    make -C /tmp/downloads/htslib && \
    rm -f /tmp/downloads/tmp.tar.gz

ENV HTSLIB /tmp/downloads/htslib

RUN curl -sSL -o tmp.tar.gz --retry 10 https://github.com/cancerit/alleleCount/archive/v2.1.2.tar.gz && \
    mkdir /tmp/downloads/alleleCount && \
    tar -C /tmp/downloads/alleleCount --strip-components 1 -zxf tmp.tar.gz && \
    cd /tmp/downloads/alleleCount/c && \
    mkdir bin && \
    make && \
    cp /tmp/downloads/alleleCount/c/bin/alleleCounter /usr/local/bin/. && \
    cd /tmp/downloads && \
    rm -rf /tmp/downloads/alleleCount /tmp/downloads/tmp.tar.gz

RUN curl -sSL -o tmp.tar.gz --retry 10 https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz && \
    mkdir /tmp/downloads/impute2 && \
    tar -C /tmp/downloads/impute2 --strip-components 1 -zxf tmp.tar.gz && \
    cp /tmp/downloads/impute2/impute2 /usr/local/bin && \
    rm -rf /tmp/downloads/impute2 /tmp/downloads/tmp.tar.gz

RUN R -q -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("devtools","RColorBrewer","ggplot2","gridExtra","readr","doParallel"))'
RUN R -q -e 'devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")'
RUN R -q -e 'devtools::install_github("Wedge-Oxford/battenberg")'

RUN curl -sSL -o battenberg_wgs.R https://raw.githubusercontent.com/Wedge-Oxford/battenberg/master/inst/example/battenberg_wgs.R
# modify paths to reference files
RUN cat battenberg_wgs.R | \
    sed 's|IMPUTEINFOFILE = \".*|IMPUTEINFOFILE = \"/opt/battenberg_reference/1000genomes_2012_v3_impute/impute_info.txt\"|' | \
    sed 's|G1000PREFIX = \".*|G1000PREFIX = \"/opt/battenberg_reference/1000genomes_2012_v3_loci/1000genomesAlleles2012_chr\"|' | \
    sed 's|G1000PREFIX_AC = \".*|G1000PREFIX_AC = \"/opt/battenberg_reference/1000genomes_2012_v3_loci/1000genomesloci2012_chr\"|' | \
    sed 's|GCCORRECTPREFIX = \".*|GCCORRECTPREFIX = \"/opt/battenberg_reference/1000genomes_2012_v3_gcContent/1000_genomes_GC_corr_chr_\"|' | \
    sed 's|PROBLEMLOCI = \".*|PROBLEMLOCI = \"/opt/battenberg_reference/battenberg_problem_loci/probloci_270415.txt.gz\"|' > /usr/local/bin/battenberg_wgs.R && rm -f battenberg_wgs.R

#RUN curl -sSL -o battenberg_snp6.R https://raw.githubusercontent.com/Wedge-Oxford/battenberg/master/inst/example/battenberg_snp6.R
#RUN cat battenberg_snp6.R | \
#    sed 's|IMPUTEINFOFILE = \".*|IMPUTEINFOFILE = \"/opt/battenberg_reference/1000genomes_2012_v3_impute/impute_info.txt\"|' | \
#    sed 's|G1000PREFIX = \".*|G1000PREFIX = \"/opt/battenberg_reference/1000genomes_2012_v3_loci/1000genomesAlleles2012_chr\"|' | \
#    sed 's|SNP6_REF_INFO_FILE = \".*|SNP6_REF_INFO_FILE = \"/opt/battenberg_reference/battenberg_snp6/snp6_ref_info_file.txt\"|' > /usr/local/bin/battenberg_snp6.R && rm -f battenberg_wgs.R
