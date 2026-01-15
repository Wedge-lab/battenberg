FROM ubuntu:24.04

USER root
ARG DEBIAN_FRONTEND=noninteractive

# 1. Install System Essentials + R + Java
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    openjdk-17-jdk \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libxml2 \
    libssl-dev \
    make \
    git \ 
    curl \
    dirmngr \
    software-properties-common \
    r-cran-rgl \
    && rm -rf /var/lib/apt/lists/*

ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH=$JAVA_HOME/bin:$PATH

## install alll the side tools: 

RUN mkdir /tmp/downloads

RUN curl -sSL -o tmp.tar.gz --retry 10 https://github.com/samtools/htslib/archive/1.7.tar.gz && \
    mkdir /tmp/downloads/htslib && \
    tar -C /tmp/downloads/htslib --strip-components 1 -zxf tmp.tar.gz && \
    make -C /tmp/downloads/htslib && \
    rm -f /tmp/downloads/tmp.tar.gz

ENV HTSLIB=/tmp/downloads/htslib

RUN curl -sSL -o tmp.tar.gz --retry 10 https://github.com/cancerit/alleleCount/archive/v4.0.0.tar.gz && \
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



# 2. Install pak (the engine for your Makefile)
RUN Rscript -e "install.packages('pak', repos = 'https://cloud.r-project.org')"

# 3. Setup work directory and copy the project
WORKDIR /opt/battenberg
COPY . .

# 4. Use the Makefile to do the heavy lifting
# This installs R deps via pak and then installs the package itself
RUN make deps
RUN make install

WORKDIR /home/ubuntu
CMD ["/bin/bash"]