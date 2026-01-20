# Stage 1: Build C dependencies
FROM ubuntu:24.04 AS builder
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    make git curl gcc g++ bzip2 zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /tmp/downloads
# Build htslib
RUN curl -sSL -o htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && \
    mkdir /tmp/htslib && \
    tar -C /tmp/htslib --strip-components 1 -xjf htslib.tar.bz2 && \
    cd /tmp/htslib && \
    ./configure && \
    make -j$(nproc) && \
    make install

# Build alleleCount
RUN curl -sSL -o allelecount.tar.gz https://github.com/cancerit/alleleCount/archive/v4.0.0.tar.gz && \
    mkdir /tmp/allelecount && \
    tar -C /tmp/allelecount --strip-components 1 -zxf allelecount.tar.gz && \
    cd /tmp/allelecount/c && \
    mkdir -p bin && \
    make bin/alleleCounter && \
    cp bin/alleleCounter /usr/local/bin/


# Stage 2: Final image
FROM ubuntu:24.04
ARG DEBIAN_FRONTEND=noninteractive

# 1. Install R and System Dependencies
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    openjdk-17-jdk \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    make \
    curl \
    git \
    r-cran-rgl \
    && rm -rf /var/lib/apt/lists/*

# 2. OPTIMIZATION: Configure Posit Binary Repository for Ubuntu Noble
# We do this AFTER R is installed so the directory exists.
RUN mkdir -p /usr/lib/R/etc && \
    echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))' >> /usr/lib/R/etc/Rprofile.site && \
    echo 'options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))' >> /usr/lib/R/etc/Rprofile.site

# 3. Copy binaries from builder stage
COPY --from=builder /usr/local/bin/alleleCounter /usr/local/bin/
# Impute2 (Static x86_64 binary)
RUN curl -sSL -o tmp.tar.gz https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz && \
    tar -C /usr/local/bin --strip-components 1 -zxf tmp.tar.gz && \
    rm tmp.tar.gz

# 4. Install pak (improved installation for Linux)
RUN Rscript -e "install.packages('pak', repos = 'https://r-lib.github.io/p/pak/stable')"

WORKDIR /opt/battenberg

# 5. OPTIMIZATION: Cache dependency installation layer
# Copy DESCRIPTION first so that changes to code don't invalidate the dependency cache.
COPY DESCRIPTION .
COPY Makefile .
RUN make deps

# 6. Copy the rest of the code and install the package
COPY . .
RUN make install

WORKDIR /home/ubuntu
CMD ["/bin/bash"]