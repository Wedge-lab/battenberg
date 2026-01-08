FROM ubuntu:24.04

USER root
ARG DEBIAN_FRONTEND=noninteractive

# 1. Install System Essentials + R + Java
RUN apt-get update && apt-get install -y \
    r-base r-base-dev \
    openjdk-8-jdk \
    libcurl4-gnutls-dev libxml2-dev libssl-dev \
    make git curl \
    && rm -rf /var/lib/apt/lists/*

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