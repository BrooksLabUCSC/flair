# 26.04 LTS (resolute)
FROM ubuntu:26.04 AS base
ENV DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Mark Diekhans <markd@ucsc.edu>"

###
# Arguments
##

# Controls if FLAIR will be installed, set to `no' to not install
# install for testing before publishing
ARG FLAIR_INSTALL=yes

# pip specification from where to install FLAIR, which can be a URL
#    ARG FLAIR_PIP_SPEC=git+https://github.com/BrooksLabUCSC/flair.git
ARG FLAIR_PIP_SPEC=flair-brookslab==v3.0.1

##
# safer bash
## 
SHELL ["/bin/bash", "-beEu", "-o", "pipefail", "-c"]

##
# deal with tzdata timezone question
# see https://dev.to/setevoy/docker-configure-tzdata-and-timezone-during-build-20bk
##
ENV TZ=UTC
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

##
# Install essentials
##
RUN apt-get update && apt-get install -y \
	build-essential \
	wget curl gcc make \
	git pkg-config \
        libc-bin libc6-dev \
        file net-tools iproute2 iputils-ping dnsutils \
	bedtools=2.31.1+dfsg-3 \
	samtools=1.22.1-1 \
        python3=3.14.3-0ubuntu2 \
	python3-pip \
        libtirpc-dev \
        libcurl4-openssl-dev \
        zlib1g-dev libbz2-dev liblzma-dev \
        r-base

# minimap2 2.30
RUN \
    MM2_VER=2.30 && MM2_LINUX=minimap2-${MM2_VER}_x64-linux && \
    MM2_URL=https://github.com/lh3/minimap2/releases/download/v${MM2_VER}/${MM2_LINUX}.tar.bz2 && \
    curl -sL ${MM2_URL} | tar -xj -C /tmp && \
    cp  /tmp/${MM2_LINUX}/{minimap2,k8,paftools.js} /usr/local/bin/ && \
    rm -rf /tmp/${MM2_LINUX}

FROM base AS unix_base

##
# R dependencies
##
RUN Rscript -e 'options(repos=c(CRAN = "https://cloud.r-project.org"), \
                        Ncpus=parallel::detectCores()); \
  update.packages(ask=FALSE); \
  install.packages(c("devtools", "BiocManager", "ggplot2", "qqman", "lazyeval", "argparse", "data.table")); \
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); \
  BiocManager::install(c("DRIMSeq", "stageR", "DESeq2", "apeglm"), ask=FALSE)'

FROM unix_base AS r_base

##
# FLAIR
##
RUN if [ "$FLAIR_INSTALL" = "yes" ]; then \
        pip install --break-system-packages ${FLAIR_PIP_SPEC} ; \
    fi

RUN mkdir /data
WORKDIR /data

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* 

