# 24.04 LTS (Noble Numbat)
FROM ubuntu:24.04 AS base
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
ARG FLAIR_PIP_SPEC=flair-brookslab==v2.2.0+master

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
	bedtools=2.31.1+dfsg-2 \
        minimap2=2.26+dfsg-1build1 \
        python3=3.12.3-0ubuntu2 \
	python3-pip=24.0+dfsg-1ubuntu1 \
        libtirpc-dev \
        libcurl4-openssl-dev \
        r-base

#  24.04.2 Noble has an old samtool package, download newer from a subsequent release
RUN \
   wget -nv -O /tmp/samtools.deb http://mirrors.kernel.org/ubuntu/pool/universe/s/samtools/samtools_1.21-1_amd64.deb && \
   wget -nv -O /tmp/libhts3t64.deb http://mirrors.kernel.org/ubuntu/pool/universe/h/htslib/libhts3t64_1.21+ds-1_amd64.deb && \
   apt-get install -y /tmp/libhts3t64.deb /tmp/samtools.deb

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

