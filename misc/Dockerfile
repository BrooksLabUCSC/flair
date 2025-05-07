# 24.04 LTS (Noble Numbat)
FROM ubuntu:24.04
LABEL maintainer="Mark Diekhans <markd@ucsc.edu>"

# Controls if FLAIR will be installed, don't install for testing before publishing
ARG FLAIR_INSTALL=yes

# deal with tzdata timezone question
# see https://dev.to/setevoy/docker-configure-tzdata-and-timezone-during-build-20bk
ENV TZ=UTC
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Install Essentials
RUN apt-get update && apt-get install -y \
	build-essential \
	wget \
        gcc \
        make \
	git \
	pkg-config \
        libc-bin \
        libc6-dev \
	bedtools=2.31.1+dfsg-2 \
        minimap2=2.26+dfsg-1build1 \
        python3=3.12.3-0ubuntu2 \
	python3-pip=24.0+dfsg-1ubuntu1 \
        r-base

#  24.04.2 Noble has an old samtool package, download newer from a subsequent release
RUN \
   wget -nv -O /tmp/samtools.deb http://mirrors.kernel.org/ubuntu/pool/universe/s/samtools/samtools_1.21-1_amd64.deb && \
   wget -nv -O /tmp/libhts3t64.deb http://mirrors.kernel.org/ubuntu/pool/universe/h/htslib/libhts3t64_1.21+ds-1_amd64.deb && \
   apt-get install -y /tmp/libhts3t64.deb /tmp/samtools.deb

RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); \
  update.packages(ask = FALSE); \
  install.packages(c("devtools", "BiocManager", "ggplot2", "qqman", "lazyeval")); \
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); \
  BiocManager::install(c("DRIMSeq", "stageR", "DESeq2", "apeglm"), ask = FALSE)'

RUN if [ "$FLAIR_INSTALL" = "yes" ]; then \
        pip install --break-system-packages flair-brookslab[diffexp]==v2.2.0 ; \
    fi

RUN mkdir /data
WORKDIR /data

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* 

