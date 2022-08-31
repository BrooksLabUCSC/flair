FROM ubuntu:22.04

MAINTAINER Jeltje van Baren, jeltje.van.baren@gmail.com

# deal with tzdata timezone question
# see https://dev.to/setevoy/docker-configure-tzdata-and-timezone-during-build-20bk
ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

#Install Essentials
RUN apt-get update && apt-get install -y \
	wget \
	build-essential \
	cmake \
	git \
	pkg-config \
	libcurl4-openssl-dev \
	libssl-dev \
	libxml2-dev \
	samtools \
	bedtools \
	minimap2 \
	r-base \
	salmon \
	python3 \
	python3-pip 

RUN pip install flair-brookslab==1.6.3 \
	rpy2 \
	pandas \
	matplotlib \
	seaborn

RUN R -e 'update.packages(ask=FALSE)'
RUN R -e 'install.packages(c("devtools", "BiocManager", "ggplot2", "qqman", "lazyeval"))'
RUN R -e 'requireNamespace("BiocManager"); BiocManager::install(c("DRIMSeq", "stageR", "DESeq2", "apeglm"))' 

RUN mkdir /data
WORKDIR /data

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* 


