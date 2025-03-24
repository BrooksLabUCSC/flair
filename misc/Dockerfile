FROM ubuntu:22.04

MAINTAINER Mark Diekhans, markd@ucsc.edu

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
	samtools \
	bedtools \
	minimap2 \
	python3 \
	python3-pip 

RUN pip install flair-brookslab==2.0.0

RUN mkdir /data
WORKDIR /data

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* 


