#
# This image contains working version of custom scripts from NIH/NCI/CCBR
# https://github.com/CCBR/Pipeliner
# 
# get_strandness.py - os,sys 
# https://github.com/stevetsa/Pipeliner/blob/master/Results-template/Scripts/get_strandness.py
# merge_rsem_results.py -  os,sys,pandas 
# https://github.com/stevetsa/Pipeliner/blob/master/Results-template/Scripts/merge_rsem_results.py
# EBSeq.py - os,sys,pandas,numpy
# https://github.com/stevetsa/Pipeliner/blob/master/Results-template/Scripts/EBSeq.py


FROM ubuntu:16.04
MAINTAINER Steve Tsang <mylagimail2004@yahoo.com>
RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
 build-essential \
 apt-utils \
 git-all \
 python \ 
 python-pip \
 wget

COPY requirements.txt /opt/.
RUN pip install --upgrade pip
RUN pip install -r /opt/requirements.txt

RUN DEBIAN_FRONTEND=noninteractive apt-get install --yes \
  python3 \
  python3-pip && pip3 install snakemake

WORKDIR /opt
RUN git clone https://github.com/stevetsa/Pipeliner.git
RUN cp /opt/Pipeliner/Results-template/Scripts/*.py /opt/.
#RUN cp /opt/*.py /usr/local/bin
