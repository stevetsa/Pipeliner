
# R Scripts from NCI/CCBR/Pipeliner - https://github.com/CCBR/Pipeliner
#

FROM r-base:3.5.1
#FROM rocker/r-devel:latest
RUN apt-get update \
 && apt-get build-dep libcurl4-openssl-dev \
 && apt-get install -y --no-install-recommends \
   libfftw3-dev \
   gcc \
   git-all \
   file \
   libapparmor1 \
   libcurl4-openssl-dev \
   libedit2 \
   libssl-dev \
   lsb-release \
   psmisc \
   procps \
   libxml2-dev \
   pkg-config \
   libcairo2-dev \
   libmagick++-dev \
   postgresql postgresql-contrib libpq-dev libudunits2-dev libx11-dev mesa-common-dev libglu1-mesa-dev libhdf5-dev gfortran && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y pandoc python3
COPY install.packages4.R /opt/.
RUN Rscript /opt/install.packages4.R

WORKDIR /opt
RUN git clone https://github.com/stevetsa/Pipeliner.git
WORKDIR /opt/Pipeliner
RUN cp Dockerfile.R4 /opt/.
RUN cp Results-template/Scripts/*.R /opt/.
RUN cp Results-template/Scripts/*.Rmd /opt/.
ENV PATH="$PATH:/opt"
