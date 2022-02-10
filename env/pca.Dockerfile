FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

RUN conda install bioconda::multiqc=1.7
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  
