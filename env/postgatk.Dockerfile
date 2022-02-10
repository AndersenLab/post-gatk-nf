FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

COPY ugh_conda.yml .
RUN conda env update -n root -f ugh_conda.yml && conda clean -a

RUN conda install -c bioconda bedtools=2.30.0
RUN conda install -c bioconda mosdepth=0.3.2
RUN conda install -c bioconda plink=1.90b6.21
RUN conda install -c bioconda quicktree=2.5
# RUN conda install -c bioconda samtools
RUN conda install -c bioconda vcflib
# RUN conda install -c bioconda eigensoft
# RUN conda install -c bioconda bioconvert

RUN apt-get --allow-releaseinfo-change update && \
   apt-get install -y procps && \
	rm -rf /var/lib/apt/lists/*