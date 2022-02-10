FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

COPY conda.yml .
RUN conda env update -n root -f conda.yml && conda clean -a

RUN conda install bioconda::quicktree
RUN conda install bioconda::bioconvert

RUN apt-get --allow-releaseinfo-change update && \
   apt-get install -y procps && \
	rm -rf /var/lib/apt/lists/*