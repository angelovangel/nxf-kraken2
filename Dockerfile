# build for nextflow-kraken2, based on miniconda3 as multiqc 1.9 needs python3
# build stage
FROM continuumio/miniconda3:4.7.12

LABEL description="Docker image containing all requirements for running kraken2, bracken and kaiju"
LABEL maintainer="Angel Angelov <aangeloo@gmail.com>"

# the base image of continuumio does not have some programs needed to compile kaiju, procps needed ny nextflow
RUN apt-get update --fix-missing && \
    apt-get install -y git gcc g++ libz-dev make procps && apt-get clean -y

# kaiju is not available on anaconda
# install bracken also from github
RUN git clone https://github.com/bioinformatics-centre/kaiju.git && \
    cd kaiju/src && \
    make
RUN git clone https://github.com/jenniferlu717/Bracken.git && \
    cd Bracken && \
    chmod a+x install_bracken.sh && ./install_bracken.sh

# make kaiju and bracken binaries available
ENV PATH /kaiju/bin:/Bracken:$PATH

# final stage
COPY environment.yml .
RUN conda env update -n root -f environment.yml && \
    conda clean -afy