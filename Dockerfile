# multistage build for nextflow-kraken2, first stage is for building kaiju only
# build stage
FROM debian:latest AS build

# the base debian image of continuumio does not have some programs needed to compile kaiju
RUN apt-get update --fix-missing && \
    apt-get install -y git gcc g++ libz-dev make && apt-get clean -y

# kaiju is not available on anaconda
RUN git clone https://github.com/bioinformatics-centre/kaiju.git && \
    cd kaiju/src && \
    make

# final stage
FROM continuumio/miniconda:4.7.12
LABEL description="Docker image containing all requirements for running kraken2, bracken and kaiju"
LABEL maintainer="Angel Angelov <aangeloo@gmail.com>"

# make kaiju binaries available
COPY --from=build /kaiju/bin /kaiju/bin
ENV PATH /kaiju/bin:$PATH

COPY environment.yml .
RUN conda env update -n root -f environment.yml && \
    conda clean -afy
# procps needed by nextflow
RUN apt-get update && apt-get install -y procps && \
    apt-get clean -y