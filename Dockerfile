FROM continuumio/miniconda:4.7.12

LABEL description="Docker image containing all requirements for running kraken2, bracken and kaiju"
LABEL maintainer="Angel Angelov <aangeloo@gmail.com>"

COPY environment.yml .
RUN conda env update -n root -f environment.yml && conda clean -afy

# the base debian image of continuumio does not have some programs needed to compile kaiju
RUN apt-get update && apt-get install -y gcc g++ libz-dev make procps && apt-get clean -y

# kaiju is not available on anaconda
RUN git clone https://github.com/bioinformatics-centre/kaiju.git && cd kaiju/src && make

# make kaiju binaries available
ENV PATH /kaiju/bin:$PATH