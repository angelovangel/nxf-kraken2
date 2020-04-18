FROM continuumio/miniconda:4.7.12

LABEL description="Docker image containing all requirements for running kraken2 and bracken"
LABEL maintainer="Angel Angelov <aangeloo@gmail.com>"

COPY environment.yml .
RUN conda env update -n root -f environment.yml && conda clean -afy
RUN apt-get update