# Dockerfile for Internal Ion Explorer
# author: Micha Birklbauer
# version: 1.2.0

FROM ubuntu:22.04

LABEL maintainer="micha.birklbauer@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y gnupg gnupg1 gnupg2 software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get update && apt-get install -y \
    curl \
    git \
    python3-distutils \
    python3-lxml \
    python3-pip \
    libglib2.0-0 \
    libsm6 \
    libxrender1 \
    libxext6 \
    r-base \
    python3-rpy2

ENV R_HOME="/usr/lib/R"

RUN mkdir internal_ions
COPY ./ internal_ions/
WORKDIR internal_ions

RUN pip3 install -r Docker.config

CMD  ["streamlit", "run", "streamlit_app.py"]
