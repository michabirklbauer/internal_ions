# Dockerfile for Internal Ion Explorer
# author: Micha Birklbauer
# version: 1.2.1

FROM ubuntu:22.04

LABEL maintainer="micha.birklbauer@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y gnupg gnupg1 gnupg2 software-properties-common
RUN apt-get update && apt-get install -y \
    curl \
    git \
    python3-distutils \
    python3-lxml \
    python3-pip \
    libglib2.0-0 \
    libsm6 \
    libxrender1 \
    libxext6

RUN mkdir internal_ions
COPY ./ internal_ions/
WORKDIR internal_ions

RUN pip3 install -r requirements.txt

CMD  ["streamlit", "run", "streamlit_app.py"]
