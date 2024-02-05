# Dockerfile for Internal Ion Explorer
# author: Micha Birklbauer
# version: 1.1.1

FROM ubuntu:22.04

LABEL maintainer="micha.birklbauer@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive
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

RUN git clone https://github.com/arthur-grimaud/fragannot.git
WORKDIR fragannot

RUN pip3 install -r requirements.txt
RUN pip3 install -r developer_requirements.txt

WORKDIR publish
RUN python3 setup.py install

WORKDIR ../gui
RUN pip3 install -r requirements.txt

CMD  ["streamlit", "run", "streamlit_app.py", "--server.maxUploadSize", "5000"]
