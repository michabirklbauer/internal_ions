# Dockerfile for Internal Ion Explorer
# author: Micha Birklbauer
# version: 1.2.1

FROM ubuntu:23.10

LABEL maintainer="micha.birklbauer@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y gnupg gnupg1 gnupg2 software-properties-common
RUN apt-get update && apt-get install -y \
    curl \
    git \
    python3-lxml \
    python3-pip \
    python3-venv \
    python3-setuptools \
    libglib2.0-0 \
    libsm6 \
    libxrender1 \
    libxext6

RUN mkdir internal_ions
COPY ./ internal_ions/
WORKDIR internal_ions

RUN python3 -m venv venv && \
    . venv/bin/activate && \
    pip install --upgrade pip && \
    pip install --upgrade setuptools && \
    pip install -r python3117.txt

CMD  ["sh", "-c", ". venv/bin/activate && streamlit run streamlit_app.py"]
