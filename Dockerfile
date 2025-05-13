# Dockerfile for Internal Ion Explorer
# author: Micha Birklbauer
# version: 1.2.2

FROM python:3.12

LABEL maintainer="micha.birklbauer@gmail.com"

RUN mkdir internal_ions
COPY ./ internal_ions/
WORKDIR internal_ions

RUN python3 -m venv venv && \
    . venv/bin/activate && \
    pip install --upgrade pip && \
    pip install --upgrade setuptools && \
    pip install --no-cache-dir -r env.txt

CMD  ["sh", "-c", ". venv/bin/activate && streamlit run streamlit_app.py"]
