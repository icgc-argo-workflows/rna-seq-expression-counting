FROM ubuntu:18.04

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/expression-counting

ENV TZ=Asia/Seoul
ENV DEBCONF_NOWARNINGS yes
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt-get update && apt-get install -y software-properties-common 
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt install -y  python3.7 python3-pip
RUN python3.7 -m pip install pip
RUN apt-get install -y tzdata build-essential wget tar cmake
RUN apt-get install -y build-essential python3-numpy python3-matplotlib 
RUN apt-get install -y python3-distutils python3-setuptools 

RUN python3.7 -m pip install pip --upgrade pip

RUN pip3 install pysam
RUN pip3 install HTSeq==2.0.0
RUN pip3 install pandas 
RUN pip3 install argparse


ENV PATH="/tools:${PATH}"

COPY *.py /tools/

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
