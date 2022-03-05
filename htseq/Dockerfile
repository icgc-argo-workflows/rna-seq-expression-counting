FROM ubuntu:18.04

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/expression-counting

ENV TZ=Asia/Seoul
ENV DEBCONF_NOWARNINGS yes
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y tzdata build-essential wget python3.7 python3-pip tar cmake

RUN pip3 install HTseq &&\
    pandas && \
    pip3 install argparse


ENV PATH="/tools:${PATH}"

COPY *.py /tools/

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
