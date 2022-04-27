FROM ubuntu:18.04

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/expression-counting

ENV TZ=Asia/Seoul
ENV DEBCONF_NOWARNINGS yes
ARG DEBIAN_FRONTEND=noninteractive

RUN mkdir tools

<<<<<<< HEAD
RUN apt-get update && apt-get install -y tzdata build-essential wget python2.7 python-pip python3.7 python3-pip tar cmake 
RUN pip install pandas 
RUN pip3 install pandas
RUN pip3 install argparse

RUN apt-get install -y zlib1g-dev
=======
RUN apt-get update && apt-get install -y tzdata build-essential wget python3.7 python3-pip tar cmake
RUN apt install zlib1g-dev 
>>>>>>> bf86503fcf918cc7607314d6312f37048fa55f26

RUN cd /tools &&\
    wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.7.tar.gz &&\
    tar -xzvf stringtie-2.1.7.tar.gz &&\
    cd stringtie-2.1.7 && \
    make release

<<<<<<< HEAD
ENV PATH="/tools:${PATH}:/tools/stringtie-2.1.7:/usr/bin/python3:/usr/bin/python"
=======
RUN pip3 install pandas && \
    pip3 install argparse


ENV PATH="/tools:${PATH}:/tools/stringtie-2.1.7"
>>>>>>> bf86503fcf918cc7607314d6312f37048fa55f26

COPY *.py /tools/

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
