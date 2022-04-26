FROM pegi3s/fastqc:0.11.9

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/expression-counting

RUN apt-get update && apt-get install -y python3-pip

RUN pip3 install pandas && \
    pip3 install argparse


ENV PATH="/tools:${PATH}"

COPY *.py /tools/

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
