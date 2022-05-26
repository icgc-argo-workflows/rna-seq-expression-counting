FROM combinelab/salmon

LABEL org.opencontainers.image.source https://github.com/icgc-argo-rna-workflows/expression-counting

ENV TZ=Asia/Seoul

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y python3-pip

RUN pip3 install pandas && \
    pip3 install argparse

RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/" >> /etc/apt/sources.list
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN apt-get update 
RUN apt-get install -y r-base 
RUN Rscript -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager');BiocManager::install('tximport')"
RUN Rscript -e "install.packages('argparse')"


COPY *.py /tools/
COPY *.R /tools/

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
