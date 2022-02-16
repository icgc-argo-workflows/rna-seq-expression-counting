FROM combinelab/salmon

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/expression-counting

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
