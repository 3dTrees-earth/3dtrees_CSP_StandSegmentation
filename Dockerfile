FROM rocker/geospatial:4.4.3

ENV DEBIAN_FRONTEND=noninteractive \
    RGL_USE_NULL=TRUE

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        libgl1-mesa-dev \
        libglu1-mesa-dev \
        libx11-dev \
        libxt-dev \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled --ncpus -1 \
        BH \
        RANN \
        RCSF \
        Rcpp \
        RcppArmadillo \
        colorspace \
        conicfit \
        data.table \
        dbscan \
        doParallel \
        foreach \
        geometry \
        igraph \
        jsonlite \
        lidR \
        magrittr \
        rgl \
        testthat

WORKDIR /opt/CspStandSegmentation
COPY . /opt/CspStandSegmentation
RUN chmod -R a+rX /opt/CspStandSegmentation \
    && R CMD INSTALL --no-multiarch --with-keep.source /opt/CspStandSegmentation

CMD ["Rscript", "/opt/CspStandSegmentation/exec/run.R", "--help"]
