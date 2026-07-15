FROM rocker/geospatial:4.4.3

ENV DEBIAN_FRONTEND=noninteractive \
    RGL_USE_NULL=TRUE

ARG RLAS_COMMIT=82cbba42f158d1dfc91efda3207923260a052564

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        libgl1-mesa-dev \
        libglu1-mesa-dev \
        libx11-dev \
        libxt-dev \
        curl \
        patch \
        python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --break-system-packages --no-cache-dir --no-deps \
        laspy==2.6.1 \
        lazrs==0.8.1

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

COPY patches/rlas-read-all-extrabytes.patch /tmp/rlas-read-all-extrabytes.patch
RUN curl -L -sS "https://github.com/r-lidar/rlas/archive/${RLAS_COMMIT}.tar.gz" -o /tmp/rlas.tar.gz \
    && mkdir -p /tmp/rlas-src \
    && tar -xzf /tmp/rlas.tar.gz -C /tmp/rlas-src --strip-components=1 \
    && patch -d /tmp/rlas-src -p1 < /tmp/rlas-read-all-extrabytes.patch \
    && R CMD INSTALL /tmp/rlas-src \
    && rm -rf /tmp/rlas-src /tmp/rlas.tar.gz /tmp/rlas-read-all-extrabytes.patch

WORKDIR /opt/CspStandSegmentation
COPY . /opt/CspStandSegmentation
RUN chmod -R a+rX /opt/CspStandSegmentation \
    && R CMD INSTALL --no-multiarch --with-keep.source /opt/CspStandSegmentation

CMD ["Rscript", "/opt/CspStandSegmentation/exec/run.R", "--help"]
