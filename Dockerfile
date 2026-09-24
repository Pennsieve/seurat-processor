FROM rocker/r-ver:4.4.2

WORKDIR /processor

# System dependencies for Seurat, Signac, Arrow, and genomics packages
RUN apt clean && apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgeos-dev \
    libglpk-dev \
    cmake \
    pkg-config \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R --version

# Install R packages from source
COPY ./requirements-src.R .
RUN Rscript requirements-src.R

# Create data directories
RUN mkdir -p /data/input /data/output

# Copy processor code
COPY ./processor /processor

# Copy and set entrypoint
COPY ./entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
