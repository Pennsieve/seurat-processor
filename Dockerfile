FROM rocker/r-ver:4.4.2

WORKDIR /app

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
    libgsl-dev \
    libfftw3-dev \
    libudunits2-dev \
    libgdal-dev \
    libproj-dev \
    libsqlite3-dev \
    cmake \
    pkg-config \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R --version

# Install R packages (copy requirements first for Docker layer caching)
COPY processor/requirements.R /app/processor/requirements.R
RUN Rscript /app/processor/requirements.R

# Copy processor code
COPY processor/ /app/processor/

# Create data directories
RUN mkdir -p /data/input /data/output

# Copy and set entrypoint
COPY entrypoint.sh /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh

ENTRYPOINT ["/app/entrypoint.sh"]
