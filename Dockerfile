FROM rocker/r-ver:4.4.2

WORKDIR /processor
RUN apt clean && apt-get update



# R program dependencies
#RUN apt-get install -y libudunits2-dev && apt-get install -y libgeos-dev && apt-get install -y libproj-dev && apt-get -y install libnlopt-dev && apt-get -y install pkg-config && apt-get -y install gdal-bin && apt-get install -y libgdal-dev
#RUN apt-get -y install libcurl4-openssl-dev libfontconfig1-dev libxml2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
#RUN apt-get -y install glpk-utils libglpk-dev glpk-doc
RUN R --version

# Install remaining packages from source
COPY ./requirements-src.R .
RUN Rscript requirements-src.R

## Clean up package registry
#RUN rm -rf /var/lib/apt/lists/*

# Add additional program specific dependencies below
RUN mkdir -p data

COPY ./processor /processor

ENTRYPOINT [ "Rscript", "main.R" ]