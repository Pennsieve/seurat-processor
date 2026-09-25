# CRAN packages
cran_pkgs <- c(
  'arrow',        # Fast parquet I/O for expanded mode
  'nanoparquet',  # Lightweight parquet I/O for legacy mode
  'data.table',   # Fast data manipulation
  'jsonlite',     # JSON serialization
  'yaml',         # YAML file I/O
  'Matrix',       # Sparse matrix support (dependency of Seurat)
  'Rcpp'          # C++ interface (dependency of many packages)
)

install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

# Seurat and SeuratObject from CRAN
install.packages(c('Seurat', 'SeuratObject'), repos = "https://cloud.r-project.org")

# Bioconductor packages (for Signac / multiome support)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c(
  'GenomicRanges',
  'GenomeInfoDb',
  'IRanges',
  'Rsamtools',
  'Biostrings',
  'BSgenome',
  'EnsDb.Hsapiens.v86',
  'biovizBase'
), ask = FALSE, update = FALSE)

# Signac for chromatin assay support
install.packages('Signac', repos = "https://cloud.r-project.org")
