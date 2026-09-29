# Helper: install and verify, stop the build if anything fails
install_and_check <- function(pkgs, ...) {
  install.packages(pkgs, ...)
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("FATAL: Package '", pkg, "' failed to install"), call. = FALSE)
    }
  }
}

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

install_and_check(cran_pkgs, repos = "https://cloud.r-project.org")

# Seurat and SeuratObject from CRAN
install_and_check(c('SeuratObject', 'Seurat'), repos = "https://cloud.r-project.org")

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

for (pkg in c('GenomicRanges', 'GenomeInfoDb', 'IRanges', 'Rsamtools',
              'Biostrings', 'BSgenome', 'EnsDb.Hsapiens.v86', 'biovizBase')) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("FATAL: Bioconductor package '", pkg, "' failed to install"), call. = FALSE)
  }
}

# Signac for chromatin assay support
install_and_check('Signac', repos = "https://cloud.r-project.org")

message("All packages installed successfully")
