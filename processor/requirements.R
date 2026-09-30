# Use Posit Package Manager for pre-built Linux binaries (no compilation)
options(repos = c(CRAN = "https://p3m.dev/cran/__linux__/noble/latest"))

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

install.packages(cran_pkgs)

# Seurat and SeuratObject
install.packages(c('Seurat', 'SeuratObject'))

# Bioconductor packages (for Signac / multiome support)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Set BiocManager to also use PPM for CRAN deps
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
install.packages('Signac')

# --- Verify critical packages installed ---
critical <- c('arrow', 'Seurat', 'SeuratObject', 'Signac',
              'GenomicRanges', 'data.table', 'jsonlite', 'yaml')
missing <- critical[!sapply(critical, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  message("WARNING: These packages failed to install: ", paste(missing, collapse = ", "))
  message("Attempting reinstall of missing packages...")
  for (pkg in missing) {
    tryCatch(install.packages(pkg), error = function(e) {
      message("  Failed to install ", pkg, ": ", e$message)
    })
  }
}

# Final check — log status but don't block build
for (pkg in critical) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message("OK: ", pkg)
  } else {
    message("MISSING: ", pkg)
  }
}
