# Use Posit Package Manager for pre-built CRAN Linux binaries
ppm <- "https://p3m.dev/cran/__linux__/noble/latest"

# CRAN packages (from PPM binaries)
install.packages(c(
  'arrow', 'nanoparquet', 'data.table', 'jsonlite', 'yaml',
  'Matrix', 'Rcpp', 'Seurat', 'SeuratObject'
), repos = ppm)

# BiocManager (from PPM)
install.packages("BiocManager", repos = ppm)

# Install ALL Bioconductor + Signac through BiocManager
# BiocManager resolves both Bioconductor and CRAN deps together
BiocManager::install(c(
  'GenomicRanges', 'GenomeInfoDb', 'IRanges', 'Rsamtools',
  'Biostrings', 'BSgenome', 'EnsDb.Hsapiens.v86', 'biovizBase',
  'Signac'
), ask = FALSE, update = FALSE)

# --- Verify ---
critical <- c('arrow', 'Seurat', 'SeuratObject', 'Signac',
              'GenomicRanges', 'GenomeInfoDb', 'data.table', 'jsonlite')
for (pkg in critical) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  message(if (ok) "OK: " else "MISSING: ", pkg)
}
