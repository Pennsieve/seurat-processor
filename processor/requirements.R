# Use Posit Package Manager for pre-built CRAN Linux binaries
ppm <- "https://p3m.dev/cran/__linux__/noble/latest"

# CRAN packages (from PPM binaries)
install.packages(c(
  'arrow', 'nanoparquet', 'data.table', 'jsonlite', 'yaml',
  'Matrix', 'Rcpp', 'Seurat', 'SeuratObject'
), repos = ppm)

# BiocManager (from PPM)
install.packages("BiocManager", repos = ppm)

# Bioconductor packages — let BiocManager use its OWN repos (not PPM)
# PPM doesn't mirror Bioconductor, so we must not override repos here
BiocManager::install(c(
  'GenomicRanges', 'GenomeInfoDb', 'IRanges', 'Rsamtools',
  'Biostrings', 'BSgenome', 'EnsDb.Hsapiens.v86', 'biovizBase'
), ask = FALSE, update = FALSE)

# Signac — depends on Bioconductor pkgs above
install.packages('Signac', repos = ppm)

# --- Verify ---
critical <- c('arrow', 'Seurat', 'SeuratObject', 'Signac',
              'GenomicRanges', 'GenomeInfoDb', 'data.table', 'jsonlite')
for (pkg in critical) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  message(if (ok) "OK: " else "MISSING: ", pkg)
}
