# ============================================================================
# legacy.R — Original processor logic (OUTPUT_MODE=consolidated)
#
# Preserved verbatim from the original main.R.
# Produces a single results.parquet + viewer_config.yml.
# ============================================================================

run_legacy <- function(seurat_obj, config) {
  library(nanoparquet)
  library(jsonlite)
  library(SeuratObject)
  library(yaml)

  log_info("Running in CONSOLIDATED (legacy) mode")

  data <- seurat_obj

  variables = c()
  for (r in names(data@reductions)) {
    if (r == 'umap' || r == 'tsne') {
      dimNames <- names(data@reductions[[r]])
      variables <- c(variables, dimNames)
    }
  }
  variables <- c(variables, names(data@meta.data))

  result <- FetchData(data, variables)

  metadata <- c()
  for (n in names(data@meta.data)) {
    result[[n]] <- as.factor(result[[n]])
    metadata <- c(metadata, setNames(c(toJSON(unique(result[[n]]))), c(n) ))
  }

  outputFolder = config$output_dir

  # Viewer Config File
  viewerInfo <- list(name = "parquet-umap-viewer", options = list())
  write_yaml(viewerInfo, file.path(outputFolder, "viewer_config.yml"), fileEncoding = "UTF-8")

  # Write Parquet File
  write_parquet(result, file.path(outputFolder, "results.parquet"),
   meta=metadata, compression= "snappy")

  log_info("Legacy mode complete: results.parquet + viewer_config.yml")
}
