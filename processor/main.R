# ============================================================================
# main.R — Seurat Processor Orchestrator
#
# Dispatches to either expanded mode (converter.R) or consolidated mode (legacy.R)
# based on the OUTPUT_MODE environment variable.
# ============================================================================

# Load modules (resolve paths whether run from repo root or processor dir)
if (file.exists("utils.R")) {
  script_dir <- "."
} else if (file.exists("processor/utils.R")) {
  script_dir <- "processor"
} else {
  stop("Cannot find processor R files. Run from repo root or processor/ directory.")
}

source(file.path(script_dir, "utils.R"))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "legacy.R"))
source(file.path(script_dir, "converter.R"))
source(file.path(script_dir, "chromatin.R"))
source(file.path(script_dir, "manifest.R"))

# Libraries
library(SeuratObject)
library(Seurat)

# --- Main --------------------------------------------------------------------

main <- function() {
  log_info("=== Seurat Processor Starting ===")

  # Load configuration
  config <- load_config()
  print_config(config)

  # Find input RDS file
  input_files <- list.files(path = config$input_dir, pattern = "\\.rds$",
                            ignore.case = TRUE, full.names = TRUE)
  if (length(input_files) == 0) {
    stop(paste0("No .rds files found in ", config$input_dir))
  }

  input_file <- input_files[1]
  log_info("Loading Seurat object from: ", input_file)

  # Load Seurat object
  seurat_obj <- readRDS(input_file)
  log_info("Loaded object of class: ", paste(class(seurat_obj), collapse = ", "))
  log_info("Cells: ", ncol(seurat_obj), " | Assays: ", paste(names(seurat_obj@assays), collapse = ", "))
  log_info("Reductions: ", paste(names(seurat_obj@reductions), collapse = ", "))
  clean_memory("RDS load")

  # Create output directory
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  # Dispatch based on output mode
  if (config$output_mode == "consolidated") {
    run_legacy(seurat_obj, config)
  } else {
    run_expanded(seurat_obj, config)
  }

  log_info("=== Seurat Processor Complete ===")
}

# --- Expanded mode orchestrator ----------------------------------------------

run_expanded <- function(seurat_obj, config) {
  log_info("Running in EXPANDED mode")

  # Detect multiome
  is_multiome <- has_chromatin_assay(seurat_obj)
  log_info("Multiome detected: ", is_multiome)

  # Determine chromatin assay name if multiome
  chromatin_assay <- NULL
  if (is_multiome) {
    chromatin_assay <- get_chromatin_assay_name(seurat_obj)
    log_info("Chromatin assay: ", chromatin_assay)
  }

  # Phase 1: Embeddings
  log_info("--- Phase 1: Embeddings ---")
  embedding_info <- write_embeddings(seurat_obj, config, is_multiome)
  clean_memory("embeddings")

  # Phase 2: Cell metadata
  log_info("--- Phase 2: Cell metadata ---")
  cell_type_col <- write_cells(seurat_obj, config)
  clean_memory("cells")

  # Phase 3: Gene expression
  log_info("--- Phase 3: Gene expression ---")
  gene_info <- write_gene_expression(seurat_obj, config)
  clean_memory("gene expression")

  # Phase 4: Chromatin (multiome only)
  chromatin_info <- list()
  if (is_multiome && !is.null(chromatin_assay)) {
    log_info("--- Phase 4: Chromatin ---")
    chromatin_info <- write_chromatin_data(seurat_obj, config, chromatin_assay)
    clean_memory("chromatin")
  }

  # Phase 5: Manifest + viewer config
  log_info("--- Phase 5: Manifest ---")
  write_manifest(config, embedding_info, gene_info, chromatin_info,
                 is_multiome, cell_type_col, ncol(seurat_obj))

  # Viewer config
  viewerInfo <- list(name = "parquet-umap-viewer", options = list())
  yaml::write_yaml(viewerInfo, file.path(config$output_dir, "viewer_config.yml"),
                   fileEncoding = "UTF-8")

  log_info("Expanded mode complete")
}

# Run
main()
