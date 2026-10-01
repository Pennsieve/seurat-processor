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
library(Signac)

# --- Main --------------------------------------------------------------------

main <- function() {
  total_timer <- timer_start()
  log_info("=== Seurat Processor Starting ===")

  # System info
  log_system_info()

  # Load configuration
  config <- load_config()
  print_config(config)

  # Check disk space
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  log_disk_space(config$output_dir)

  # Find input RDS file
  input_files <- list.files(path = config$input_dir, pattern = "\\.rds$",
                            ignore.case = TRUE, full.names = TRUE)
  if (length(input_files) == 0) {
    log_error("No .rds files found in ", config$input_dir)
    log_error("Directory contents: ", paste(list.files(config$input_dir), collapse = ", "))
    stop(paste0("No .rds files found in ", config$input_dir))
  }

  input_file <- input_files[1]
  log_input_file(input_file)

  # Load Seurat object
  log_info("Loading Seurat object (this may take several minutes for large files)...")
  load_timer <- timer_start()
  seurat_obj <- tryCatch({
    readRDS(input_file)
  }, error = function(e) {
    log_error("Failed to load RDS file: ", e$message)
    log_error("This often means the file is corrupt or the system ran out of memory.")
    stop(e)
  })
  timer_log(load_timer, "RDS load")

  # Validate it's a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    log_error("Loaded object is not a Seurat object. Class: ", paste(class(seurat_obj), collapse = ", "))
    stop("Input file does not contain a Seurat object")
  }

  log_info("Loaded object of class: ", paste(class(seurat_obj), collapse = ", "))
  log_info("Cells: ", ncol(seurat_obj))
  log_info("Assays: ", paste(names(seurat_obj@assays), collapse = ", "))
  log_info("Reductions: ", paste(names(seurat_obj@reductions), collapse = ", "))
  log_info("Metadata columns: ", paste(names(seurat_obj@meta.data), collapse = ", "))
  clean_memory("RDS load")

  # Drop unused assays to free memory — keep RNA, chromatin assay, and chromvar
  # See README for details on which assays are retained and why
  keep_assays <- c("RNA", "chromvar")
  chromatin_name <- get_chromatin_assay_name(seurat_obj)
  if (!is.null(chromatin_name)) keep_assays <- c(keep_assays, chromatin_name)
  keep_assays <- intersect(keep_assays, names(seurat_obj@assays))  # only keep what exists
  drop_assays <- setdiff(names(seurat_obj@assays), keep_assays)
  if (length(drop_assays) > 0) {
    log_info("Dropping unused assays to free memory: ", paste(drop_assays, collapse = ", "))
    for (a in drop_assays) {
      seurat_obj[[a]] <- NULL
    }
    clean_memory("drop unused assays")
  }

  # Dispatch based on output mode
  if (config$output_mode == "consolidated") {
    run_legacy(seurat_obj, config)
  } else {
    run_expanded(seurat_obj, config)
  }

  # Output summary
  log_output_summary(config$output_dir)
  timer_log(total_timer, "=== Seurat Processor")
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
  t1 <- timer_start()
  embedding_info <- write_embeddings(seurat_obj, config, is_multiome)
  timer_log(t1, "Phase 1 (Embeddings)")
  # Drop reductions to free memory — no longer needed
  for (r in names(seurat_obj@reductions)) {
    seurat_obj@reductions[[r]] <- NULL
  }
  clean_memory("embeddings + drop reductions")

  # Phase 2: Cell metadata
  log_info("--- Phase 2: Cell metadata ---")
  t2 <- timer_start()
  n_cells <- ncol(seurat_obj)  # save before we start stripping
  cell_type_col <- write_cells(seurat_obj, config)
  timer_log(t2, "Phase 2 (Cell metadata)")
  clean_memory("cells")

  # Phase 3: Gene expression
  log_info("--- Phase 3: Gene expression ---")
  t3 <- timer_start()
  gene_info <- tryCatch({
    write_gene_expression(seurat_obj, config)
  }, error = function(e) {
    log_error("Gene expression extraction failed: ", e$message)
    log_error("This may indicate insufficient memory or an unsupported assay format.")
    stop(e)
  })
  timer_log(t3, "Phase 3 (Gene expression)")
  # Drop RNA assay to free memory before chromatin phase
  if ("RNA" %in% names(seurat_obj@assays)) {
    seurat_obj[["RNA"]] <- NULL
  }
  clean_memory("gene expression + drop RNA")

  # Phase 4: Chromatin (multiome only)
  chromatin_info <- list()
  if (is_multiome && !is.null(chromatin_assay)) {
    log_info("--- Phase 4: Chromatin ---")
    t4 <- timer_start()
    chromatin_info <- write_chromatin_data(seurat_obj, config, chromatin_assay)
    timer_log(t4, "Phase 4 (Chromatin)")
    clean_memory("chromatin")
  } else {
    log_info("--- Phase 4: Chromatin — skipped (not multiome) ---")
  }

  # Phase 5: Manifest + viewer config
  log_info("--- Phase 5: Manifest ---")
  write_manifest(config, embedding_info, gene_info, chromatin_info,
                 is_multiome, cell_type_col, n_cells)

  # Viewer config
  viewerInfo <- list(name = "parquet-umap-viewer", options = list())
  yaml::write_yaml(viewerInfo, file.path(config$output_dir, "viewer_config.yml"),
                   fileEncoding = "UTF-8")

  log_info("Expanded mode complete")
}

# Run with top-level error handling
tryCatch({
  main()
}, error = function(e) {
  log_error("=== PROCESSOR FAILED ===")
  log_error("Error: ", e$message)
  log_error("Call: ", paste(deparse(e$call), collapse = " "))
  clean_memory("at failure")
  quit(status = 1)
})
