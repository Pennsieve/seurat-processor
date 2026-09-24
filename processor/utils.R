# ============================================================================
# utils.R — Shared utility functions for seurat-processor
# ============================================================================

library(arrow)

# --- Logging -----------------------------------------------------------------

log_info <- function(...) {
  message(paste0("[INFO] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " — ", ...))
}

log_warn <- function(...) {
  message(paste0("[WARN] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " — ", ...))
}

log_error <- function(...) {
  message(paste0("[ERROR] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " — ", ...))
}

# --- Parquet helpers ---------------------------------------------------------

#' Write a data.frame to parquet with snappy compression
write_parquet_safe <- function(df, filepath) {
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(df, filepath, compression = "snappy")
  log_info("Wrote ", filepath, " (", nrow(df), " rows, ", ncol(df), " cols)")
}

# --- Detection helpers -------------------------------------------------------

#' Check if a Seurat object contains a ChromatinAssay (multiome indicator)
has_chromatin_assay <- function(seurat_obj) {
  for (assay_name in names(seurat_obj@assays)) {
    assay <- seurat_obj@assays[[assay_name]]
    if (inherits(assay, "ChromatinAssay")) {
      return(TRUE)
    }
  }
  return(FALSE)
}

#' Get the name of the ChromatinAssay (prefer ATAC_macs2, then ATAC, then first found)
get_chromatin_assay_name <- function(seurat_obj) {
  chromatin_assays <- c()
  for (assay_name in names(seurat_obj@assays)) {
    assay <- seurat_obj@assays[[assay_name]]
    if (inherits(assay, "ChromatinAssay")) {
      chromatin_assays <- c(chromatin_assays, assay_name)
    }
  }

  if (length(chromatin_assays) == 0) return(NULL)

  # Prefer ATAC_macs2, then ATAC, then first found
  if ("ATAC_macs2" %in% chromatin_assays) return("ATAC_macs2")
  if ("ATAC" %in% chromatin_assays) return("ATAC")
  return(chromatin_assays[1])
}

#' Find the best cell type column in metadata
#' Preference order: final_cluster_ids > cell_type > seurat_clusters > first factor column
get_cell_type_column <- function(metadata) {
  preferred <- c("final_cluster_ids", "cell_type", "celltype", "CellType",
                 "cell_type_annotation", "seurat_clusters")

  for (col in preferred) {
    if (col %in% names(metadata)) {
      return(col)
    }
  }

  # Fallback: first factor/character column that looks like cell types
  for (col in names(metadata)) {
    if (is.factor(metadata[[col]]) || is.character(metadata[[col]])) {
      n_unique <- length(unique(metadata[[col]]))
      if (n_unique > 1 && n_unique < 200) {
        return(col)
      }
    }
  }

  return(NULL)
}

#' Get the default embedding name for the object
#' For multiome: prefer wnn.umap, then umap
#' For scRNA: prefer umap, then tsne
get_default_embedding <- function(seurat_obj, is_multiome) {
  reductions <- names(seurat_obj@reductions)

  if (is_multiome) {
    if ("wnn.umap" %in% reductions) return("wnn.umap")
    if ("wnn_umap" %in% reductions) return("wnn_umap")
  }

  if ("umap" %in% reductions) return("umap")
  if ("tsne" %in% reductions) return("tsne")

  # Return first available
  if (length(reductions) > 0) return(reductions[1])
  return(NULL)
}

#' Memory-friendly garbage collection with logging
clean_memory <- function(label = "") {
  gc_result <- gc(verbose = FALSE)
  used_mb <- sum(gc_result[, 2])
  log_info("GC after ", label, " — ", round(used_mb, 1), " MB in use")
}
