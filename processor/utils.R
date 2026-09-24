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

# --- System info -------------------------------------------------------------

#' Log system memory and R version info at startup
log_system_info <- function() {
  log_info("R version: ", R.version.string)
  log_info("Platform: ", R.version$platform)

  # System memory (macOS and Linux)
  total_ram <- tryCatch({
    if (Sys.info()["sysname"] == "Darwin") {
      raw <- system("sysctl -n hw.memsize", intern = TRUE)
      paste0(round(as.numeric(raw) / 1024^3, 1), " GB")
    } else if (file.exists("/proc/meminfo")) {
      raw <- system("grep MemTotal /proc/meminfo | awk '{print $2}'", intern = TRUE)
      paste0(round(as.numeric(raw) / 1024^2, 1), " GB")
    } else {
      "unknown"
    }
  }, error = function(e) "unknown")

  available_ram <- tryCatch({
    if (Sys.info()["sysname"] == "Darwin") {
      # Page size * free pages
      vm_stat <- system("vm_stat", intern = TRUE)
      page_size <- 16384  # Apple Silicon default
      free_line <- grep("Pages free", vm_stat, value = TRUE)
      inactive_line <- grep("Pages inactive", vm_stat, value = TRUE)
      free_pages <- as.numeric(gsub("[^0-9]", "", free_line))
      inactive_pages <- as.numeric(gsub("[^0-9]", "", inactive_line))
      paste0(round((free_pages + inactive_pages) * page_size / 1024^3, 1), " GB (free+inactive)")
    } else if (file.exists("/proc/meminfo")) {
      raw <- system("grep MemAvailable /proc/meminfo | awk '{print $2}'", intern = TRUE)
      paste0(round(as.numeric(raw) / 1024^2, 1), " GB")
    } else {
      "unknown"
    }
  }, error = function(e) "unknown")

  log_info("System RAM: ", total_ram, " total, ", available_ram, " available")

  # R memory limit
  r_max_vsize <- Sys.getenv("R_MAX_VSIZE", "not set")
  log_info("R_MAX_VSIZE: ", r_max_vsize)
}

#' Log input file info
log_input_file <- function(filepath) {
  file_size_gb <- round(file.size(filepath) / 1024^3, 2)
  log_info("Input file: ", basename(filepath), " (", file_size_gb, " GB)")
}

#' Log available disk space in output directory
log_disk_space <- function(output_dir) {
  tryCatch({
    if (Sys.info()["sysname"] == "Darwin") {
      raw <- system(paste0("df -g '", output_dir, "' | tail -1 | awk '{print $4}'"), intern = TRUE)
      log_info("Disk space available: ", raw, " GB in ", output_dir)
    } else {
      raw <- system(paste0("df -BG '", output_dir, "' | tail -1 | awk '{print $4}'"), intern = TRUE)
      log_info("Disk space available: ", raw, " in ", output_dir)
    }
  }, error = function(e) {
    log_warn("Could not check disk space: ", e$message)
  })
}

#' Log a summary of all output files
log_output_summary <- function(output_dir) {
  all_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  n_files <- length(all_files)
  total_size_mb <- round(sum(file.size(all_files)) / 1024^2, 1)
  log_info("Output summary: ", n_files, " files, ", total_size_mb, " MB total")

  # List each top-level file/dir with size
  top_level <- list.files(output_dir, full.names = TRUE)
  for (f in top_level) {
    if (file.info(f)$isdir) {
      dir_files <- list.files(f, recursive = TRUE, full.names = TRUE)
      dir_size <- round(sum(file.size(dir_files)) / 1024^2, 1)
      log_info("  ", basename(f), "/ — ", length(dir_files), " files, ", dir_size, " MB")
    } else {
      f_size <- round(file.size(f) / 1024^2, 2)
      log_info("  ", basename(f), " — ", f_size, " MB")
    }
  }
}

# --- Timing ------------------------------------------------------------------

#' Start a timer, returns the start time
timer_start <- function() {
  proc.time()
}

#' Log elapsed time since timer_start
timer_log <- function(start_time, label) {
  elapsed <- proc.time() - start_time
  log_info(label, " completed in ", round(elapsed["elapsed"], 1), "s")
}

# --- Parquet helpers ---------------------------------------------------------

#' Write a data.frame to parquet with snappy compression
write_parquet_safe <- function(df, filepath) {
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(df, filepath, compression = "snappy")
  file_size_kb <- round(file.size(filepath) / 1024, 1)
  log_info("Wrote ", filepath, " (", nrow(df), " rows, ", ncol(df), " cols, ", file_size_kb, " KB)")
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
