# ============================================================================
# converter.R — Core conversion: embeddings, cells, gene expression → parquets
#
# Outputs files matching PrecisionDashWidgets' dataManager.js expectations.
# ============================================================================

library(arrow)
library(data.table)
library(Matrix)

# =============================================================================
# EMBEDDINGS
# =============================================================================

#' Write embedding parquet files
#'
#' Each embedding gets its own file with standardized columns:
#'   cell_id, umap_1, umap_2
#'
#' Returns a list of embedding info for the manifest.
write_embeddings <- function(seurat_obj, config, is_multiome) {
  reductions <- names(seurat_obj@reductions)
  log_info("Available reductions: ", paste(reductions, collapse = ", "))

  embedding_info <- list(
    files = list(),
    default_embedding = NULL
  )

  # Define which embeddings to extract and their output filenames
  embedding_map <- list()

  if (is_multiome) {
    # Multiome: extract WNN UMAP, ATAC UMAP, RNA UMAP
    # WNN UMAP
    wnn_key <- NULL
    for (k in c("wnn.umap", "wnn_umap", "wnnumap")) {
      if (k %in% reductions) { wnn_key <- k; break }
    }
    if (!is.null(wnn_key)) {
      embedding_map[["wnn_umap_complete.parquet"]] <- wnn_key
    }

    # ATAC UMAP
    atac_key <- NULL
    for (k in c("atac.umap", "atac_umap", "atacumap")) {
      if (k %in% reductions) { atac_key <- k; break }
    }
    if (!is.null(atac_key)) {
      embedding_map[["atac_umap_complete.parquet"]] <- atac_key
    }

    # RNA UMAP — may need to compute from PCA
    rna_key <- NULL
    for (k in c("rna.umap", "rna_umap", "rnaumap", "umap")) {
      if (k %in% reductions) { rna_key <- k; break }
    }

    if (is.null(rna_key)) {
      # Try to compute RNA UMAP from PCA — returns modified seurat_obj
      compute_result <- try_compute_rna_umap(seurat_obj)
      if (!is.null(compute_result)) {
        seurat_obj <- compute_result$seurat_obj
        rna_key <- compute_result$key
      }
    }

    if (!is.null(rna_key)) {
      embedding_map[["rna_umap_complete.parquet"]] <- rna_key
    }

    # Default embedding: WNN > ATAC > RNA > first available
    default_key <- if (!is.null(wnn_key)) wnn_key else if (!is.null(atac_key)) atac_key else rna_key
  } else {
    # scRNA: extract UMAP and/or tSNE
    if ("umap" %in% reductions) {
      embedding_map[["umap_complete.parquet"]] <- "umap"
    }
    if ("tsne" %in% reductions) {
      embedding_map[["tsne_complete.parquet"]] <- "tsne"
    }
    default_key <- if ("umap" %in% reductions) "umap" else if ("tsne" %in% reductions) "tsne" else NULL
  }

  # Write each embedding
  for (filename in names(embedding_map)) {
    reduction_key <- embedding_map[[filename]]
    log_info("Extracting embedding '", reduction_key, "' -> ", filename)

    embeddings <- Embeddings(seurat_obj, reduction = reduction_key)

    # Take first 2 dimensions only
    if (ncol(embeddings) < 2) {
      log_warn("Skipping embedding '", reduction_key, "': fewer than 2 dimensions")
      next
    }

    # Standardize column names
    embed_df <- data.frame(
      cell_id = rownames(embeddings),
      umap_1 = as.numeric(embeddings[, 1]),
      umap_2 = as.numeric(embeddings[, 2]),
      stringsAsFactors = FALSE
    )

    filepath <- file.path(config$output_dir, filename)
    write_parquet_safe(embed_df, filepath)
    embedding_info$files[[filename]] <- reduction_key
  }

  # Write default embedding as umap_complete.parquet (for multiome)
  if (is_multiome && !is.null(default_key)) {
    log_info("Writing default embedding (", default_key, ") as umap_complete.parquet")
    embeddings <- Embeddings(seurat_obj, reduction = default_key)
    embed_df <- data.frame(
      cell_id = rownames(embeddings),
      umap_1 = as.numeric(embeddings[, 1]),
      umap_2 = as.numeric(embeddings[, 2]),
      stringsAsFactors = FALSE
    )
    filepath <- file.path(config$output_dir, "umap_complete.parquet")
    write_parquet_safe(embed_df, filepath)
    embedding_info$files[["umap_complete.parquet"]] <- default_key
    embedding_info$default_embedding <- default_key
  } else if (!is_multiome) {
    embedding_info$default_embedding <- default_key
  }

  log_info("Wrote ", length(embedding_info$files), " embedding files")
  return(embedding_info)
}

#' Try to compute RNA UMAP from PCA if it's available
try_compute_rna_umap <- function(seurat_obj) {
  pca_key <- NULL
  for (k in c("pca", "rna.pca", "rna_pca")) {
    if (k %in% names(seurat_obj@reductions)) { pca_key <- k; break }
  }

  if (is.null(pca_key)) {
    log_warn("No PCA found — cannot compute RNA UMAP")
    return(NULL)
  }

  log_info("Computing RNA UMAP from PCA reduction '", pca_key, "'")

  tryCatch({
    # Set RNA as default assay for UMAP computation
    DefaultAssay(seurat_obj) <- "RNA"

    seurat_obj <- RunUMAP(seurat_obj,
                          reduction = pca_key,
                          dims = 1:min(30, ncol(Embeddings(seurat_obj, pca_key))),
                          reduction.name = "rna_umap_computed",
                          reduction.key = "rnaumap_",
                          verbose = FALSE)

    log_info("RNA UMAP computed successfully")
    return(list(seurat_obj = seurat_obj, key = "rna_umap_computed"))
  }, error = function(e) {
    log_warn("Failed to compute RNA UMAP: ", e$message)
    return(NULL)
  })
}

# =============================================================================
# CELL METADATA
# =============================================================================

#' Write cells.parquet with all metadata columns
#'
#' Factors are converted to strings. Returns the detected cell type column name.
write_cells <- function(seurat_obj, config) {
  metadata <- seurat_obj@meta.data

  # Create cells dataframe with cell_id as first column
  cells_df <- data.frame(
    cell_id = rownames(metadata),
    stringsAsFactors = FALSE
  )

  # Add all metadata columns, converting factors to strings
  for (col_name in names(metadata)) {
    val <- metadata[[col_name]]
    if (is.factor(val)) {
      cells_df[[col_name]] <- as.character(val)
    } else {
      cells_df[[col_name]] <- val
    }
  }

  filepath <- file.path(config$output_dir, "cells.parquet")
  write_parquet_safe(cells_df, filepath)

  # Detect and report cell type column
  cell_type_col <- get_cell_type_column(metadata)
  if (!is.null(cell_type_col)) {
    n_types <- length(unique(metadata[[cell_type_col]]))
    log_info("Cell type column: '", cell_type_col, "' (", n_types, " types)")
  } else {
    log_warn("No cell type column detected")
  }

  return(cell_type_col)
}

# =============================================================================
# GENE EXPRESSION
# =============================================================================

#' Write all gene expression files:
#'   - genes.parquet (gene_name, gene_id)
#'   - gene_stats.parquet (gene_name, mean_expr, total_expr, pct_cells, n_cells_expressing)
#'   - gene_locations.parquet (gene_name, location, is_precomputed)
#'   - chunks/chunk_NNNNN.parquet for all genes (cell_id, expression, gene_id)
#'
#' Returns gene_info list for manifest.
write_gene_expression <- function(seurat_obj, config) {
  # Use RNA assay
  DefaultAssay(seurat_obj) <- "RNA"

  # Get the expression matrix (genes x cells, sparse)
  expr_matrix <- GetAssayData(seurat_obj, layer = "data")
  gene_names <- rownames(expr_matrix)
  cell_ids <- colnames(expr_matrix)
  n_genes <- length(gene_names)
  n_cells <- length(cell_ids)

  log_info("Expression matrix: ", n_genes, " genes x ", n_cells, " cells")

  # --- genes.parquet ---
  genes_df <- data.frame(
    gene_name = gene_names,
    gene_id = seq(0, n_genes - 1),  # 0-indexed
    stringsAsFactors = FALSE
  )
  write_parquet_safe(genes_df, file.path(config$output_dir, "genes.parquet"))

  # --- gene_stats.parquet ---
  log_info("Computing gene statistics...")
  gene_stats <- compute_gene_stats(expr_matrix, gene_names, n_cells)
  write_parquet_safe(gene_stats, file.path(config$output_dir, "gene_stats.parquet"))
  clean_memory("gene stats")

  # --- Create output directory ---
  chunks_dir <- file.path(config$output_dir, "chunks")
  dir.create(chunks_dir, recursive = TRUE, showWarnings = FALSE)

  # --- gene_locations tracking ---
  gene_locations <- data.frame(
    gene_name = character(0),
    location = character(0),
    is_precomputed = logical(0),
    stringsAsFactors = FALSE
  )

  # --- Write all genes into chunk files ---
  log_info("Writing ", n_genes, " genes in chunks of ", config$chunk_size, "...")
  chunk_idx <- 0

  for (start in seq(1, n_genes, by = config$chunk_size)) {
    end <- min(start + config$chunk_size - 1, n_genes)
    batch_indices <- start:end
    chunk_idx <- chunk_idx + 1

    chunk_filename <- sprintf("chunks/chunk_%05d.parquet", chunk_idx)

    # Build chunk data
    chunk_rows <- list()
    for (gene_idx in batch_indices) {
      expr_col <- expr_matrix[gene_idx, ]
      expressing <- which(expr_col != 0)

      if (length(expressing) > 0) {
        chunk_rows[[length(chunk_rows) + 1]] <- data.frame(
          cell_id = cell_ids[expressing],
          expression = as.numeric(expr_col[expressing]),
          gene_id = rep(gene_idx - 1L, length(expressing)),  # 0-indexed
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(chunk_rows) > 0) {
      chunk_df <- do.call(rbind, chunk_rows)
    } else {
      chunk_df <- data.frame(
        cell_id = character(0),
        expression = numeric(0),
        gene_id = integer(0),
        stringsAsFactors = FALSE
      )
    }

    write_parquet_safe(chunk_df, file.path(config$output_dir, chunk_filename))

    # Add gene locations for this chunk
    for (gene_idx in batch_indices) {
      gene_locations <- rbind(gene_locations, data.frame(
        gene_name = gene_names[gene_idx],
        location = chunk_filename,
        is_precomputed = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    if (chunk_idx %% 10 == 0) {
      log_info("  Written ", chunk_idx, " chunk files")
      clean_memory(paste0("chunk batch ", chunk_idx))
    }
  }

  log_info("Wrote ", chunk_idx, " chunk files total")

  # --- gene_locations.parquet ---
  write_parquet_safe(gene_locations, file.path(config$output_dir, "gene_locations.parquet"))

  # Explicitly free the expression matrix
  rm(expr_matrix)
  clean_memory("gene expression complete")

  gene_info <- list(
    n_genes = n_genes,
    n_chunks = chunk_idx
  )

  return(gene_info)
}

#' Compute per-gene statistics from sparse expression matrix
compute_gene_stats <- function(expr_matrix, gene_names, n_cells) {
  # Use Matrix sparse operations for efficiency
  # Row sums = total expression per gene
  total_expr <- Matrix::rowSums(expr_matrix)

  # Number of cells expressing each gene (non-zero entries per row)
  n_expressing <- Matrix::rowSums(expr_matrix != 0)

  # Mean expression = total / n_cells
  mean_expr <- total_expr / n_cells

  # Percent cells expressing
  pct_cells <- (n_expressing / n_cells) * 100

  data.frame(
    gene_name = gene_names,
    mean_expr = as.numeric(mean_expr),
    total_expr = as.numeric(total_expr),
    pct_cells = as.numeric(pct_cells),
    n_cells_expressing = as.integer(n_expressing),
    stringsAsFactors = FALSE
  )
}
