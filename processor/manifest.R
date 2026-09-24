# ============================================================================
# manifest.R — Generate manifest.json for PrecisionDashWidgets
# ============================================================================

library(jsonlite)

#' Write manifest.json describing the output dataset
write_manifest <- function(config, embedding_info, gene_info, chromatin_info,
                           is_multiome, cell_type_col, n_cells) {

  # Build modalities list
  modalities <- list("rna")
  if (is_multiome) {
    modalities <- list("rna", "atac")
  }

  # Build embeddings section
  embeddings <- list()
  for (filename in names(embedding_info$files)) {
    reduction_key <- embedding_info$files[[filename]]
    embeddings[[length(embeddings) + 1]] <- list(
      key = reduction_key,
      file = filename,
      columns = list("cell_id", "umap_1", "umap_2")
    )
  }

  # Build files section
  files <- list(
    umap = "umap_complete.parquet",
    cells = "cells.parquet",
    genes = "genes.parquet",
    gene_stats = "gene_stats.parquet",
    gene_locations = "gene_locations.parquet"
  )

  # Add optional embedding files
  if ("wnn_umap_complete.parquet" %in% names(embedding_info$files)) {
    files$wnn_umap <- "wnn_umap_complete.parquet"
  }
  if ("atac_umap_complete.parquet" %in% names(embedding_info$files)) {
    files$atac_umap <- "atac_umap_complete.parquet"
  }
  if ("rna_umap_complete.parquet" %in% names(embedding_info$files)) {
    files$rna_umap <- "rna_umap_complete.parquet"
  }
  if ("tsne_complete.parquet" %in% names(embedding_info$files)) {
    files$tsne <- "tsne_complete.parquet"
  }

  # Add chromatin files if present
  if (!is.null(chromatin_info$has_tracks) && chromatin_info$has_tracks) {
    files$chromatin_tracks <- "chromatin_tracks.parquet"
  }
  if (!is.null(chromatin_info$has_gene_coords) && chromatin_info$has_gene_coords) {
    files$gene_coordinates <- "gene_coordinates.parquet"
  }
  if (!is.null(chromatin_info$has_linked_peaks) && chromatin_info$has_linked_peaks) {
    files$linked_peaks <- "linked_peaks.parquet"
  }

  manifest <- list(
    version = "1.0",
    processor = "seurat-processor",
    output_mode = "expanded",
    modalities = modalities,
    n_cells = n_cells,
    cell_type_column = if (!is.null(cell_type_col)) cell_type_col else NULL,
    default_embedding = embedding_info$default_embedding,
    embeddings = embeddings,
    files = files,
    gene_expression = list(
      n_genes = gene_info$n_genes,
      n_chunks = gene_info$n_chunks,
      chunks_dir = "chunks"
    ),
    chromatin = list(
      has_tracks = isTRUE(chromatin_info$has_tracks),
      has_gene_coordinates = isTRUE(chromatin_info$has_gene_coords),
      has_linked_peaks = isTRUE(chromatin_info$has_linked_peaks)
    )
  )

  filepath <- file.path(config$output_dir, "manifest.json")
  write(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), filepath)
  log_info("Wrote ", filepath)
}
