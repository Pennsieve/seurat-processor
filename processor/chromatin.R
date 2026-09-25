# ============================================================================
# chromatin.R — Multiome-specific: pseudobulk tracks, gene coords, peak-gene links
#
# All operations are wrapped in tryCatch — best-effort extraction.
# If Signac is not available or data is missing, these gracefully skip.
# ============================================================================

#' Write chromatin-related parquet files for multiome data
#'
#' Returns a list of chromatin info for the manifest.
write_chromatin_data <- function(seurat_obj, config, chromatin_assay) {
  chromatin_info <- list(
    has_tracks = FALSE,
    has_gene_coords = FALSE,
    has_linked_peaks = FALSE
  )

  # Try to load Signac
  signac_available <- tryCatch({
    library(Signac)
    TRUE
  }, error = function(e) {
    log_warn("Signac not available — skipping chromatin extraction: ", e$message)
    FALSE
  })

  if (!signac_available) return(chromatin_info)

  # --- Chromatin tracks (pseudobulk from peak matrix) ---
  tryCatch({
    chromatin_info$has_tracks <- write_chromatin_tracks(
      seurat_obj, config, chromatin_assay
    )
  }, error = function(e) {
    log_warn("Failed to write chromatin tracks: ", e$message)
  })
  clean_memory("chromatin tracks")

  # --- Gene coordinates ---
  tryCatch({
    chromatin_info$has_gene_coords <- write_gene_coordinates(
      seurat_obj, config, chromatin_assay
    )
  }, error = function(e) {
    log_warn("Failed to write gene coordinates: ", e$message)
  })
  clean_memory("gene coords")

  # --- Peak-gene links (only if pre-computed) ---
  tryCatch({
    chromatin_info$has_linked_peaks <- write_linked_peaks(
      seurat_obj, config, chromatin_assay
    )
  }, error = function(e) {
    log_warn("Failed to write linked peaks: ", e$message)
  })
  clean_memory("linked peaks")

  return(chromatin_info)
}

# =============================================================================
# CHROMATIN TRACKS (pseudobulk signal from peak matrix)
# =============================================================================

write_chromatin_tracks <- function(seurat_obj, config, chromatin_assay) {
  log_info("Generating pseudobulk chromatin tracks from peak matrix...")

  # Get cell type column
  cell_type_col <- get_cell_type_column(seurat_obj@meta.data)
  if (is.null(cell_type_col)) {
    log_warn("No cell type column found — skipping chromatin tracks")
    return(FALSE)
  }

  # Switch to chromatin assay
  DefaultAssay(seurat_obj) <- chromatin_assay

  # Get the peak matrix
  peak_matrix <- GetAssayData(seurat_obj, layer = "counts")
  peak_names <- rownames(peak_matrix)
  cell_types <- seurat_obj@meta.data[[cell_type_col]]

  log_info("Peak matrix: ", nrow(peak_matrix), " peaks x ", ncol(peak_matrix), " cells")
  log_info("Cell types: ", length(unique(cell_types)), " types from '", cell_type_col, "'")

  # Parse peak names (format: chr-start-end or chr:start-end)
  peak_coords <- parse_peak_names(peak_names)
  if (is.null(peak_coords)) {
    log_warn("Could not parse peak coordinates — skipping chromatin tracks")
    return(FALSE)
  }

  # Compute pseudobulk signal per cell type
  unique_types <- unique(cell_types)
  track_rows <- list()

  for (ct in unique_types) {
    ct_cells <- which(cell_types == ct)
    if (length(ct_cells) == 0) next

    # Sum peak counts across cells of this type, then normalize by cell count
    ct_signal <- Matrix::rowSums(peak_matrix[, ct_cells, drop = FALSE]) / length(ct_cells)

    # Only keep peaks with non-zero signal
    nonzero <- which(ct_signal > 0)
    if (length(nonzero) == 0) next

    track_rows[[length(track_rows) + 1]] <- data.frame(
      cell_type = rep(as.character(ct), length(nonzero)),
      chrom = peak_coords$chrom[nonzero],
      position_start = peak_coords$start[nonzero],
      position_end = peak_coords$end[nonzero],
      signal_value = as.numeric(ct_signal[nonzero]),
      stringsAsFactors = FALSE
    )
  }

  if (length(track_rows) == 0) {
    log_warn("No chromatin track data generated")
    return(FALSE)
  }

  tracks_df <- do.call(rbind, track_rows)
  write_parquet_safe(tracks_df, file.path(config$output_dir, "chromatin_tracks.parquet"))

  # Free peak matrix
  rm(peak_matrix, track_rows, tracks_df)
  clean_memory("tracks write")

  return(TRUE)
}

#' Parse peak names in chr-start-end format
parse_peak_names <- function(peak_names) {
  # Try chr-start-end format first
  parts <- strsplit(peak_names, "[-:]")

  # Validate first few entries
  valid <- sapply(head(parts, 10), function(p) length(p) >= 3)
  if (sum(valid) < 5) {
    log_warn("Peak names don't match expected chr-start-end format. Sample: ",
             paste(head(peak_names, 3), collapse = ", "))
    return(NULL)
  }

  data.frame(
    chrom = sapply(parts, `[`, 1),
    start = as.integer(sapply(parts, `[`, 2)),
    end = as.integer(sapply(parts, `[`, 3)),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# GENE COORDINATES
# =============================================================================

write_gene_coordinates <- function(seurat_obj, config, chromatin_assay) {
  log_info("Extracting gene coordinates...")

  # Try to get gene annotations from the ChromatinAssay
  assay_obj <- seurat_obj@assays[[chromatin_assay]]

  annotations <- tryCatch({
    Signac::Annotation(assay_obj)
  }, error = function(e) {
    log_warn("No annotations in ChromatinAssay: ", e$message)
    NULL
  })

  if (is.null(annotations) || length(annotations) == 0) {
    # Try EnsDb as fallback
    annotations <- tryCatch({
      library(EnsDb.Hsapiens.v86)
      genes <- ensembldb::genes(EnsDb.Hsapiens.v86)
      genes
    }, error = function(e) {
      log_warn("Could not load gene annotations: ", e$message)
      NULL
    })
  }

  if (is.null(annotations) || length(annotations) == 0) {
    log_warn("No gene annotations available — skipping gene_coordinates.parquet")
    return(FALSE)
  }

  # Filter to gene-level features
  if ("type" %in% names(GenomicRanges::mcols(annotations))) {
    gene_annot <- annotations[annotations$type == "gene"]
  } else {
    gene_annot <- annotations
  }

  if (length(gene_annot) == 0) {
    log_warn("No gene-level annotations found")
    return(FALSE)
  }

  # Extract gene names
  gene_name_col <- NULL
  mcol_names <- names(GenomicRanges::mcols(gene_annot))
  for (col in c("gene_name", "symbol", "SYMBOL", "gene_id")) {
    if (col %in% mcol_names) { gene_name_col <- col; break }
  }

  if (is.null(gene_name_col)) {
    log_warn("No gene name column found in annotations")
    return(FALSE)
  }

  # Build gene coordinates dataframe
  gene_names <- GenomicRanges::mcols(gene_annot)[[gene_name_col]]
  chroms <- as.character(GenomicRanges::seqnames(gene_annot))
  starts <- GenomicRanges::start(gene_annot)
  ends <- GenomicRanges::end(gene_annot)
  strands <- as.character(GenomicRanges::strand(gene_annot))

  # Try to get exon information
  exon_starts <- rep("", length(gene_annot))
  exon_ends <- rep("", length(gene_annot))

  tryCatch({
    # Get exons grouped by gene
    if ("type" %in% names(GenomicRanges::mcols(annotations))) {
      exon_annot <- annotations[annotations$type == "exon"]
      if (length(exon_annot) > 0 && gene_name_col %in% names(GenomicRanges::mcols(exon_annot))) {
        exon_by_gene <- split(exon_annot, GenomicRanges::mcols(exon_annot)[[gene_name_col]])

        for (i in seq_along(gene_names)) {
          gn <- gene_names[i]
          if (gn %in% names(exon_by_gene)) {
            exons <- exon_by_gene[[gn]]
            exon_starts[i] <- paste(GenomicRanges::start(exons), collapse = ",")
            exon_ends[i] <- paste(GenomicRanges::end(exons), collapse = ",")
          }
        }
      }
    }
  }, error = function(e) {
    log_warn("Could not extract exon info: ", e$message)
  })

  coords_df <- data.frame(
    gene_name = gene_names,
    chrom = chroms,
    start = starts,
    end = ends,
    strand = strands,
    exon_starts = exon_starts,
    exon_ends = exon_ends,
    stringsAsFactors = FALSE
  )

  # Remove duplicates (keep first occurrence)
  coords_df <- coords_df[!duplicated(coords_df$gene_name), ]

  write_parquet_safe(coords_df, file.path(config$output_dir, "gene_coordinates.parquet"))
  log_info("Wrote gene coordinates for ", nrow(coords_df), " genes")

  return(TRUE)
}

# =============================================================================
# PEAK-GENE LINKS (only if pre-computed in Seurat object)
# =============================================================================

write_linked_peaks <- function(seurat_obj, config, chromatin_assay) {
  log_info("Checking for pre-computed peak-gene links...")

  assay_obj <- seurat_obj@assays[[chromatin_assay]]

  # Check if links exist in the assay
  links <- tryCatch({
    Signac::Links(assay_obj)
  }, error = function(e) {
    NULL
  })

  if (is.null(links) || length(links) == 0) {
    log_info("No pre-computed peak-gene links found — skipping linked_peaks.parquet")
    return(FALSE)
  }

  log_info("Found ", length(links), " pre-computed peak-gene links")

  # Extract link data
  link_df <- data.frame(
    gene_name = links$gene,
    peak_chrom = as.character(GenomicRanges::seqnames(links)),
    peak_start = GenomicRanges::start(links),
    peak_end = GenomicRanges::end(links),
    correlation = if ("score" %in% names(GenomicRanges::mcols(links))) links$score else NA_real_,
    p_value = if ("pvalue" %in% names(GenomicRanges::mcols(links))) links$pvalue else NA_real_,
    link_type = "peak_gene",
    distance_to_tss = NA_integer_,
    stringsAsFactors = FALSE
  )

  write_parquet_safe(link_df, file.path(config$output_dir, "linked_peaks.parquet"))

  return(TRUE)
}
