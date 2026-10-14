# ============================================================================
# chromatin.R — Multiome-specific: pseudobulk tracks, gene coords, peak-gene links
#
# All operations are wrapped in tryCatch — best-effort extraction.
# If Signac is not available or data is missing, these gracefully skip.
# ============================================================================

#' Write chromatin-related parquet files for multiome data
#'
#' Peak-gene links are written earlier (write_linked_peaks, before the RNA assay
#' is dropped); pass whether that succeeded as `has_linked_peaks`.
#' Returns a list of chromatin info for the manifest.
write_chromatin_data <- function(seurat_obj, config, chromatin_assay,
                                 has_linked_peaks = FALSE) {
  chromatin_info <- list(
    has_tracks = FALSE,
    has_gene_coords = FALSE,
    has_linked_peaks = isTRUE(has_linked_peaks)
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
# GENE ANNOTATIONS
# =============================================================================

#' Gene annotations for the chromatin assay, in Signac's format
#' (one row per exon/cds/utr/gap, with gene_name and type columns).
#' Falls back to Ensembl v86 (hg38) when the object has none. Cached per run.
get_gene_annotations <- function(seurat_obj, chromatin_assay) {
  if (!is.null(.annotation_cache$annotations)) return(.annotation_cache$annotations)

  annotations <- tryCatch(
    Signac::Annotation(seurat_obj[[chromatin_assay]]),
    error = function(e) NULL
  )

  if (is.null(annotations) || length(annotations) == 0) {
    log_info("No annotations in ChromatinAssay — using EnsDb.Hsapiens.v86 (hg38)")
    annotations <- tryCatch({
      suppressMessages(library(EnsDb.Hsapiens.v86))
      annot <- suppressWarnings(Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86))
      # Ensembl names chromosomes "1"; peaks use "chr1"
      GenomeInfoDb::seqlevelsStyle(annot) <- "UCSC"
      annot
    }, error = function(e) {
      log_warn("Could not load gene annotations: ", e$message)
      NULL
    })
  }

  if (!is.null(annotations) && length(annotations) > 0) {
    log_info("Gene annotations: ", length(annotations), " ranges, types: ",
             paste(unique(as.character(annotations$type)), collapse = ", "))
  }
  .annotation_cache$annotations <- annotations
  annotations
}
.annotation_cache <- new.env()

# =============================================================================
# GENE COORDINATES
# =============================================================================

write_gene_coordinates <- function(seurat_obj, config, chromatin_assay) {
  log_info("Extracting gene coordinates...")

  annotations <- get_gene_annotations(seurat_obj, chromatin_assay)
  if (is.null(annotations) || length(annotations) == 0) {
    log_warn("No gene annotations available — skipping gene_coordinates.parquet")
    return(FALSE)
  }

  coords_df <- build_gene_coordinates(annotations)
  if (is.null(coords_df) || nrow(coords_df) == 0) {
    log_warn("No gene coordinates could be built from annotations")
    return(FALSE)
  }

  write_parquet_safe(coords_df, file.path(config$output_dir, "gene_coordinates.parquet"))
  log_info("Wrote gene coordinates for ", nrow(coords_df), " genes")

  return(TRUE)
}

#' One row per gene: span over all its annotation ranges, plus exons merged
#' across transcripts (comma-separated starts/ends).
#'
#' Signac annotations have no "gene" rows — only exon/cds/utr/gap per
#' transcript — so the gene span is built from those.
build_gene_coordinates <- function(annotations) {
  mc <- GenomicRanges::mcols(annotations)
  name_col <- intersect(c("gene_name", "symbol", "SYMBOL", "gene_id"), names(mc))[1]
  if (is.na(name_col)) {
    log_warn("No gene name column found in annotations")
    return(NULL)
  }

  dt <- data.table::data.table(
    gene_name = as.character(mc[[name_col]]),
    chrom = as.character(GenomicRanges::seqnames(annotations)),
    start = GenomicRanges::start(annotations),
    end = GenomicRanges::end(annotations),
    strand = as.character(GenomicRanges::strand(annotations)),
    type = if ("type" %in% names(mc)) as.character(mc$type) else NA_character_
  )
  dt <- dt[!is.na(gene_name) & gene_name != ""]

  # Keep one chromosome per gene (e.g. PAR genes are annotated on X and Y)
  first_chrom <- dt[, .(chrom = chrom[1]), by = gene_name]
  dt <- dt[first_chrom, on = c("gene_name", "chrom")]

  genes <- dt[, .(start = min(start), end = max(end), strand = strand[1]),
              by = .(gene_name, chrom)]

  # Merge overlapping exons from all transcripts into one gene-level model
  exons <- dt[type == "exon", .(gene_name, start, end)]
  if (nrow(exons) > 0) {
    data.table::setorder(exons, gene_name, start)
    exons[, block := cumsum(start > data.table::shift(cummax(end), fill = 0L)),
          by = gene_name]
    merged <- exons[, .(start = min(start), end = max(end)), by = .(gene_name, block)]
    exon_strings <- merged[, .(
      exon_starts = paste(start, collapse = ","),
      exon_ends = paste(end, collapse = ",")
    ), by = gene_name]
    genes <- exon_strings[genes, on = "gene_name"]
  } else {
    genes[, `:=`(exon_starts = NA_character_, exon_ends = NA_character_)]
  }
  genes[is.na(exon_starts), `:=`(exon_starts = "", exon_ends = "")]

  as.data.frame(genes[, .(gene_name, chrom, start, end, strand, exon_starts, exon_ends)])
}

# =============================================================================
# PEAK-GENE LINKS
# =============================================================================

#' Write linked_peaks.parquet. Uses links already in the object if present,
#' otherwise computes them with Signac::LinkPeaks (correlating each gene's
#' expression with accessibility of peaks within 500 kb).
#'
#' Needs the RNA assay, so it runs before RNA is dropped. Never fails the run:
#' any error is logged and the file is skipped.
write_linked_peaks <- function(seurat_obj, config, chromatin_assay) {
  tryCatch({
    links <- Signac::Links(seurat_obj[[chromatin_assay]])
    if (length(links) > 0) {
      log_info("Found ", length(links), " pre-computed peak-gene links")
    } else if (!config$link_peaks) {
      log_info("No pre-computed peak-gene links and LINK_PEAKS=false — skipping")
      return(FALSE)
    } else {
      links <- compute_linked_peaks(seurat_obj, config, chromatin_assay)
    }

    if (is.null(links) || length(links) == 0) {
      log_info("No peak-gene links found — skipping linked_peaks.parquet")
      return(FALSE)
    }

    write_parquet_safe(links_to_df(links),
                       file.path(config$output_dir, "linked_peaks.parquet"))
    TRUE
  }, error = function(e) {
    log_warn("Failed to write peak-gene links: ", e$message)
    FALSE
  })
}

#' Run LinkPeaks in batches of genes, logging progress and stopping when the
#' time budget runs out (links found so far are kept).
compute_linked_peaks <- function(seurat_obj, config, chromatin_assay) {
  annotations <- get_gene_annotations(seurat_obj, chromatin_assay)
  if (is.null(annotations) || length(annotations) == 0) {
    log_warn("No gene annotations — cannot compute peak-gene links")
    return(NULL)
  }
  # LinkPeaks reads gene positions from the assay's annotation
  if (length(Signac::Annotation(seurat_obj[[chromatin_assay]])) == 0) {
    Signac::Annotation(seurat_obj[[chromatin_assay]]) <- annotations
  }
  annotated_genes <- unique(annotations$gene_name)

  seurat_obj <- add_peak_gc_content(seurat_obj, chromatin_assay)
  if (is.null(seurat_obj)) return(NULL)

  # LinkPeaks checks the expression layer against the default assay
  DefaultAssay(seurat_obj) <- "RNA"

  # Annotated genes, most widely expressed first, capped at LINK_MAX_GENES
  expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  n_expressing <- tabulate(expr@i + 1L, nbins = nrow(expr))
  names(n_expressing) <- rownames(expr)
  rm(expr)
  candidates <- n_expressing[names(n_expressing) %in% annotated_genes &
                               n_expressing > 10]
  genes <- names(sort(candidates, decreasing = TRUE))
  if (config$link_max_genes > 0) genes <- head(genes, config$link_max_genes)
  if (length(genes) == 0) {
    log_warn("No expressed, annotated genes to link")
    return(NULL)
  }

  batches <- split(genes, ceiling(seq_along(genes) / config$link_batch_size))
  log_info("Linking peaks for ", length(genes), " genes in ", length(batches),
           " batches (time budget: ",
           if (config$link_max_minutes > 0) paste(config$link_max_minutes, "min") else "none", ")")

  started <- proc.time()[["elapsed"]]
  results <- list()
  for (b in seq_along(batches)) {
    elapsed_min <- (proc.time()[["elapsed"]] - started) / 60
    if (config$link_max_minutes > 0 && elapsed_min > config$link_max_minutes) {
      log_warn("Time budget reached after ", b - 1, " / ", length(batches),
               " batches — keeping links found so far")
      break
    }

    linked <- NULL
    batch_result <- tryCatch({
      linked <- suppressMessages(Signac::LinkPeaks(
        object = seurat_obj,
        peak.assay = chromatin_assay,
        expression.assay = "RNA",
        genes.use = batches[[b]],
        verbose = FALSE
      ))
      Signac::Links(linked[[chromatin_assay]])
    }, error = function(e) {
      log_warn("LinkPeaks batch ", b, " failed: ", e$message)
      NULL
    })
    rm(linked)

    if (length(batch_result) > 0) results[[length(results) + 1]] <- batch_result
    log_info("  Linked batch ", b, " / ", length(batches), " — ",
             length(batch_result), " links, ",
             round((proc.time()[["elapsed"]] - started) / 60, 1), " min elapsed")
    clean_memory(paste0("link batch ", b))
  }

  if (length(results) == 0) return(NULL)
  do.call(c, unname(results))
}

#' LinkPeaks compares each peak against background peaks with similar GC
#' content, which needs the genome sequence. Skips links if unavailable.
add_peak_gc_content <- function(seurat_obj, chromatin_assay) {
  meta <- seurat_obj[[chromatin_assay]][[]]
  if (all(c("GC.percent", "sequence.length") %in% colnames(meta))) {
    return(seurat_obj)
  }

  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
    log_warn("BSgenome.Hsapiens.UCSC.hg38 not installed — cannot compute peak-gene links")
    return(NULL)
  }

  log_info("Computing peak GC content (RegionStats, hg38)...")
  tryCatch(
    Signac::RegionStats(seurat_obj, assay = chromatin_assay,
                        genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38),
    error = function(e) {
      log_warn("RegionStats failed: ", e$message)
      NULL
    }
  )
}

#' Flatten Signac links to the linked_peaks.parquet schema.
#' Signac stores each link as a range from the gene TSS to the peak midpoint;
#' the peak itself is in the `peak` column ("chr-start-end").
links_to_df <- function(links) {
  peaks <- Signac::StringToGRanges(links$peak, sep = c("-", "-"))
  peak_mid <- GenomicRanges::start(peaks) + GenomicRanges::width(peaks) / 2
  link_start <- GenomicRanges::start(links)
  link_end <- GenomicRanges::end(links)
  tss <- ifelse(abs(link_start - peak_mid) < abs(link_end - peak_mid), link_end, link_start)
  mc <- GenomicRanges::mcols(links)

  data.frame(
    gene_name = as.character(links$gene),
    peak_chrom = as.character(GenomicRanges::seqnames(peaks)),
    peak_start = GenomicRanges::start(peaks),
    peak_end = GenomicRanges::end(peaks),
    correlation = if ("score" %in% names(mc)) as.numeric(mc$score) else NA_real_,
    zscore = if ("zscore" %in% names(mc)) as.numeric(mc$zscore) else NA_real_,
    p_value = if ("pvalue" %in% names(mc)) as.numeric(mc$pvalue) else NA_real_,
    link_type = "peak_gene",
    distance_to_tss = as.integer(round(peak_mid - tss)),
    stringsAsFactors = FALSE
  )
}
