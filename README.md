# Pennsieve Seurat Processor

Pennsieve workflow component that processes Seurat Objects (.rds) and exports UMAP embeddings, cell metadata, and gene expression data to Parquet files for web visualization via PrecisionDashWidgets.

## Output Modes

### Expanded Mode (default, `OUTPUT_MODE=expanded`)

Produces multiple Parquet files optimized for lazy-loading in the browser:

| File | Description |
|------|-------------|
| `umap_complete.parquet` | Default embedding (WNN UMAP for multiome, UMAP for scRNA) |
| `wnn_umap_complete.parquet` | WNN UMAP coordinates (multiome only) |
| `atac_umap_complete.parquet` | ATAC UMAP coordinates (multiome only) |
| `rna_umap_complete.parquet` | RNA UMAP coordinates (multiome only, computed from PCA if missing) |
| `cells.parquet` | Cell metadata (all columns, factors converted to strings) |
| `genes.parquet` | Gene name to gene_id (0-indexed) mapping |
| `gene_stats.parquet` | Per-gene statistics (mean, total, pct_cells, n_expressing) |
| `gene_locations.parquet` | Maps each gene to its chunk file |
| `chunks/chunk_NNNNN.parquet` | Gene expression data grouped by chunk |
| `chromatin_tracks.parquet` | Pseudobulk chromatin signal per cell type (multiome only) |
| `gene_coordinates.parquet` | Gene span, strand, and merged exons per gene (multiome only) |
| `linked_peaks.parquet` | Peak-gene links: pre-computed in the object, or computed with Signac `LinkPeaks` (multiome only). Columns: gene, peak location, correlation, z-score, p-value, `link_type` (`promoter` within 2 kb of the TSS, `intronic` inside the gene, else `distal`), `distance_to_tss` (negative = upstream), `top_cell_type` and its mean signal |
| `manifest.json` | Dataset manifest describing all output files |
| `viewer_config.yml` | Pennsieve platform viewer configuration |

### Consolidated Mode (`OUTPUT_MODE=consolidated`)

Original single-file output: one `results.parquet` + `viewer_config.yml`.

## Data Provenance: What Is Copied vs. Calculated

Most outputs are copied straight from the Seurat object. The values below are
**calculated by this processor** (or chosen by its rules) and are not part of
the submitted data. They depend on our choices and can differ from what the
researchers would get with their own analysis. Dashboards should label them.

### Copied from the object (no calculation)

- WNN and ATAC UMAP coordinates (`wnn.umap`, `atac.umap`)
- Cell metadata (`cells.parquet`), including cell type labels
- Gene expression values: the RNA `data` layer (log-normalized), unchanged
- Gene locations, strands and exons: from the gene annotation stored in the chromatin assay.
  The processor falls back to Ensembl v86 (hg38) **only** if the object has none, and logs it
- Peak locations (the chromatin assay's peak set)

### Calculated by the processor

| Output | What it is | How it is made | Room for error |
|--------|-----------|----------------|----------------|
| RNA UMAP (`rna_umap_complete.parquet`) | Cell map from RNA only | `RunUMAP` on the object's PCA (first 30 dims, Seurat defaults), **only when the object has no RNA UMAP** | Not the researchers' own RNA UMAP; layout differs between runs and from theirs |
| Which cell type column is used | Groups for tracks and top cell type | First of `final_cluster_ids`, `cell_type`, `celltype`, ... found in the metadata | Wrong column if the object names its labels differently |
| Which chromatin assay is used | Peak set for tracks and links | `ATAC_macs2` if present, else `ATAC` | Results reflect that peak set only |
| `gene_stats.parquet` | Mean, total, % of cells expressing | Summed from the log-normalized `data` layer | Simple sums; means are of log values |
| Gene span and exons (`gene_coordinates.parquet`) | One row per gene | Span = outermost transcript ends; exons merged across all transcripts | Shows a union of transcripts, not one isoform; genes on several chromosomes keep the first |
| Chromatin tracks (`chromatin_tracks.parquet`) | Accessibility per cell type | Mean peak counts per cell in each cell type (pseudobulk); zero-signal peaks omitted | Not normalized for sequencing depth (no TF-IDF / CPM); per peak, not base-pair coverage |
| **Peak-gene links (`linked_peaks.parquet`)** | Regions linked to a gene | **Computed with Signac `LinkPeaks` when the object has no links.** Defaults: peaks within 500 kb of the TSS, Pearson correlation of peak counts vs. RNA `data`, GC-matched background (hg38 genome assumed), kept if p < 0.05 and \|r\| > 0.05 | Statistical association, not proof of regulation. Only the genes chosen by `LINK_GENES` / most variable (`LINK_MAX_GENES`) are linked; other genes show no links even if links exist. Correlations in single-cell data are low (typically 0.05–0.4) |
| `correlation`, `zscore`, `p_value` | Link strength | From `LinkPeaks` (above) | p-values are per link, not corrected across all links |
| `link_type` | promoter / intronic / distal | **Our rule:** `promoter` if the peak is within 2 kb of the TSS; `intronic` if its midpoint is inside the gene span; else `distal`. Peaks overlapping exons count as intronic | Cut-offs are a convention, not from the data. Uses the gene span (all transcripts), so a peak in one isoform's intron may be another's exon |
| `distance_to_tss` | Distance from gene start to peak | Peak midpoint minus the TSS used by `LinkPeaks`, in bp, **in the gene's direction**: negative = upstream (before the gene starts), positive = downstream (into or past the gene) | TSS is the start of the longest transcript; other isoforms may start elsewhere |
| `top_cell_type`, `top_cell_type_signal` | Cell type where the peak is most open | Highest mean peak count per cell among cell types | Ties go to the first type alphabetically; small cell types can be noisy |

If the object already contains links (`Links()` on the chromatin assay), those
are used as-is and only `link_type`, `distance_to_tss` and `top_cell_type` are
calculated.

### Provenance in `manifest.json`

Because this depends on the input, each run records it in the manifest:

- Each entry in `embeddings` has `source`: `"object"` (copied) or `"computed"` (with `method`).
- `provenance` has one entry per output (`gene_coordinates`, `chromatin_tracks`,
  `linked_peaks`, `gene_stats`, `cell_type_column`, `chromatin_assay`) with
  `source` (`object`, `computed` or `selected`) and how it was made.
  `linked_peaks` also lists the `LinkPeaks` parameters, how genes were chosen,
  whether the run finished (`complete`), and `columns` with the source and rule
  for each column.

Dashboards should read these fields to mark calculated values, not assume them.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `OUTPUT_MODE` | `expanded` | `expanded` (multi-file) or `consolidated` (legacy single file) |
| `CHUNK_SIZE` | `500` | Number of genes per chunk file in expanded mode |
| `LINK_PEAKS` | `true` | Compute peak-gene links with `LinkPeaks` when the object has none |
| `LINK_GENES` | _(empty)_ | Genes to always link, comma-separated (e.g. `SCN10A,TAC1`); linked first |
| `LINK_MAX_GENES` | `5000` | Link the N most variable genes (`0` = all, about 50 min per 1000 genes on a 66k-cell dataset) |
| `LINK_BATCH_SIZE` | `1000` | Genes per `LinkPeaks` call (progress is logged per batch) |
| `LINK_MAX_MINUTES` | `0` | Time budget for linking (`0` = no limit); links found so far are kept when it runs out |

## Memory Management

The processor handles large Seurat objects (30+ GB) with aggressive memory optimization:

### Assay Dropping

After loading the RDS file, the processor **drops unused assays from memory** to free RAM for processing. This does NOT modify the source .rds file.

**Retained assays:**
- `RNA` — used for gene expression extraction (Phase 3)
- `chromvar` — motif activity scores, retained for future dashboard features
- Chromatin assay (`ATAC_macs2` or `ATAC`) — used for chromatin track extraction (Phase 4)

**Dropped assays:**
- `SCT` — SCTransform normalization (redundant with RNA `data` layer for visualization)
- `integrated` — batch-corrected values (used for clustering, not individual gene queries)
- `ATAC` — original peaks (superseded by ATAC_macs2 when present)

### Phase-Based Cleanup

The processor also strips data between processing phases:
- After Phase 1 (Embeddings): reductions are removed from memory
- After Phase 3 (Gene expression) and Phase 3b (Peak-gene links): RNA assay is removed from memory. Linking runs before the drop because it correlates RNA expression with ATAC accessibility
- `gc()` is called between every phase with memory usage logging

### Why This Is Necessary

R's `readRDS()` loads the entire object into memory at once — there is no way to stream or partially read an RDS file. A 31 GB multiome RDS file consumes ~57 GB of RAM after deserialization, leaving very little headroom on a 62 GB machine. Without assay dropping, gene expression processing would OOM.

## Multiome Detection

Multiome objects (RNA + ATAC) are auto-detected by checking for `ChromatinAssay` class in the assay list. When detected:
- Multiple embeddings are extracted (WNN, ATAC, RNA UMAP)
- Chromatin tracks are generated from the peak matrix
- Gene coordinates are built from the ChromatinAssay annotation (or Ensembl v86 / hg38 if it has none)
- Peak-gene links are exported if pre-computed, otherwise computed with `LinkPeaks` (needs `BSgenome.Hsapiens.UCSC.hg38` for peak GC content; skipped with a warning if unavailable)

## Requirements

### Input
1. An .rds file containing a Seurat object

### Output
Depends on `OUTPUT_MODE` — see Output Modes above.

## Notes
- The container requires 64 GB memory for large multiome objects
- Processing a 31 GB RDS file takes approximately 20-30 minutes
- The processor runs on the Pennsieve workflow platform as a Docker container
