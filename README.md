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
| `gene_coordinates.parquet` | Gene genomic coordinates (multiome only) |
| `linked_peaks.parquet` | Peak-gene linkages if pre-computed (multiome only) |
| `manifest.json` | Dataset manifest describing all output files |
| `viewer_config.yml` | Pennsieve platform viewer configuration |

### Consolidated Mode (`OUTPUT_MODE=consolidated`)

Original single-file output: one `results.parquet` + `viewer_config.yml`.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `OUTPUT_MODE` | `expanded` | `expanded` (multi-file) or `consolidated` (legacy single file) |
| `CHUNK_SIZE` | `500` | Number of genes per chunk file in expanded mode |

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
- After Phase 3 (Gene expression): RNA assay is removed from memory
- `gc()` is called between every phase with memory usage logging

### Why This Is Necessary

R's `readRDS()` loads the entire object into memory at once — there is no way to stream or partially read an RDS file. A 31 GB multiome RDS file consumes ~57 GB of RAM after deserialization, leaving very little headroom on a 62 GB machine. Without assay dropping, gene expression processing would OOM.

## Multiome Detection

Multiome objects (RNA + ATAC) are auto-detected by checking for `ChromatinAssay` class in the assay list. When detected:
- Multiple embeddings are extracted (WNN, ATAC, RNA UMAP)
- Chromatin tracks are generated from the peak matrix
- Gene coordinates and peak-gene linkages are extracted if available

## Requirements

### Input
1. An .rds file containing a Seurat object

### Output
Depends on `OUTPUT_MODE` — see Output Modes above.

## Notes
- The container requires 64 GB memory for large multiome objects
- Processing a 31 GB RDS file takes approximately 20-30 minutes
- The processor runs on the Pennsieve workflow platform as a Docker container
