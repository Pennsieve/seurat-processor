# Pennsieve Seurat Processor
Pennsieve workflow component that processes Seurat Objects in a .rds file and exports the UMAP and metadata to a Parquet file.

## Description
This processor extracts viewer assets from Seurat objects that are stored using RDSWrite in R.
It stores the UMAP and TSNE reductions and the Seurate Metadata to a new Parquet file.

Metadata values are stored as categorical values in the Parquet file and a dictionary with unique values
is stored in the Parquet Metadata Key Value pairs. 

This file can be leveraged to provide a responsive UMAP and TSNE viewer in a web-application or for
downstream analysis without needing to download the original Seurat file. 

## Requirements
### Input
1. Input file with extension .Rds with a Seurat Object

### Output
1. A single Parquet file with the reductions and metadata 

## Notes
- Depending on the size of the .Rds file, the container might require significant memory as R only provides a way to load the entire Rds dataset in memory.

