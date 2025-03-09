library(nanoparquet)
library(jsonlite)
library(SeuratObject)

data <- readRDS("/data/input/DRG_nonneurons_release.rds")

variables = c()
for (r in names(data@reductions)) {
  if (r == 'umap' || r == 'tsne') {
    dimNames <- names(data@reductions[[r]])
    variables <- c(variables, dimNames)
  }
}
variables <- c(variables, names(data@meta.data))


result <- FetchData(data, variables)

metadata <- c()
for (n in names(data@meta.data)) {
  result[[n]] <- as.factor(result[[n]])
  metadata <- c(metadata, setNames(c(toJSON(unique(result[[n]]))), c(n) ))
}

write_parquet(result, "/data/output/DRG_nonneurons_release.parquet",
 meta=metadata, compression= "snappy")