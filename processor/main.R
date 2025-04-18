library(nanoparquet)
library(jsonlite)
library(SeuratObject)
library(yaml)

inputFolder = Sys.getenv('INPUT_DIR', "/data/input")
print(inputFolder)
inputFile = list.files(path = inputFolder)
print(inputFile)
data <- readRDS(file.path(inputFolder,inputFile[1]))

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

outputFolder = Sys.getenv('OUTPUT_DIR',"/data/output")

# Viewer Config File
viewerInfo <- list(name = "parquet-umap-viewer", options = list())
write_yaml(viewerInfo, file.path(outputFolder, "viewer_config.yml"), fileEncoding = "UTF-8")

# Write Parquet File
write_parquet(result, file.path(outputFolder, "results.parquet"),
 meta=metadata, compression= "snappy")