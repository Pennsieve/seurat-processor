library(nanoparquet)
library(jsonlite)
library(SeuratObject)

data <- readRDS("/data/input/TG_neurons_release.Rds")

variables = c()
for (r in names(data@reductions)) {
  dimNames <- names(data@reductions[[r]])
  variables <- c(variables, dimNames)
}
variables <- c(variables, names(data@meta.data))


result <- FetchData(data, variables)

metadata <- c()
for (n in names(data@meta.data)) {
  result[[n]] <- as.factor(result[[n]])
  metadata <- c(metadata, setNames(c(toJSON(unique(result[[n]]))), c(n) ))
}

write_parquet(result, "/data/output/TG_neurons_release.parquet",
 meta=metadata, compression= "snappy")