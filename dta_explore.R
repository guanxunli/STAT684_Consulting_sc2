## load data
dta <- readRDS("data set/dta.rds")

## data explore
library(Seurat)
dim(dta)
dta_seurat <- CreateSeuratObject(counts = dta, project = "count", min.cells = 3, min.features = 200)
dta_seurat[["percent.mt"]] <- PercentageFeatureSet(dta_seurat, pattern = "^MT-")
dta_seurat <- subset(dta_seurat, subset = nFeature_RNA > 200  & percent.mt < 5)
dta_seurat <- NormalizeData(dta_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
## variable selection
dta_seurat <- FindVariableFeatures(dta_seurat, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(dta_seurat), 10)
all.genes <- rownames(dta_seurat)
dta_seurat <- ScaleData(dta_seurat, features = all.genes)
## dimension reduction
dta_seurat <- RunPCA(dta_seurat, features = VariableFeatures(object = dta_seurat), verbose = FALSE)
DimPlot(dta_seurat, reduction = "pca")
## Do cluster
dta_seurat <- FindNeighbors(dta_seurat, dims = 1:10)
dta_seurat <- FindClusters(dta_seurat, resolution = 0.5)
## dimension reduction
dta_seurat <- RunPCA(dta_seurat, features = VariableFeatures(object = dta_seurat), verbose = FALSE)
DimPlot(dta_seurat, reduction = "pca")
##Run non-linear dimensional reduction (UMAP/tSNE)
dta_seurat <- RunTSNE(dta_seurat, reduction = "pca", dims = 1:50)
DimPlot(dta_seurat, reduction = "tsne")
dta_seurat <- RunUMAP(dta_seurat, reduction = "pca", dims = 1:50)
DimPlot(dta_seurat, reduction = "umap")


