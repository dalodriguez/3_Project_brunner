setwd("C:/Users/dalod/Desktop/Projects/7_SepalAI/3_Project_brunner/Tanycytes/data/Combined/outs/count/filtered_feature_bc_matrix")

library(Seurat)
library(DropletUtils)
library(scuttle)
getwd()

sce <- read10xCounts(samples=getwd(), col.names=TRUE)
set.seed(100)
e.out <- emptyDrops(counts(sce))
sce <- sce[, which(e.out$FDR <= 0.01)]

sce <- addPerCellQC(sce, subsets = list(Mito = grep("^mt-", rowData(sce)$Symbol, value = FALSE)))

libsize.drop <- isOutlier(sce$sum, nmads = 3, type = "lower", log = TRUE)
feature.drop <- isOutlier(sce$detected, nmads = 3, type = "lower", log = TRUE)
mito.drop <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher")

sce <- sce[, !(libsize.drop | feature.drop | mito.drop)]
colnames(sce) <- sce$Barcode
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

seurat <- as.Seurat(sce)
seurat <- RenameAssays(seurat, originalexp = "RNA")
seurat <- SCTransform(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
seurat <- FindNeighbors(seurat)
seurat <- FindClusters(seurat)
seurat <- RunUMAP(seurat, dims = 1:30)
saveRDS(seurat, file = "C:/Users/dalod/Desktop/Projects/7_SepalAI/3_Project_brunner/Tanycytes/data/seurat_object.rds")
DimPlot(seurat)
