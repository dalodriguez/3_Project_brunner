#!/bin/bash

# QUESTION 1
# How many features were detected in the four samples?
SRR_FILES=(SRR28910565 SRR28910566 SRR28910567 SRR28910568)
for srr in "${SRR_FILES[@]}"; do
    zcat "${srr}/outs/filtered_feature_bc_matrix/features.tsv.gz" | wc -l
done

# QUESTION 2
# Which sample has the highest number of detected cells in the raw feature count matrix?
for srr in "${SRR_FILES[@]}"; do
    echo $srr $(zcat "${srr}/outs/raw_feature_bc_matrix/barcodes.tsv.gz" | wc -l)
done

# QUESTION 3
# How many cells are conserved after Seurat quality control and preprocessing? 
Rscript -e '
seurat <- readRDS("C:/Users/dalod/Desktop/Projects/7_SepalAI/3_Project_brunner/Tanycytes/data/seurat_object.rds")
print(dim(seurat))
'

# QUESTION 4
# How many clusters were identified in the Seurat object after clustering using default resolution?
Rscript -e '
seurat <- readRDS("C:/Users/dalod/Desktop/Projects/7_SepalAI/3_Project_brunner/Tanycytes/data/seurat_object.rds")
print(length(unique(seurat$seurat_clusters)))
'

# QUESTION 5
# Which feature has the highest expression level in the final seurat object?
Rscript -e '
seurat <- readRDS("C:/Users/dalod/Desktop/Projects/7_SepalAI/3_Project_brunner/Tanycytes/data/seurat_object.rds")
print(which.max(rowMeans(seurat@assays$RNA@counts)))
'
