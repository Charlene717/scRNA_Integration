# 1. Ensure the Harmony package is installed and loaded
if(!require('harmony')) { install.packages('harmony'); library(harmony) }

# 2. Perform PCA before batch correction
seuratObject_Sample_RunHar <- RunPCA(seuratObject_Sample_RunHar, npcs = 50, verbose = FALSE)
# Harmony adjusts data in the PCA space, so PCA must be computed first.
# `npcs = 50` specifies the number of principal components to calculate.
# `verbose = FALSE` suppresses excessive output messages.

# 3. Apply batch correction using Harmony
seuratObject_Sample_RunHar <- RunHarmony(seuratObject_Sample_RunHar, group.by.vars = "FullName")
# `group.by.vars = "FullName"` specifies the metadata variable for batch correction.
# Harmony adjusts the PCA space to make data from different batches more comparable.

# 4. Perform UMAP for visualization after batch correction
seuratObject_Sample_RunHar <- RunUMAP(seuratObject_Sample_RunHar, reduction = "harmony", dims = 1:50)
# `reduction = "harmony"` ensures UMAP uses Harmony-corrected PCA.
# `dims = 1:50` selects the first 50 principal components for UMAP calculation.

# 5. Compute nearest neighbors to prepare for clustering
seuratObject_Sample_RunHar <- FindNeighbors(seuratObject_Sample_RunHar, reduction = "harmony", dims = 1:50)
# `FindNeighbors()` computes cell-to-cell similarity, a prerequisite for clustering.

# 6. Perform clustering using the Seurat Louvain/Leiden algorithm
seuratObject_Sample_RunHar <- FindClusters(seuratObject_Sample_RunHar, resolution = 0.5)
# `resolution = 0.5` controls clustering granularity: higher values yield more clusters.

# 7. Visualize UMAP to examine clustering and batch effects
DimPlot(seuratObject_Sample_RunHar, group.by = "orig.ident", reduction = "umap")  # Original sample
DimPlot(seuratObject_Sample_RunHar, group.by = "FullName", reduction = "umap")   # Batch information
DimPlot(seuratObject_Sample_RunHar, group.by = "Tissue", reduction = "umap")     # Tissue type
DimPlot(seuratObject_Sample_RunHar, group.by = "seurat_clusters", reduction = "umap") # Seurat clustering result

# 8. Compare tissue-based clustering with Seurat's clustering
DimPlot(seuratObject_Sample_RunHar, group.by = "Tissue", reduction = "umap") +
  DimPlot(seuratObject_Sample_RunHar, group.by = "seurat_clusters", reduction = "umap")


## Comparison: Seurat vs Harmony Integration
library(patchwork) # Ensure the patchwork package is installed for combining plots

p1 <- DimPlot(seuratObject_Sample, group.by = "FullName", reduction = "umap") + 
  ggtitle("Seurat Integration")  # Title for Seurat-based integration result

p2 <- DimPlot(seuratObject_Sample_RunHar, group.by = "FullName", reduction = "umap") + 
  ggtitle("Harmony Batch Correction (Aligned PCA)") # Title for Harmony-based integration result

p1 + p2  # Combine the two UMAP plots for comparison


## Export RData
save.image("C:/Charlene/Code_GitHub_BioInport2025/KGD_CellTypeAnnot_Skin/20250304_EPPK_GSE202352_Harmony.RData")

## Export rds
saveRDS(seuratObject_Sample, file = "20250304_EPPK_GSE202352_Harmony.rds")

