##### Presetting ######
rm(list = ls()) # Clean all variables ##* Comment out if running the entire script

## Speed up (Old Version)
memory.limit(150000)

## Speed up
if(!require('future')) install.packages('future'); library(future)
## https://github.com/immunogenomics/presto
if(!require('presto')) devtools::install_github("immunogenomics/presto"); library(presto) # Speeds up FindAllMarkers
plan(multicore, workers = 20)
options(future.globals.maxSize = 2048*100 * 1024^2) # Set memory limit to ~204.8 GB

#### Set parameters ####
# source("Set_Parameter.R")
# 
# if (!require('readr')) {install.packages('readr'); library(readr)}
# writeLines(readLines("Set_Parameter.R"), 
#            paste0(Name_ExportFolder,"/", Name_Export,"_Set_Parameter_Content.txt"))


#### Load Packages ####
if(!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if(!require('tidyverse')) {install.packages('tidyverse'); library(tidyverse)}

if(!require('devtools')) {install.packages('devtools'); library(devtools)}
if(!require('scMRMA')) {devtools::install_github("JiaLiVUMC/scMRMA"); library(scMRMA)}

#### Function ####
library(Seurat)
Seurat_Prepocessing <- function(seurat_obj, Num_PCA = 50 ,Set_nfeatures = 2000 ) {
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- seurat_obj  %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = Set_nfeatures) %>%
    ScaleData() %>%
    RunPCA(npcs = Num_PCA) %>%
    FindNeighbors(dims = 1:Num_PCA) %>%
    FindClusters() %>% # resolution = 0.5
    # RunTSNE(dims = 1:Num_PCA) %>%
    RunUMAP(dims = 1:Num_PCA)
}
# seuratObject_Sample <- Seurat_Prepocessing(seuratObject_Sample)



#### Load Data ####
## EPPK
load("C:/Charlene/Code_GitHub/KGD_10x_PreProcessing/Export_2025011315RGZ_EPPK_RefselfDB_SetClt_seurat_clusters/CellTypeAnnot/2025011315RGZ_ROGUE.RData")
rm(list = setdiff(ls(), c("seuratObject_Sample", "Seurat_Prepocessing")))

DimPlot(seuratObject_Sample, reduction = "umap")
DimPlot(seuratObject_Sample, group.by = "seurat_clusters", reduction = "umap")
DimPlot(seuratObject_Sample, group.by = "orig.ident", reduction = "umap")


seuratObject_Sample_EPPK <- seuratObject_Sample[, seuratObject_Sample$orig.ident == "HS23378-2862_EPPK"]
rm(seuratObject_Sample)

seuratObject_Sample_EPPK$FullName <- "KGD_2862_EPPK_sole"
seuratObject_Sample_EPPK$Tissue <- "sole_EPPK"

seuratObject_Sample_EPPK <- Seurat_Prepocessing(seuratObject_Sample_EPPK)

DimPlot(seuratObject_Sample_EPPK, group.by = "orig.ident", reduction = "umap")

## Normal
seuratObject_Sample_Normal_List <- readRDS(file = "C:/Charlene/Code_GitHub_BioInport2025/KGD_CellTypeAnnot_Skin/Input/GSE202352_combined_seurat_list.rds")

#### Integration ####
# library(Seurat)  # 請確保已在函式外部載入 Seurat
seuratObject_Sample_EPPK <- JoinLayers(seuratObject_Sample_EPPK)
Layers(seuratObject_Sample_EPPK)



# 1. 載入 Seurat 套件
library(Seurat)

# 2. 建立一個 List，放入欲整合的物件
seurat.list <- c(
  EPPK = seuratObject_Sample_EPPK,
  seuratObject_Sample_Normal_List
)

# 3. 選擇合適的整合特徵（integration features）
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)

# 4. 找到整合的錨點（Integration Anchors）
anchors <- FindIntegrationAnchors(
  object.list = seurat.list,
  anchor.features = features,
  dims = 1:50
)


# 5. 進行批次效應校正並整合
seuratObject_Sample <- IntegrateData(anchorset = anchors, dims = 1:50)

# 6. 整合後的資料會存放在 "integrated" assay，後續分析記得切換至該 assay
DefaultAssay(seuratObject_Sample) <- "integrated"

# 7. 可以進一步進行 ScaleData、RunPCA、RunUMAP 等標準下游分析
seuratObject_Sample <- ScaleData(seuratObject_Sample)
seuratObject_Sample <- RunPCA(seuratObject_Sample)
seuratObject_Sample <- FindNeighbors(seuratObject_Sample, dims = 1:50) # dims可視需要調整
seuratObject_Sample <- FindClusters(seuratObject_Sample, resolution = 0.5)
seuratObject_Sample <- RunUMAP(seuratObject_Sample, dims = 1:50, min.dist = 0.1)


DimPlot(seuratObject_Sample, group.by = "orig.ident", reduction = "umap")
DimPlot(seuratObject_Sample, group.by = "FullName", reduction = "umap")
DimPlot(seuratObject_Sample, group.by = "Tissue", reduction = "umap") +
DimPlot(seuratObject_Sample, group.by = "seurat_clusters", reduction = "umap")

# save.image("C:/Charlene/Code_GitHub_BioInport2025/KGD_CellTypeAnnot_Skin/20250304_EPPK_GSE202352.RData")
# saveRDS(seuratObject_Sample, file = "20250304_EPPK_GSE202352.rds")

rm(seurat.list,anchors, seuratObject_Sample_EPPK,seuratObject_Sample_Normal_List)

seuratObject_Sample <- JoinLayers(seuratObject_Sample, assay = "RNA")
DefaultAssay(seuratObject_Sample) <- "RNA"
# DefaultAssay(seuratObject_Sample) <- "integrated"
Layers(seuratObject_Sample)



#### ROGUE ####
if(!require("ROGUE")) {devtools::install_github("PaulingLiu/ROGUE"); library(ROGUE)}

# Extract the expression matrix
expr <- seuratObject_Sample[["RNA"]]$data %>% as.data.frame()
meta <- seuratObject_Sample@meta.data

## Filter low-abundance genes and low-quality cells
expr <- matr.filter(expr, min.cells = 10, min.genes = 10)

## Calculate expression entropy
ent.res <- SE_fun(expr)

## Perform ROGUE calculation based on the specified column
rogue.value <- CalculateRogue(ent.res, platform = "UMI")

## Calculate ROGUE values for each cluster
cell_type_col <- "seurat_clusters"
rogue.res <- rogue(expr,
                   labels  = meta[[cell_type_col]],  # Use the specified column as labels
                   samples = meta[[cell_type_col]],  # Optionally, replace with another column
                   platform = "UMI",
                   span = 0.6)

DimPlot(seuratObject_Sample, reduction = "umap", label = TRUE, group.by = "seurat_clusters") + rogue.boxplot(rogue.res)


#### Cell Type Annotation by scMRMA ####
Set_scMRMA_P <- 0.05
Set_SubClust_k <- 20

# source("DBMarker_Skin_ChatGPT_Pro_Keratinocyte.R")
# source("DBMarker_Skin_Expert.R")
# source("DBMarker_Skin_Combine.R")

source("DBMarker_Skin_Combine_Lab_M2_Fib.R")

DB_Skin_Marker_Test <- df_all_markers_ChatGPT

DefaultAssay(seuratObject_Sample) <- "RNA"
# Call scMRMA
result_scMRMA <- scMRMA(
  input = seuratObject_Sample,
  p = Set_scMRMA_P, # 0.05,
  normalizedData = F,
  selfDB = DB_Skin_Marker_Test,
  selfClusters = seuratObject_Sample$seurat_clusters,
  # selfClusters = NULL,
  k = Set_SubClust_k # 20
)

#### Visualization ####
source("FUNPlot_Bar.R")

# Add scMRMA results to meta.data
seuratObject_Sample@meta.data$label_scMRMA_Level1 <- result_scMRMA[["multiR"]][["annotationResult"]]$Level1
seuratObject_Sample@meta.data$label_scMRMA_Level2 <- result_scMRMA[["multiR"]][["annotationResult"]]$Level2
seuratObject_Sample@meta.data$label_scMRMA_Level3 <- result_scMRMA[["multiR"]][["annotationResult"]]$Level3
seuratObject_Sample@meta.data$label_scMRMA_Level4 <- result_scMRMA[["multiR"]][["annotationResult"]]$Level4
DimPlot(seuratObject_Sample, group.by = "label_scMRMA_Level4", reduction = "umap", label = T)

# Process Level 1 to Level 4 using a for loop
plots <- list() # Save all Level plots in a list
for (level in 1:4) {
  label_column <- paste0("label_scMRMA_Level", level)
  
  # Combine plots
  combined_plot <- patchwork::wrap_plots(
    DimPlot(seuratObject_Sample, reduction = "umap", label = TRUE, group.by = label_column) +
      ggtitle(paste("UMAP - Cell Type - Level", level)),
    DimPlot(seuratObject_Sample, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
      ggtitle(paste("UMAP - Cluster - Level", level)),
    plot_seurat_meta_data(seurat_object = seuratObject_Sample, 
                          x_var = "seurat_clusters", 
                          fill_var = label_column)$proportion_plot,
    plot_seurat_meta_data(seurat_object = seuratObject_Sample, 
                          x_var = "seurat_clusters", 
                          fill_var = label_column)$count_plot,
    ncol = 2
  )
  
  # Display the plot immediately
  print(combined_plot)
  
  # Save the plot to the list
  plots[[level]] <- combined_plot
}



################################################################################
#### Visualization ####
seuratObject_Sample$label_scMRMA_Level4[
  seuratObject_Sample$label_scMRMA_Level4 == "Unassigned"
] <- "Melanocytes"

# seuratObject_Sample@meta.data$Cell_Type <- result_scMRMA[["multiR"]][["annotationResult"]]$Level4
seuratObject_Sample@meta.data$Cell_Type <- seuratObject_Sample$label_scMRMA_Level4

DimPlot(seuratObject_Sample, group.by = "seurat_clusters", reduction = "umap", label = T) +
  DimPlot(seuratObject_Sample, group.by = "Cell_Type", reduction = "umap", label = T)


Set_Com_broad_cell_clusters <- TRUE
if(Set_Com_broad_cell_clusters){
  source("FUN_Annotate_broad_cell_clusters.R")
  seuratObject_Sample <- annotate_broad_cell_clusters(
    seurat_object = seuratObject_Sample,
    broad_cell_type_col = "Cell_Type",
    seurat_cluster_col = "seurat_clusters"
  )
} 


DimPlot(seuratObject_Sample, group.by = "orig.ident", reduction = "umap", label = T) +
  DimPlot(seuratObject_Sample, group.by = "seurat_clusters", reduction = "umap", label = T) +
  DimPlot(seuratObject_Sample, group.by = "label_scMRMA_Level4", reduction = "umap", label = T) +
  DimPlot(seuratObject_Sample, group.by = "BroadCellTypeAnnot_SeuratClusters", reduction = "umap", label = T)

DimPlot(seuratObject_Sample, group.by = "label_scMRMA_Level4", reduction = "umap", label = T)
DimPlot(seuratObject_Sample, group.by = "BroadCellTypeAnnot_SeuratClusters", reduction = "umap", label = T)


######################################################################################
#### Heatmap ####

# 安裝並載入必要的 R 套件
if (!require('Seurat')) { install.packages('Seurat'); library(Seurat) }
if (!require('pheatmap')) { install.packages('pheatmap'); library(pheatmap) }
if (!require('dplyr')) { install.packages('dplyr'); library(dplyr) }

# 確保 `seuratObject_Sample` 具有細胞類型標註
Idents(seuratObject_Sample) <- seuratObject_Sample$label_scMRMA_Level4

# **1. 找出各 Cell Type 的 Marker Genes**
marker_genes <- FindAllMarkers(
  object = seuratObject_Sample,
  only.pos = TRUE,   # 只取表達上調的基因
  min.pct = 0.25,    # 至少在 25% 的細胞中表達
  logfc.threshold = 0.25  # 最低 log fold change
)

# **2. 選取每個 cell type 的前 10 個 marker genes**
top_markers <- marker_genes %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# **3. 擷取這些基因的表達矩陣**
expr_matrix <- FetchData(
  seuratObject_Sample, 
  vars = unique(top_markers$gene)
)

# **4. 轉換矩陣格式 (確保 row 為 gene，col 為 cell)**
expr_matrix <- t(expr_matrix)

# **5. 繪製 Heatmap**
pheatmap(
  expr_matrix,
  scale = "row",  # 標準化每列
  clustering_method = "ward.D2", 
  show_rownames = TRUE, 
  show_colnames = FALSE, 
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Cell Type Marker Gene Heatmap"
)



DoHeatmap(
  object = seuratObject_Sample, 
  features = top_markers$gene, 
  group.by = "label_scMRMA_Level4"
) + 
  scale_fill_gradient2(
    low = "#004B97", 
    mid = "white", 
    high = "#D9006C", 
    midpoint = 0
  ) + 
  theme(
    axis.text.y = element_text(size = 8)
  )




######################################################################################

DefaultAssay(seuratObject_Sample) <- "RNA"

library(Seurat)
library(ggplot2)

# 設定目標基因
target_genes <- c("EPGN", "KRT1",  "KRT5", 
                  "EGFR","KRT9","KRT14")

target_genes <- c("EPGN", "EGFR", "KRT1", 
                  "KRT9", "KRT5", "KRT14")

FeaturePlot(seuratObject_Sample, features = target_genes, ncol = 3)



# 小提琴圖 (Violin Plot) - 分開顯示 Normal vs EPPK
p1 <- VlnPlot(seuratObject_Sample, features = target_genes, 
              group.by = "BroadCellTypeAnnot_SeuratClusters", 
              split.by = "Type",  # 這行會讓 Normal 和 EPPK 分開顯示
              pt.size = 0, # 不顯示點
              combine = TRUE) 
# 顯示圖表
p1

library(ggplot2)

p1 <- VlnPlot(seuratObject_Sample, features = target_genes, 
              group.by = "BroadCellTypeAnnot_SeuratClusters", 
              split.by = "Type",  
              pt.size = 0, 
              combine = TRUE) + 
  theme(legend.position = "right")  # 顯示圖例在右側

# 顯示圖表
p1


# 氣泡圖 (Bubble Plot) - 分開顯示 Normal vs EPPK
# 安裝並載入必要的套件
if (!require('ggplot2')) { install.packages('ggplot2'); library(ggplot2) }
if (!require('dplyr')) { install.packages('dplyr'); library(dplyr) }
if (!require('tidyr')) { install.packages('tidyr'); library(tidyr) }
if (!require('stringr')) { install.packages('stringr'); library(stringr) }

# Step 1: 繪製 Seurat DotPlot
p2 <- DotPlot(seuratObject_Sample, features = target_genes, 
              group.by = "BroadCellTypeAnnot_SeuratClusters", 
              split.by = "Type", cols="RdBu") +
  theme_minimal() + 
  scale_size_continuous(range = c(0,5)) +  
  ggtitle("Bubble Plot of Target Genes Expression (EPPK vs Normal)")

p2

# Step 2: 取得 `p2` (DotPlot) 的 Y 軸排序
dotplot_data <- DotPlot(seuratObject_Sample, features = target_genes, 
                        group.by = "BroadCellTypeAnnot_SeuratClusters", 
                        split.by = "Type")$data

# 取得 `id` (細胞類型) 的排序順序
id_order <- levels(dotplot_data$id)  # 提取 DotPlot 的 Y 軸順序

# Step 3: 取得 `p2_2` 的資料
expr_data <- dotplot_data

# 確保 `id` 為字串類型
expr_data <- expr_data %>%
  mutate(id = as.character(id))

# Step 4: 提取 `Type` (Normal / EPPK)
expr_data <- expr_data %>%
  mutate(Type = word(id, -1, sep = "_")) %>%  # 取得 `id` 欄位最後一部分
  filter(Type %in% c("Normal", "EPPK"))  # 只保留 Normal / EPPK

# Step 5: 設定 `id` 為 Factor，按照 `p2` 的順序排序
expr_data$id <- factor(expr_data$id, levels = id_order)

# Step 6: 設定 Normal & EPPK 的顏色
custom_colors <- c("Normal" = "#90EE90", "EPPK" = "#FFB6C1")  # Normal = 粉綠，EPPK = 粉紅

# Step 7: 繪製 **直式** 氣泡圖 (Bubble Plot)，並加上黑色外框
p2_2 <- ggplot(expr_data, aes(x = features.plot, y = id, size = pct.exp, fill = Type)) +
  geom_point(shape = 21, stroke = 0.8, color = "black") +  # 加入黑色邊框
  scale_fill_manual(values = custom_colors) +  # 設定 Normal & EPPK 顏色
  scale_size_continuous(range = c(0,5)) +  # 控制氣泡大小
  theme_minimal() + 
  labs(x = "Marker Genes", y = "Cell Type", size = "% Expression", fill = "Sample Type") +
  ggtitle("Bubble Plot of Target Genes Expression (EPPK vs Normal)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # X 軸標籤旋轉以便閱讀

# Step 8: 顯示氣泡圖
p2
p2_2



################################################################################
# 先建立一個布林向量，表示哪些細胞的註解含有 "keratinocyte"
keratinocyte_idx <- grepl("keratinocyte", seuratObject_Sample$BroadCellTypeAnnot_SeuratClusters)

# 使用 [] 索引方式，只保留這些細胞
seuratObject_keratinocyte <- seuratObject_Sample[, keratinocyte_idx]


DimPlot(seuratObject_keratinocyte, group.by = "label_scMRMA_Level4", reduction = "umap", label = T)
DimPlot(seuratObject_keratinocyte, group.by = "BroadCellTypeAnnot_SeuratClusters", reduction = "umap", label = T)

