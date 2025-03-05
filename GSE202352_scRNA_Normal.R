##### Presetting #####
rm(list = ls()) 
memory.limit(150000)

#### Load Packages ####
if(!require('Seurat')) {install.packages('Seurat'); library(Seurat)}
if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}
if(!require('stringr')) {install.packages('stringr'); library(stringr)}
if(!require('R.utils')) {install.packages('R.utils'); library(R.utils)}

data_dir <- "C:/Charlene/Dataset_Online/Healthy_Skin/GSE202352_RAW_Ori"
tar_files <- list.files(data_dir, pattern = "\\.tar$|\\.tar\\.gz$", full.names = TRUE)
if(length(tar_files) == 0) {
  stop("No .tar or .tar.gz files found. Check directory or file extension.")
}


#### 1) Read each tarball ####
extract_metadata <- function(filename) {
  base_name <- basename(filename)
  # Remove .tar/.tar.gz plus the 'filtered_feature_bc_matrix' suffix, if present
  base_name <- sub("\\.filtered_feature_bc_matrix\\.tar\\.gz$","",base_name)
  base_name <- sub("\\.filtered_feature_bc_matrix\\.tar$","",base_name)
  
  parts <- unlist(str_split(base_name, "_"))
  list(
    FullName = base_name,
    GSM_ID = parts[1],
    Sample_Subject = paste(parts[2], parts[3], parts[4], sep = "_"),
    Tissue = parts[5]
  )
}


seurat_list <- list()

for (file in tar_files) {
  metadata <- extract_metadata(file)
  
  # Create a folder to extract into
  extract_dir <- file.path(data_dir, metadata$FullName)
  dir.create(extract_dir, showWarnings = FALSE)
  
  # Untar
  untar(file, exdir = extract_dir)
  
  # After extraction, we may have one of these scenarios:
  # 1) matrix.mtx(.gz) etc. are directly in extract_dir
  # 2) The files are inside a subfolder named "filtered_feature_bc_matrix"
  # 3) The files might be uncompressed rather than .gz
  
  # Try to find either extract_dir OR a subfolder named 'filtered_feature_bc_matrix'
  possible_dirs <- c(extract_dir, file.path(extract_dir, "filtered_feature_bc_matrix"))
  data_found <- FALSE
  
  for(dsub in possible_dirs) {
    if(dir.exists(dsub)) {
      # Check for presence of any matrix .mtx(.gz) + features + barcodes
      # If "matrix.mtx" or "matrix.mtx.gz" is in dsub, we let Read10X handle it.
      # We simply attempt to read; if it fails, we move on.
      tryCatch({
        message("Attempting to Read10X from: ", dsub)
        counts_data <- Read10X(data.dir = dsub)
        
        # If successful, create Seurat object
        seurat_obj <- CreateSeuratObject(
          counts = counts_data,
          # project = "GSE202352"  # uniform short project name
        )
        
        # Add minimal metadata
        for (nm in names(metadata)) {
          seurat_obj[[nm]] <- metadata[[nm]]
          # seurat_obj@meta.data[[nm]] <- metadata[[nm]]
        }
        seurat_obj <- RenameCells(seurat_obj, add.cell.id = metadata$GSM_ID)
        
        # Store
        seurat_list[[metadata$FullName]] <- seurat_obj
        
        data_found <- TRUE
        break  # we can break out of the loop once successful
      }, error = function(e) {
        # do nothing, just try the next possibility
      })
    }
    if(data_found) break
  } # end for(dsub in possible_dirs)
  
  if(!data_found) {
    warning("Could not find valid matrix/barcodes/features under ", file)
    next
  }
}

#### 2) Normalize and find variable features for each object ####
seurat_list <- lapply(seurat_list, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  return(obj)
})

saveRDS(seurat_list, file = "GSE202352_combined_seurat_list.rds")

# 批次效應整合
if(length(seurat_list) > 0) {
  # 使用 FindIntegrationAnchors 和 IntegrateData 來進行批次效應處理
  features <- SelectIntegrationFeatures(object.list = seurat_list)
  anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
  combined_seurat <- IntegrateData(anchorset = anchors)
  
  # 正規化和降維
  DefaultAssay(combined_seurat) <- "integrated"
  combined_seurat <- ScaleData(combined_seurat)
  combined_seurat <- RunPCA(combined_seurat)
  combined_seurat <- FindNeighbors(combined_seurat, dims = 1:10)
  combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)
  combined_seurat <- RunUMAP(combined_seurat, dims = 1:10)
  
  # 更新 Tissue 資訊
  combined_seurat$Tissue <- ifelse(
    combined_seurat$GSM_ID %in% c("GSM6111848", "GSM6111852"), "hip_dermis",
    ifelse(
      combined_seurat$GSM_ID %in% c("GSM6111849", "GSM6111853"), "palm_dermis",
      ifelse(
        combined_seurat$GSM_ID %in% c("GSM6111850", "GSM6111854"), "hip_epi",
        ifelse(combined_seurat$GSM_ID %in% c("GSM6111851", "GSM6111855"), "palm_epi", combined_seurat$Tissue)
      )
    )
  )
  
  # 保存整合後的 Seurat 物件
  saveRDS(combined_seurat, file = "GSE202352_combined_seurat_integrated.rds")
  
} else {
  stop("No valid Seurat objects were created. Please check file structure and contents.")
}

# 後續分析
DefaultAssay(combined_seurat) <- "RNA"

Marker.lt <- c("EPGN", "EGFR", "KRT1", "KRT9", 
               "KIT", "EP300", "NF1", "COL1A1")
VlnPlot(combined_seurat, features = Marker.lt, ncol = 4)
# VlnPlot(combined_seurat, features = Marker.lt[1:4], ncol = 2, group.by = "blueprint.main")
# VlnPlot(combined_seurat, features = Marker.lt[5:8], ncol = 2, group.by = "blueprint.main")
FeaturePlot(combined_seurat, features = Marker.lt, ncol = 4)

DimPlot(combined_seurat, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "GSM_ID") + NoLegend()
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "GSM_ID")
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "Tissue" ) 
DimPlot(combined_seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = "seurat_clusters" ) 
