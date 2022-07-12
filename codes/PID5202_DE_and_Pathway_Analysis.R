###
#   File name : PID5202_DE_and_Pathway_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : Jul 7, 2022
#   Email     : hyunjin.kim2@regeneron.com
#   Purpose   : Run DE and Pathway analysis on PID5202 data (convert h5ad to seurat object first)
#               The h5ad was generated from Cooper's scPIPE pre-processing
#
#   Instruction
#               1. Source("PID5202_DE_and_Pathway_Analysis.R")
#               2. Run the function "some_analyses" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_PID5202_DE_and_Pathway_Analysis.R/PID5202_DE_and_Pathway_Analysis.R")
#               > some_analyses(h5ad_file_path="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_gex.h5ad",
#                               Seurat_RObj_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/PID5202/PID5202_concat_gex.rds",
#                               outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/")
###

some_analyses <- function(h5ad_file_path="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_gex.h5ad",
                          Seurat_RObj_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/PID5202/PID5202_concat_gex.rds",
                          outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(SeuratDisk, quietly = TRUE)) {
    remotes::install_github("mojaveazure/seurat-disk")
    require(SeuratDisk, quietly = TRUE)
  }
  if(!require(data.table, quietly = TRUE)) {
    install.packages("data.table")
    require(data.table, quietly = TRUE)
  }
  
  ### convert h5ad to h5seurat and load
  SeuratDisk::Convert(source = h5ad_file_path, dest = "h5seurat", overwrite = TRUE)
  
  ### load the converted h5ad object
  seurat_obj <- LoadH5Seurat(paste0(substr(h5ad_file_path, 1, nchar(h5ad_file_path)-5), ".h5seurat"), misc = FALSE)
  
  ### because the process above takes some memory I ran it on JupyterHub
  seurat_obj <- readRDS("/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_gex_new.rds")
  
  ### load the meta data
  meta_data <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_metadata.csv",
                     header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)

  ### change data.table to data.frame
  meta_data <- setDF(meta_data)

  ### set gene names as row names
  rownames(meta_data) <- meta_data[,1]
  meta_data <- meta_data[,-1]
  
  ### attach metadata
  seurat_obj@meta.data <- meta_data
  
  ### load the raw counts
  count_mat <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_raw_counts.csv",
                     header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)

  ### change data.table to data.frame
  count_mat <- setDF(count_mat)

  ### set gene names as row names
  rownames(count_mat) <- count_mat[,1]
  count_mat <- count_mat[,-1]

  ### transverse the count matrix (rows: genes, columns: cells)
  count_mat <- t(count_mat)
  
  ### check the seurat_obj@assays$RNA@counts and count_mat have the same values
  print(as.vector(as.matrix(seurat_obj@assays$RNA@counts[1:5, 1:5])))
  print(as.vector(as.matrix(count_mat[1:5, 1:5])))
  print(as.vector(as.matrix(seurat_obj@assays$RNA@counts[1000:1005, 2000:2005])))
  print(as.vector(as.matrix(count_mat[1000:1005, 2000:2005])))
  print(as.character(seurat_obj@assays$RNA@counts[5, 3]))
  print(as.character(count_mat[5, 3]))
  
  ### remove loaded count_mat
  rowNames <- rownames(count_mat)
  colNames <- colnames(count_mat)
  rm(count_mat)
  rm(meta_data)
  gc()
  
  ### copy gene names and the cell names
  ### checked manually that they have the same values
  rownames(seurat_obj@assays$RNA@counts) <- rowNames
  colnames(seurat_obj@assays$RNA@counts) <- colNames
  rownames(seurat_obj@assays$RNA@data) <- rowNames
  colnames(seurat_obj@assays$RNA@data) <- colNames
  
  ### now load the scaled data
  scaled_mat <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_scaled_counts.csv",
                      header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)

  ## change data.table to data.frame
  scaled_mat <- setDF(scaled_mat)

  ### set gene names as row names
  rownames(scaled_mat) <- scaled_mat[,1]
  scaled_mat <- scaled_mat[,-1]

  ### transverse the count matrix (rows: genes, columns: cells)
  scaled_mat <- t(scaled_mat)
  
  ### check the seurat_obj@assays$RNA@scale.data and scaled_mat have the same values
  print(as.vector(as.matrix(seurat_obj@assays$RNA@scale.data[1:5, 1:5])))
  print(as.vector(as.matrix(scaled_mat[1:5, 1:5])))
  print(as.vector(as.matrix(seurat_obj@assays$RNA@scale.data[1000:1005, 1000:1005])))
  print(as.vector(as.matrix(scaled_mat[1000:1005, 1000:1005])))
  print(as.character(seurat_obj@assays$RNA@scale.data[5, 5]))
  print(as.character(scaled_mat[5, 5]))
  
  ### remove loaded scaled_mat
  rowNames <- rownames(scaled_mat)
  colNames <- colnames(scaled_mat)
  rm(scaled_mat)
  gc()
  
  ### copy gene names and the cell names
  ### checked manually that they have the same values
  rownames(seurat_obj@assays$RNA@scale.data) <- rowNames
  colnames(seurat_obj@assays$RNA@scale.data) <- colNames
  
  ### rownames in the meta.data should be in the same order as colnames in the counts
  print(identical(rownames(seurat_obj@assays$RNA@counts), rownames(seurat_obj@assays$RNA@data)))
  print(identical(rownames(seurat_obj@assays$RNA@counts), rownames(seurat_obj@assays$RNA@scale.data)))
  print(identical(colnames(seurat_obj@assays$RNA@counts), colnames(seurat_obj@assays$RNA@data)))
  print(identical(colnames(seurat_obj@assays$RNA@counts), colnames(seurat_obj@assays$RNA@scale.data)))
  seurat_obj@meta.data <- seurat_obj@meta.data[colnames(seurat_obj@assays$RNA@counts),]
  print(identical(rownames(seurat_obj@meta.data), colnames(seurat_obj@assays$RNA@counts)))
  
  ### garbage collection
  gc()
  
  ### load PCA/UMAP
  pca_data <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_pca.csv",
                    header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE)
  pca_data2 <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_pca_harmony.csv",
                     header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
  umap_data <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_umap.csv",
                     header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
  umap_data2 <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_umap_raw.csv",
                      header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)

  ### change data.table to data.frame
  pca_data <- setDF(pca_data)
  pca_data2 <- setDF(pca_data2)
  umap_data <- setDF(umap_data)
  umap_data2 <- setDF(umap_data2)

  ### set gene names as row names
  rownames(pca_data) <- pca_data[,1]
  pca_data <- pca_data[,-1]
  rownames(pca_data2) <- pca_data2[,1]
  pca_data2 <- pca_data2[,-1]
  rownames(umap_data) <- umap_data[,1]
  umap_data <- umap_data[,-1]
  rownames(umap_data2) <- umap_data2[,1]
  umap_data2 <- umap_data2[,-1]
  
  ### check the seurat_obj@assays$RNA@scale.data and scaled_mat have the same values
  print(identical(seurat_obj@reductions$umap@cell.embeddings[1:5,1:2], umap_data[1:5,1:2]))
  print(seurat_obj@reductions$umap@cell.embeddings[1:5,1:2])
  print(umap_data[1:5,1:2])
  print(seurat_obj@reductions$pca@cell.embeddings[100:105,1:2])
  print(pca_data[100:105,1:2])
  print(seurat_obj@reductions$pca_harmony@cell.embeddings[100:105,1:2])
  print(pca_data2[100:105,1:2])
  print(seurat_obj@reductions$umap_raw@cell.embeddings[200:205,1:2])
  print(umap_data2[200:205,1:2])
  
  ### copy cell names
  rownames(seurat_obj@reductions$pca@cell.embeddings) <- rownames(pca_data)
  rownames(seurat_obj@reductions$pca_harmony@cell.embeddings) <- rownames(pca_data2)
  rownames(seurat_obj@reductions$umap@cell.embeddings) <- rownames(umap_data)
  rownames(seurat_obj@reductions$umap_raw@cell.embeddings) <- rownames(umap_data2)
  
  ### remove dim objects
  rm(pca_data)
  rm(pca_data2)
  rm(umap_data)
  rm(umap_data2)
  gc()
  
  ### attach new annotation
  ### time, treatment, MSR1+/-
  seurat_obj$time <- "NO_LPS"
  seurat_obj$time[grep(pattern = "2h", seurat_obj$meta_tissue)] <- "2HR_PRIOR_TO_LPS"
  seurat_obj$time[grep(pattern = "24h", seurat_obj$meta_tissue)] <- "24HR_PRIOR_TO_LPS"
  
  seurat_obj$treatment <- "PBS"
  seurat_obj$treatment[grep(pattern = "M483", seurat_obj$meta_tissue)] <- "Steroid"
  seurat_obj$treatment[grep(pattern = "REGN3892", seurat_obj$meta_tissue)] <- "Ctrl_ncADC"
  seurat_obj$treatment[grep(pattern = "REGN4323", seurat_obj$meta_tissue)] <- "MSR1_ncADC"
  
  plot(density(seurat_obj@assays$RNA@scale.data["Msr1",]),
       main = "Scaled Expression of Msr1",
       xlab = "Scaled Expression")
  
  seurat_obj$Msr1 <- "Msr1_Neg"
  seurat_obj@meta.data[colnames(seurat_obj@assays$RNA@scale.data)[which(seurat_obj@assays$RNA@scale.data["Msr1",] > 0)],"Msr1"] <- "Msr1_Pos"
  
  ### test the umap plot
  DimPlot(seurat_obj,
          group.by = "meta_tissue",
          raster = TRUE)
  
  ### save the new seurat object
  saveRDS(seurat_obj,
          file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_gex_new.rds")
  
  
  
  
  # ### load the raw counts
  # count_mat <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_raw_counts.csv",
  #                    header = TRUE, 
  #                    stringsAsFactors = FALSE, check.names = FALSE)
  # 
  # ### change data.table to data.frame
  # count_mat <- setDF(count_mat)
  # 
  # ### set gene names as row names
  # rownames(count_mat) <- count_mat[,1]
  # count_mat <- count_mat[,-1]
  # 
  # ### transverse the count matrix (rows: genes, columns: cells)
  # count_mat <- t(count_mat)
  # 
  # ### load the meta data
  # meta_data <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_metadata.csv",
  #                    header = TRUE, 
  #                    stringsAsFactors = FALSE, check.names = FALSE)
  # 
  # ### change data.table to data.frame
  # meta_data <- setDF(meta_data)
  # 
  # ### set gene names as row names
  # rownames(meta_data) <- meta_data[,1]
  # meta_data <- meta_data[,-1]
  # 
  # ### garbage collection
  # gc()
  # 
  # ### create seurat object
  # seurat_obj <- CreateSeuratObject(count_mat, project = "PID5202", assay = "RNA",
  #                                  min.cells = 0, min.features = 0, names.field = 1,
  #                                  names.delim = "_", meta.data = meta_data)
  # 
  # ### active assay = "RNA"
  # seurat_obj@active.assay <- "RNA"
  # 
  # ### remove trash items and garbage collection
  # rm(count_mat)
  # rm(meta_data)
  # gc()
  # 
  # ### now load the scaled data
  # scaled_mat <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_scaled_counts.csv",
  #                     header = TRUE, 
  #                     stringsAsFactors = FALSE, check.names = FALSE)
  # 
  # ## change data.table to data.frame
  # scaled_mat <- setDF(scaled_mat)
  # 
  # ### set gene names as row names
  # rownames(scaled_mat) <- scaled_mat[,1]
  # scaled_mat <- scaled_mat[,-1]
  # 
  # ### transverse the count matrix (rows: genes, columns: cells)
  # scaled_mat <- t(scaled_mat)
  # 
  # ### set scaled data to the seurat object
  # seurat_obj@assays$RNA@scale.data <- scaled_mat
  # 
  # ### remove trash items and garbage collection
  # rm(scaled_mat)
  # gc()
  # 
  # ### load PCA/UMAP
  # pca_data <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_pca.csv",
  #                   header = TRUE, 
  #                   stringsAsFactors = FALSE, check.names = FALSE)
  # pca_data2 <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_pca_harmony.csv",
  #                    header = TRUE, 
  #                    stringsAsFactors = FALSE, check.names = FALSE)
  # umap_data <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_umap.csv",
  #                    header = TRUE, 
  #                    stringsAsFactors = FALSE, check.names = FALSE)
  # umap_data2 <- fread(file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_umap_raw.csv",
  #                     header = TRUE, 
  #                     stringsAsFactors = FALSE, check.names = FALSE)
  # 
  # ### change data.table to data.frame
  # pca_data <- setDF(pca_data)
  # pca_data2 <- setDF(pca_data2)
  # umap_data <- setDF(umap_data)
  # umap_data2 <- setDF(umap_data2)
  # 
  # ### set gene names as row names
  # rownames(pca_data) <- pca_data[,1]
  # pca_data <- pca_data[,-1]
  # rownames(pca_data2) <- pca_data2[,1]
  # pca_data2 <- pca_data2[,-1]
  # rownames(umap_data) <- umap_data[,1]
  # umap_data <- umap_data[,-1]
  # rownames(umap_data2) <- umap_data2[,1]
  # umap_data2 <- umap_data2[,-1]
  # 
  # ### garbage collection
  # gc()
  
  
  ### the converting with SeuratDisk performed on JupyterHub
  ### but metadata is not added
  ### so we load the partial output and attach the metadata
  ### line 45 of this script
  
  
  
}
