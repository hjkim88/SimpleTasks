###
#   File name : PID5137_Lineage_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : May 31, 2022
#   Email     : hyunjin.kim2@regeneron.com
#   Purpose   : Run monocle2 on the given single cell dataset.
#
#   Instruction
#               1. Source("PID5137_Lineage_Analysis")
#               2. Run the function "lineage_analysis" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_PID5137_Lineage_Analysis/PID5137_Lineage_Analysis")
#               > lineage_analysis(Seurat_RObj_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/tx_fil_dim20_res250.Robj",
#                                  outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5137/")
###

lineage_analysis <- function(Seurat_RObj_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/tx_fil_dim20_res250.Robj",
                             outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5137/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(devtools, quietly = TRUE)) {
    install.packages("devtools")
    require(devtools, quietly = TRUE)
  }
  if(!require(monocle3, quietly = TRUE)) {
    devtools::install_github('cole-trapnell-lab/monocle3')
    require(monocle, quietly = TRUE)
  }
  
  ### loads an RData file, and returns it
  loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  
  ### load R object
  seurat_obj <- loadRData(Seurat_RObj_path)
  gc()
  
  ### Construct a monocle cds
  cds <- new_cell_data_set(as(seurat_obj@assays$RNA@counts, 'sparseMatrix'),
                           cell_metadata = seurat_obj@meta.data,
                           gene_metadata = data.frame(gene_short_name = row.names(seurat_obj@assays$RNA@counts),
                                                      row.names = row.names(seurat_obj@assays$RNA@counts),
                                                      stringsAsFactors = FALSE, check.names = FALSE))
  
  ### pre-process the data
  cds <- preprocess_cds(cds, num_dim = 50)
  
  ### dimensionality reduction
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  
  ### cluster the data
  cds <- cluster_cells(cds)
  
  ### learn trajectory 
  cds <- learn_graph(cds)
  
  ### plot trajectory
  plot_cells(cds,
             color_cells_by = "seurat_clusters",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  
  target_idx <- intersect(intersect(which(cds@colData$sample.tis == "lung"),
                                    which(cds@colData$sample.hdm == "saline")),
                          intersect(which(cds@colData$sample.cel == "macrophages"),
                                    which(cds@colData$sample.trt == "15 week saline")))
  
  cds@colData$target <- "Non-Target"
  cds@colData$target[target_idx] <- "Target"
  
  plot_cells(cds,
             color_cells_by = "target",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  
  ### cluster 47
  
  
  ###
  cds <- order_cells(cds, root_cells = rownames(cds@colData)[which(cds@colData$seurat_clusters == "47")])
  
  ###
  plot_cells(cds,
             color_cells_by = "pseudotime",
             label_cell_groups=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             graph_label_size=1.5)
  
  ### divide the object and look into those deeply
  
  
  
  
  
  
}