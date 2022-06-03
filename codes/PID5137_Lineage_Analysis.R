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
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
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
  
  ### find all possible partitions
  all_partitions <- unique(cds@clusters$UMAP$partitions)
  all_partitions <- all_partitions[all_partitions != "1"]
  
  ### set all partitions to 1
  cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions %in% all_partitions] <- "1"
  
  ### learn trajectory 
  cds <- learn_graph(cds)
  
  ### plot trajectory
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "seurat_clusters",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_clusters.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "sample.cel",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  label_cell_groups = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_sample_cel.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "sample.tis",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  label_cell_groups = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_sample_tis.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "sample.trt",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  label_cell_groups = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_sample_trt.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "sample.hdm",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  label_cell_groups = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_sample_hdm.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### phenotype that tends to be a starting point
  target_idx <- intersect(intersect(which(cds@colData$sample.tis == "lung"),
                                    which(cds@colData$sample.hdm == "saline")),
                          intersect(which(cds@colData$sample.cel == "macrophages"),
                                    which(cds@colData$sample.trt == "15 week saline")))
  
  ### check where the phenotype is located
  cds@colData$target <- "Non-Target"
  cds@colData$target[target_idx] <- "Target"
  
  plot_cells(cds,
             color_cells_by = "target",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  
  ### based on the plot, cluster 47 seems to be the starting point
  ### order cells based on the cluster 47 cells
  cds <- order_cells(cds, root_cells = rownames(cds@colData)[which(cds@colData$seurat_clusters == "47")])
  
  ### see how the pseudotime looks like
  p <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  label_cell_groups=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE,
                  graph_label_size=1.5)
  ggsave(paste0(outputDir, "Monocle3_trajectory_pseudotime.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### divide the object and look into those deeply
  ### 1. based on sample.cel
  ### 2. based on sample.trt
  
  
  ### check whether the orders are the same
  print(identical(rownames(seurat_obj@meta.data), colnames(seurat_obj@assays$RNA@counts)))
  print(identical(names(Idents(object = seurat_obj)), rownames(seurat_obj@meta.data)))
  
  ### set coloring variables
  col_var <- c("seurat_clusters", "sample.cel", "sample.tis", "sample.trt", "sample.hdm", "pseudotime")
  
  ### 1. based on sample.cel
  for(cel in unique(seurat_obj$sample.cel)) {
    
    ### create new output dir
    outputDir2 <- paste0(outputDir, "/", cel, "/")
    dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
    
    ### subset
    subset_seurat_obj <- subset(seurat_obj,
                                cells = rownames(seurat_obj@meta.data)[which(seurat_obj$sample.cel == cel)])
    
    ### construct a monocle cds
    sub_cds <- new_cell_data_set(as(subset_seurat_obj@assays$RNA@counts, 'sparseMatrix'),
                                 cell_metadata = subset_seurat_obj@meta.data,
                                 gene_metadata = data.frame(gene_short_name = row.names(subset_seurat_obj@assays$RNA@counts),
                                                            row.names = row.names(subset_seurat_obj@assays$RNA@counts),
                                                            stringsAsFactors = FALSE, check.names = FALSE))
    
    ### pre-process the data
    sub_cds <- preprocess_cds(sub_cds, num_dim = 50)
    
    ### dimensionality reduction
    sub_cds <- reduce_dimension(sub_cds, reduction_method = "UMAP")
    
    ### cluster the data
    sub_cds <- cluster_cells(sub_cds)
    
    ### find all possible partitions
    all_partitions <- unique(sub_cds@clusters$UMAP$partitions)
    all_partitions <- all_partitions[all_partitions != "1"]
    
    ### set all partitions to 1
    sub_cds@clusters$UMAP$partitions[sub_cds@clusters$UMAP$partitions %in% all_partitions] <- "1"
    
    ### learn trajectory 
    sub_cds <- learn_graph(sub_cds)
    
    ### pick a most bottom-left cell's cluster
    umap_sum <- apply(sub_cds@reduce_dim_aux$UMAP$model$umap_model$embedding, 1, sum)
    root_cluster <- as.character(sub_cds@colData[names(umap_sum)[which(umap_sum == min(umap_sum))[1]],"seurat_clusters"])
    
    ### order the cells
    sub_cds <- order_cells(sub_cds, root_cells = rownames(sub_cds@colData)[which(sub_cds@colData$seurat_clusters == root_cluster)])
    
    ### set pseudotime info to the seurat object
    subset_seurat_obj$pseudotime <- NA
    subset_seurat_obj@meta.data[names(pseudotime(sub_cds)),"pseudotime"] <- pseudotime(sub_cds)
    
    ### print trajectory & violin plots
    for(var in setdiff(col_var, "sample.cel")) {
      
      ### trajectory plots
      if(var == "seurat_clusters") {
        p <- plot_cells(sub_cds,
                        label_roots = FALSE,
                        color_cells_by = var,
                        label_groups_by_cluster = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = FALSE,
                        group_label_size = 7)
      } else {
        p <- plot_cells(sub_cds,
                        label_roots = FALSE,
                        color_cells_by = var,
                        label_groups_by_cluster = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = FALSE,
                        label_cell_groups = FALSE,
                        group_label_size = 7)
      }
      ggsave(paste0(outputDir2, "Monocle3_trajectory_", cel, "_subset_", var, ".pdf"), plot = p, width = 15, height = 10, dpi = 350)
      
      ### violin plots
      if((var != "seurat_clusters") && (var != "pseudotime")) {
        ### set ident with the persistency info
        subset_seurat_obj <- SetIdent(object = subset_seurat_obj,
                                      cells = rownames(subset_seurat_obj@meta.data),
                                      value = subset_seurat_obj@meta.data[,var])
        
        ### draw violin plot
        p <- VlnPlot(subset_seurat_obj, features = "pseudotime",
                     pt.size = 0) +
                theme_classic(base_size = 40) +
                theme(legend.position = "none",
                      legend.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
                      legend.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
                      plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35, color = "black", face = "bold"),
                      plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, size = 25, color = "black", face = "bold"),
                      axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
                      axis.title.x = element_blank(),
                      axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
        
        ### save the violin plot
        ggsave(file = paste0(outputDir2, "Monocle3_pseudotime_box_", cel, "_subset_", var, ".pdf"), plot = p, width = 15, height = 10, dpi = 350)
      }
      
    }
  }
  
  
  ### 2. based on sample.trt
  for(trt in unique(seurat_obj$sample.trt)) {
    
    ### create new output dir
    outputDir2 <- paste0(outputDir, "/", trt, "/")
    dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
    
    ### subset
    subset_seurat_obj <- subset(seurat_obj,
                                cells = rownames(seurat_obj@meta.data)[which(seurat_obj$sample.trt == trt)])
    
    ### construct a monocle cds
    sub_cds <- new_cell_data_set(as(subset_seurat_obj@assays$RNA@counts, 'sparseMatrix'),
                                 cell_metadata = subset_seurat_obj@meta.data,
                                 gene_metadata = data.frame(gene_short_name = row.names(subset_seurat_obj@assays$RNA@counts),
                                                            row.names = row.names(subset_seurat_obj@assays$RNA@counts),
                                                            stringsAsFactors = FALSE, check.names = FALSE))
    
    ### pre-process the data
    sub_cds <- preprocess_cds(sub_cds, num_dim = 50)
    
    ### dimensionality reduction
    sub_cds <- reduce_dimension(sub_cds, reduction_method = "UMAP")
    
    ### cluster the data
    sub_cds <- cluster_cells(sub_cds)
    
    ### find all possible partitions
    all_partitions <- unique(sub_cds@clusters$UMAP$partitions)
    all_partitions <- all_partitions[all_partitions != "1"]
    
    ### set all partitions to 1
    sub_cds@clusters$UMAP$partitions[sub_cds@clusters$UMAP$partitions %in% all_partitions] <- "1"
    
    ### learn trajectory 
    sub_cds <- learn_graph(sub_cds)
    
    ### pick a most bottom-left cell's cluster
    umap_sum <- apply(sub_cds@reduce_dim_aux$UMAP$model$umap_model$embedding, 1, sum)
    root_cluster <- as.character(sub_cds@colData[names(umap_sum)[which(umap_sum == min(umap_sum))[1]],"seurat_clusters"])
    
    ### order the cells
    sub_cds <- order_cells(sub_cds, root_cells = rownames(sub_cds@colData)[which(sub_cds@colData$seurat_clusters == root_cluster)])
    
    ### set pseudotime info to the seurat object
    subset_seurat_obj$pseudotime <- NA
    subset_seurat_obj@meta.data[names(pseudotime(sub_cds)),"pseudotime"] <- pseudotime(sub_cds)
    
    ### print trajectory & violin plots
    for(var in setdiff(col_var, "sample.trt")) {
      
      ### trajectory plots
      if(var == "seurat_clusters") {
        p <- plot_cells(sub_cds,
                        label_roots = FALSE,
                        color_cells_by = var,
                        label_groups_by_cluster = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = FALSE,
                        group_label_size = 7)
      } else {
        p <- plot_cells(sub_cds,
                        label_roots = FALSE,
                        color_cells_by = var,
                        label_groups_by_cluster = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = FALSE,
                        label_cell_groups = FALSE,
                        group_label_size = 7)
      }
      ggsave(paste0(outputDir2, "Monocle3_trajectory_", trt, "_subset_", var, ".pdf"), plot = p, width = 15, height = 10, dpi = 350)
      
      ### violin plots
      if((var != "seurat_clusters") && (var != "pseudotime")) {
        ### set ident with the persistency info
        subset_seurat_obj <- SetIdent(object = subset_seurat_obj,
                                      cells = rownames(subset_seurat_obj@meta.data),
                                      value = subset_seurat_obj@meta.data[,var])
        
        ### draw violin plot
        p <- VlnPlot(subset_seurat_obj, features = "pseudotime",
                     pt.size = 0) +
          theme_classic(base_size = 40) +
          theme(legend.position = "none",
                legend.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
                legend.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
                plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35, color = "black", face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, size = 25, color = "black", face = "bold"),
                axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
                axis.title.x = element_blank(),
                axis.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
        
        ### save the violin plot
        ggsave(file = paste0(outputDir2, "Monocle3_pseudotime_box_", trt, "_subset_", var, ".pdf"), plot = p, width = 15, height = 10, dpi = 350)
      }
      
    }
  }
  
  
  
}