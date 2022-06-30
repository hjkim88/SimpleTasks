###
#   File name : PID4604_Lineage_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : May 31, 2022
#   Email     : hyunjin.kim2@regeneron.com
#   Purpose   : Run monocle3 on the given single cell dataset.
#
#   Instruction
#               1. Source("PID4604_Lineage_Analysis.R")
#               2. Run the function "lineage_analysis" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_PID4604_Lineage_Analysis.R/PID4604_Lineage_Analysis.R")
#               > lineage_analysis(Seurat_RObj_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/tx_fil_dim20_res250.Robj",
#                                  outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID4604/")
###

lineage_analysis <- function(Seurat_RObj_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/tx_fil_dim20_res250.Robj",
                             outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID4604/") {
  
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
    require(monocle3, quietly = TRUE)
  }
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
  }
  if(!require(igraph, quietly = TRUE)) {
    install.packages("igraph")
    require(igraph, quietly = TRUE)
  }
  if(!require(philentropy, quietly = TRUE)) {
    install.packages("philentropy")
    require(philentropy, quietly = TRUE)
  }
  if(!require(rgdal, quietly = TRUE)) {
    install.packages("rgdal")
    require(rgdal, quietly = TRUE)
  }
  
  ### loads an RData file, and returns it
  loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  
  ### load R object
  seurat_obj <- loadRData(Seurat_RObj_path)
  gc()
  
  ### load the new annotation info
  ### id - cluster#, label - annotation info
  new_anno <- read.table(file = "./data/cluster_annot.txt", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(new_anno) <- new_anno$id
  
  ### annotate
  seurat_obj$new_anno <- NA
  for(clstr in unique(as.character(seurat_obj$seurat_clusters))) {
    seurat_obj$new_anno[which(seurat_obj$seurat_clusters == clstr)] <- new_anno[clstr,"label"]
  }
  
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
  
  
  ### using existing UMAP cords for the trajectory analysis
  
  ### overwrite the current CDS UMAP with original seurat UMAP
  cds@int_colData@listData$reducedDims$UMAP <- seurat_obj@reductions$umap@cell.embeddings
  
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
  ggsave(paste0(outputDir, "Monocle3_trajectory_clusters_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "sample.cel",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  label_cell_groups = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_sample_cel_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "sample.tis",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  label_cell_groups = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_sample_tis_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "sample.trt",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  label_cell_groups = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_sample_trt_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "sample.hdm",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  label_cell_groups = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_sample_hdm_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### phenotype that tends to be a starting point
  target_idx <- intersect(intersect(which(cds@colData$sample.tis == "blood"),
                                    which(cds@colData$sample.hdm == "saline")),
                          intersect(which(cds@colData$sample.cel == "monocytes"),
                                    which(cds@colData$sample.trt == "15 week saline")))
  
  ### check where the phenotype is located
  cds@colData$target <- "Non-Target"
  cds@colData$target[target_idx] <- "Target"
  
  plot_cells(cds,
             color_cells_by = "target",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  
  ### based on the plot, cluster 5 seems to be the starting point
  ### order cells based on the cluster 5 cells
  cds <- order_cells(cds, root_cells = rownames(cds@colData)[which(cds@colData$seurat_clusters == "5")])
  
  ### see how the pseudotime looks like
  p <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  label_cell_groups=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE,
                  label_roots = FALSE,
                  trajectory_graph_color = "green",
                  trajectory_graph_segment_size = 1.5)
  ggsave(paste0(outputDir, "Monocle3_trajectory_pseudotime_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  #
  ### MONOCYTES -> DC2 SUB CLUSTERS
  #
  
  ### check whether the orders are the same
  print(identical(rownames(seurat_obj@meta.data), colnames(seurat_obj@assays$RNA@counts)))
  print(identical(names(Idents(object = seurat_obj)), rownames(seurat_obj@meta.data)))
  
  ### subset Mgl2+ cells
  subset_seurat_obj <- subset(seurat_obj,
                              cells = colnames(seurat_obj@assays$RNA@counts)[which(seurat_obj@assays$RNA@counts["Mgl2",] > 10)])
  
  ### check the UMAP
  DimPlot(object = subset_seurat_obj, reduction = "umap",
          group.by = "seurat_clusters",
          pt.size = 1,
          label = TRUE) +
    theme(legend.position = "none")
  
  ### construct a monocle cds
  sub_cds <- new_cell_data_set(as(subset_seurat_obj@assays$RNA@counts, 'sparseMatrix'),
                               cell_metadata = subset_seurat_obj@meta.data,
                               gene_metadata = data.frame(gene_short_name = row.names(subset_seurat_obj@assays$RNA@counts),
                                                          row.names = row.names(subset_seurat_obj@assays$RNA@counts),
                                                          stringsAsFactors = FALSE, check.names = FALSE))
  
  ### pre-process the data
  sub_cds <- preprocess_cds(sub_cds, num_dim = 50)
  
  ### overwrite the current CDS UMAP with original seurat UMAP
  sub_cds@int_colData@listData$reducedDims$UMAP <- subset_seurat_obj@reductions$umap@cell.embeddings
  
  ### cluster the data
  sub_cds <- cluster_cells(sub_cds,
                           k = 100,
                           partition_qval = 0.01)
  
  ### learn trajectory 
  sub_cds <- learn_graph(sub_cds)
  
  ### plot trajectory
  p <- plot_cells(sub_cds,
                  label_roots = FALSE,
                  color_cells_by = "seurat_clusters",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  group_label_size = 0,
                  cell_size = 0.75)
  ggsave(paste0(outputDir, "Monocle3_trajectory_Mgl2_subset_clusters_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### now divide into saline, 4wk HDM, and 15wk HDM
  ### and do the same analysis
  
  for(trt in unique(subset_seurat_obj$sample.trt)) {
    
    ### subset
    subset_seurat_obj2 <- subset(subset_seurat_obj,
                                 cells = rownames(subset_seurat_obj@meta.data)[which(subset_seurat_obj$sample.trt == trt)])
    
    ### construct a monocle cds
    sub_cds2 <- new_cell_data_set(as(subset_seurat_obj2@assays$RNA@counts, 'sparseMatrix'),
                                  cell_metadata = subset_seurat_obj2@meta.data,
                                  gene_metadata = data.frame(gene_short_name = row.names(subset_seurat_obj2@assays$RNA@counts),
                                                            row.names = row.names(subset_seurat_obj2@assays$RNA@counts),
                                                            stringsAsFactors = FALSE, check.names = FALSE))
    
    ### pre-process the data
    sub_cds2 <- preprocess_cds(sub_cds2, num_dim = 50)
    
    ### overwrite the current CDS UMAP with original seurat UMAP
    sub_cds2@int_colData@listData$reducedDims$UMAP <- subset_seurat_obj2@reductions$umap@cell.embeddings
    
    ### cluster the data
    sub_cds2 <- cluster_cells(sub_cds2,
                              k = 50,
                              partition_qval = 0.01)
    
    ### learn trajectory 
    sub_cds2 <- learn_graph(sub_cds2)
    
    ### plot trajectory
    p <- plot_cells(sub_cds2,
                    label_roots = FALSE,
                    color_cells_by = "seurat_clusters",
                    label_groups_by_cluster = FALSE,
                    label_leaves = FALSE,
                    label_branch_points = FALSE,
                    group_label_size = 0,
                    cell_size = 0.75)
    ggsave(paste0(outputDir, "Monocle3_trajectory_Mgl2_", trt, "_subset_clusters_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
    
  }
  
  
  #
  ### macrophage trajectory
  #
  
  ### subset lung cells
  subset_seurat_obj <- subset(seurat_obj,
                              cells = rownames(seurat_obj@meta.data)[which(seurat_obj$sample.cel == "macrophages")])
  
  ### check the UMAP
  DimPlot(object = subset_seurat_obj, reduction = "umap",
          group.by = "seurat_clusters",
          pt.size = 1,
          label = TRUE) +
    theme(legend.position = "none")
  
  ### construct a monocle cds
  sub_cds <- new_cell_data_set(as(subset_seurat_obj@assays$RNA@counts, 'sparseMatrix'),
                               cell_metadata = subset_seurat_obj@meta.data,
                               gene_metadata = data.frame(gene_short_name = row.names(subset_seurat_obj@assays$RNA@counts),
                                                          row.names = row.names(subset_seurat_obj@assays$RNA@counts),
                                                          stringsAsFactors = FALSE, check.names = FALSE))
  
  ### pre-process the data
  sub_cds <- preprocess_cds(sub_cds, num_dim = 50)
  
  ### overwrite the current CDS UMAP with original seurat UMAP
  sub_cds@int_colData@listData$reducedDims$UMAP <- subset_seurat_obj@reductions$umap@cell.embeddings
  
  ### cluster the data
  sub_cds <- cluster_cells(sub_cds,
                           k = 300,
                           partition_qval = 0.01)
  
  ### learn trajectory 
  sub_cds <- learn_graph(sub_cds)
  
  ### plot trajectory
  p <- plot_cells(sub_cds,
                  label_roots = FALSE,
                  color_cells_by = "seurat_clusters",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  group_label_size = 0,
                  cell_size = 0.75)
  ggsave(paste0(outputDir, "Monocle3_trajectory_Macrophage_subset_clusters_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### now divide into saline, 4wk HDM, and 15wk HDM
  ### and do the same analysis
  
  for(trt in unique(subset_seurat_obj$sample.trt)) {
    
    ### subset
    subset_seurat_obj2 <- subset(subset_seurat_obj,
                                 cells = rownames(subset_seurat_obj@meta.data)[which(subset_seurat_obj$sample.trt == trt)])
    
    ### construct a monocle cds
    sub_cds2 <- new_cell_data_set(as(subset_seurat_obj2@assays$RNA@counts, 'sparseMatrix'),
                                  cell_metadata = subset_seurat_obj2@meta.data,
                                  gene_metadata = data.frame(gene_short_name = row.names(subset_seurat_obj2@assays$RNA@counts),
                                                             row.names = row.names(subset_seurat_obj2@assays$RNA@counts),
                                                             stringsAsFactors = FALSE, check.names = FALSE))
    
    ### pre-process the data
    sub_cds2 <- preprocess_cds(sub_cds2, num_dim = 50)
    
    ### overwrite the current CDS UMAP with original seurat UMAP
    sub_cds2@int_colData@listData$reducedDims$UMAP <- subset_seurat_obj2@reductions$umap@cell.embeddings
    
    ### cluster the data
    sub_cds2 <- cluster_cells(sub_cds2,
                              k = 100,
                              partition_qval = 0.01)
    
    ### learn trajectory 
    sub_cds2 <- learn_graph(sub_cds2)
    
    ### plot trajectory
    p <- plot_cells(sub_cds2,
                    label_roots = FALSE,
                    color_cells_by = "seurat_clusters",
                    label_groups_by_cluster = FALSE,
                    label_leaves = FALSE,
                    label_branch_points = FALSE,
                    group_label_size = 0,
                    cell_size = 0.75)
    ggsave(paste0(outputDir, "Monocle3_trajectory_Macrophage_", trt, "_subset_clusters_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
    
  }
  
  
  ### load the new annotation info
  ### id - cluster#, label - annotation info
  new_anno <- read.table(file = "./data/cluster_annot.txt", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
  rownames(new_anno) <- new_anno$id
  
  ### annotate
  seurat_obj$new_anno <- NA
  for(clstr in unique(as.character(seurat_obj$seurat_clusters))) {
    seurat_obj$new_anno[which(seurat_obj$seurat_clusters == clstr)] <- new_anno[clstr,"label"]
  }
  
  ### check the order of cds@colData and seurat_obj@meta.data are the same
  identical(as.character(cds@colData$seurat_clusters), as.character(seurat_obj$seurat_clusters))
  identical(as.character(cds@colData$sample.cel), as.character(seurat_obj$sample.cel))
  identical(as.character(cds@colData$sample.trt), as.character(seurat_obj$sample.trt))
  
  ### if true, attach the new_anno to the cds object
  cds@colData$new_anno <- as.character(seurat_obj$new_anno)
  
  ### plot trajectory
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "new_anno",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_annotation_ORIGINAL_UMAP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### tree representation of the overall trajectory
  ### root cluster = 5
  
  ###
  ### Make a manual tree visualization function from Monocle3 object
  ###
  ### obj                 : Monocle3 object
  ### trace_base          : Basis of the trajectory; Should be one of colnames(obj@colData)
  ### root                : Root state of the trajectory; Should be one of unique(obj@colData[,trace_base])
  ### reduced_dim_method  : Dimensionality reduction method of the Monocle3; Should be either UMAP or PCA
  ### color_scheme        : Color scheme for the trace_base; If NULL, default colors will be used
  ### node_shape          : Shape of the node in the tree graph
  ### result_path         : The path for the result file; Should be ended with .pdf (pdf file)
  ###
  plot_trajectory_as_tree <- function(obj,
                                      trace_base,
                                      root,
                                      reduced_dim_method = c("UMAP", "PCA"),
                                      color_scheme = NULL,
                                      node_shape = c("circle", "rectangle"),
                                      result_path = "./trajectory_tree.pdf") {
    
    ### set igraph graph
    if(as.character(class(obj)) == "cell_data_set") {
      igraph_grp <- obj@principal_graph[[reduced_dim_method[1]]]
      original_adj_mat <- as.data.frame(obj@principal_graph_aux[[reduced_dim_method[1]]]$stree)
      rownames(original_adj_mat) <- colnames(obj@principal_graph_aux[[reduced_dim_method[1]]]$dp_mst)
      colnames(original_adj_mat) <- colnames(obj@principal_graph_aux[[reduced_dim_method[1]]]$dp_mst)
    } else {
      stop("ERROR: obj should be cell_data_set (Monocle3 object).")
    }
    
    ### set unique trace base
    if(is.null(levels(obj@colData[,trace_base]))) {
      unique_trace_base <- unique(as.character(obj@colData[,trace_base]))
    } else {
      unique_trace_base <- levels(obj@colData[,trace_base])
    }
    
    ### if a vertex is given, find the nearest trace_base of the vertex
    ### reduced_dim cords of the trace_bases (median center of each trace_base)
    trace_base_cords <- sapply(unique_trace_base, function(x) {
      target_cell_list <- rownames(obj@colData)[which(obj@colData[,trace_base] == x)]
      target_cell_cords <- obj@int_colData@listData$reducedDims[[reduced_dim_method[1]]][target_cell_list,]
      mean_x <- median(as.numeric(target_cell_cords[,1]))
      mean_y <- median(as.numeric(target_cell_cords[,2]))
      
      return(c(mean_x, mean_y))
    })
    
    ### find the closest base of the given vertex
    find_the_closest_base <- function(vtx) {
      vtx_cords <- obj@principal_graph_aux[[reduced_dim_method[1]]]$dp_mst[,vtx]
      dist_result <- sapply(colnames(trace_base_cords), function(x) {
        return(euclidean(vtx_cords, trace_base_cords[,x], FALSE))
      })
      closest_base <- names(dist_result)[which(dist_result == min(dist_result))]
      return(closest_base)
    }
    
    ### vertex - base mapping
    vertex_base_map <- data.frame(vertex=colnames(obj@principal_graph_aux[[reduced_dim_method[1]]]$dp_mst),
                                  base=NA,
                                  stringsAsFactors = FALSE, check.names = FALSE)
    rownames(vertex_base_map) <- vertex_base_map$vertex
    for(vtx in vertex_base_map$vertex) {
      vertex_base_map[vtx,"base"] <- find_the_closest_base(vtx)
    }
    
    ### re-construct the igraph based on the trace base
    new_adj_mat <- matrix(0, nrow = length(unique_trace_base), ncol = length(unique_trace_base))
    rownames(new_adj_mat) <- unique_trace_base
    colnames(new_adj_mat) <- unique_trace_base
    
    for(tr_bs in unique_trace_base) {
      ### what are the vertices that are affiliated to the given base
      target_vertex_list <- unique(vertex_base_map$vertex[which(vertex_base_map$base == tr_bs)])
      
      ### get vertices that are linked to the given vertices
      linked_vertex_list <- names(which(apply(original_adj_mat[target_vertex_list,], 2, sum) > 0))
      
      ### linked bases
      linked_base_list <- unique(vertex_base_map[linked_vertex_list,"base"])
      
      ### remove the self-link
      linked_base_list <- linked_base_list[which(linked_base_list != tr_bs)]
      
      ### fill into the new adj mat
      new_adj_mat[tr_bs,linked_base_list] <- 1
      new_adj_mat[linked_base_list, tr_bs] <- 1
    }
    
    ### diagonal = 0 (no self-linking)
    diag(new_adj_mat) <- 0
    
    ### remove nodes with no connections
    keep_idx <- which(apply(new_adj_mat, 1, sum) > 0)
    new_adj_mat <- new_adj_mat[keep_idx,keep_idx]
    
    ### build igraph with the new adj mat
    grp <- graph_from_adjacency_matrix(new_adj_mat)
    
    #
    ### set params for visualization
    #
    ### set node color
    if(is.null(color_scheme)) {
      V(grp)$color <- "gold"
    } else {
      V(grp)$color <- color_scheme[names(V(grp))]  
    }
    ### set node shape
    V(grp)$shape <- node_shape[1]
    ### set node size
    if(node_shape == "circle") {
      V(grp)$size <- 10
    } else if(node_shape == "rectangle") {
      V(grp)$size <- 25
    }
    ### set node label size
    V(grp)$label.cex <- 3
    ### set node label color
    V(grp)$label.color <- "black"
    ### set node label font
    V(grp)$label.font <- 2
    ### set edge width
    E(grp)$width <- 3
    ### set edge arrow size
    E(grp)$arrow.size <- 0
    ### set edge color
    E(grp)$color <- "black"
    
    ### trajectory tree plot
    pdf(result_path,
        width = 25, height = 30)
    plot.igraph(grp,
                layout = layout.reingold.tilford(grp, root=root),
                asp = 1.5,
                margin = 0)
    title(paste0("Monocle3 Trajectory Based on ", trace_base),
          cex.main=3, col.main="black")
    dev.off()
    
    
  }
  
  
  ### get the same color scheme as the UMAP
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "seurat_clusters",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  group_label_size = 7)
  g <- ggplot_build(p)
  color_sch <- unique(g$data[[2]]$colour)
  names(color_sch) <- unique(g$plot$data$cell_color)
  color_sch <- color_sch[order(as.numeric(names(color_sch)))]
  
  ### trajectory tree plot based on seurat_clusters
  plot_trajectory_as_tree(obj = cds,
                          trace_base = "seurat_clusters",
                          root = "5",
                          reduced_dim_method = "UMAP",
                          color_scheme = color_sch,
                          node_shape = "circle",
                          result_path = paste0(outputDir, "Monocle3_Trajectory_Tree_Clusters.pdf"))
  
  
  ### get the same color scheme as the UMAP
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "new_anno",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  group_label_size = 7)
  g <- ggplot_build(p)
  color_sch <- unique(g$data[[2]]$colour)
  names(color_sch) <- unique(g$plot$data$cell_color)
  color_sch <- color_sch[order(names(color_sch))]
  
  ### trajectory tree plot based on the new anno
  plot_trajectory_as_tree(obj = cds,
                          trace_base = "new_anno",
                          root = "Mo-2",
                          reduced_dim_method = "UMAP",
                          color_scheme = color_sch,
                          node_shape = "rectangle",
                          result_path = paste0(outputDir, "Monocle3_Trajectory_Tree_New_Anno.pdf"))
  
  
  ### remove "-" cluster and do the analysis again
  seurat_obj <- subset(seurat_obj,
                       cells = rownames(seurat_obj@meta.data)[which(seurat_obj$new_anno != "-")])
  
  ### construct a monocle cds
  cds <- new_cell_data_set(as(seurat_obj@assays$RNA@counts, 'sparseMatrix'),
                           cell_metadata = seurat_obj@meta.data,
                           gene_metadata = data.frame(gene_short_name = row.names(seurat_obj@assays$RNA@counts),
                                                      row.names = row.names(seurat_obj@assays$RNA@counts),
                                                      stringsAsFactors = FALSE, check.names = FALSE))
  
  ### overwrite the current CDS UMAP with original seurat UMAP
  cds@int_colData@listData$reducedDims$UMAP <- seurat_obj@reductions$umap@cell.embeddings
  
  ### pre-process the data
  cds <- preprocess_cds(cds, num_dim = 50)
  
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
  ggsave(paste0(outputDir, "Monocle3_trajectory_clusters_ORIGINAL_UMAP2.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### plot trajectory
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "new_anno",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  group_label_size = 7)
  ggsave(paste0(outputDir, "Monocle3_trajectory_annotation_ORIGINAL_UMAP2.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  ### get the same color scheme as the UMAP
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "seurat_clusters",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  group_label_size = 7)
  g <- ggplot_build(p)
  color_sch <- unique(g$data[[2]]$colour)
  names(color_sch) <- unique(g$plot$data$cell_color)
  color_sch <- color_sch[order(as.numeric(names(color_sch)))]
  
  ### trajectory tree plot based on seurat_clusters
  plot_trajectory_as_tree(obj = cds,
                          trace_base = "seurat_clusters",
                          root = "5",
                          reduced_dim_method = "UMAP",
                          color_scheme = color_sch,
                          node_shape = "circle",
                          result_path = paste0(outputDir, "Monocle3_Trajectory_Tree_Clusters2.pdf"))
  
  
  ### get the same color scheme as the UMAP
  p <- plot_cells(cds,
                  label_roots = FALSE,
                  color_cells_by = "new_anno",
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  group_label_size = 7)
  g <- ggplot_build(p)
  color_sch <- unique(g$data[[2]]$colour)
  names(color_sch) <- unique(g$plot$data$cell_color)
  color_sch <- color_sch[order(names(color_sch))]
  
  ### trajectory tree plot based on the new anno
  plot_trajectory_as_tree(obj = cds,
                          trace_base = "new_anno",
                          root = "Mo-2",
                          reduced_dim_method = "UMAP",
                          color_scheme = color_sch,
                          node_shape = "rectangle",
                          result_path = paste0(outputDir, "Monocle3_Trajectory_Tree_New_Anno2.pdf"))
  
  
  
  
  ### a function that finds genes that change over specific trajectory
  ### and generates a feature plot on UMAP and a heatmap across the trajectory
  plot_trajectory_gexp <- function(obj = cds,
                                   trace_base,
                                   specifix_bases = NULL,
                                   root,
                                   reduced_dim_method,
                                   result_path = "./") {
    
  }
  
  ### or just use slingshot function
  
  ### get slingshot object
  slingshot_obj_mnn <- slingshot(mnn_map,
                                 clusterLabels = mnn_map_time,
                                 start.clus = "GMP",
                                 end.clus = "3mo")
  slingshot_obj_mnn <- as.SlingshotDataSet(slingshot_obj_mnn)
  
  ### get colors for the clustering result
  cell_colors_clust <- cell_pal(unique(mnn_map_time), hue_pal())
  
  ### Trajectory inference
  png(paste0(outputDir2, "CARpos_Trajectory_Inference_Time_MNN.png"), width = 5000, height = 3000, res = 350)
  par(mar=c(7, 7, 7, 1), mgp=c(4,1,0))
  plot(reducedDim(slingshot_obj_mnn),
       main=paste("CAR+ Trajectory Inference"),
       col = cell_colors_clust[as.character(mnn_map_time)],
       pch = 19, cex = 1, cex.lab = 3, cex.main = 3, cex.axis = 2)
  lines(slingshot_obj_mnn, lwd = 4, type = "lineages", col = "black",
        show.constraints = TRUE, constraints.col = cell_colors_clust)
  legend("bottomleft", legend = names(cell_colors_clust), col = cell_colors_clust,
         pch = 19, cex = 1.5)
  dev.off()
  
  ### find genes that change their expression over the course of development
  ###  If "consecutive", then consecutive points along each lineage will be used as contrasts
  sce <- fitGAM(counts = mnn_map_exp, sds = slingshot_obj_mnn)
  ATres <- associationTest(sce, contrastType = "consecutive")
  
  ### give FDR
  ATres <- ATres[order(ATres$pvalue),]
  ATres$FDR <- p.adjust(p = ATres$pvalue, method = "BH")
  
  ### save the result
  write.xlsx2(data.frame(Gene=rownames(ATres),
                         ATres,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir2, "/CARpos_Trajectory_Inference_Pseudotime_DEGs_Slingshot.xlsx"),
              sheetName = "CARpos_Trajectory_Inference_Pseudotime_DEGs_Slingshot", row.names = FALSE)
  
  ### get those genes from the Slingshot
  topgenes <- rownames(ATres[order(ATres$FDR), ])[1:100]
  pst.ord <- order(sce$slingshot$pseudotime.Lineage1, na.last = NA)
  heatdata <- assays(sce)$counts[topgenes, pst.ord]
  heatclus <- JCC_Seurat_Obj@meta.data[colnames(heatdata),"time2"]
  
  ### draw the heatmap
  png(paste0(outputDir2, "CARpos_Trajectory_Inference_Pseudotime_DEGs_Slingshot_Heatmap.png"),
      width = 3000, height = 3000, res = 350)
  par(oma=c(0,3,0,3))
  heatmap.2(log1p(heatdata), col = colorpanel(24, low = "blue", high = "red"),
            scale = "none", dendrogram = "row", trace = "none",
            cexRow = 0.5, key.title = "", main = "Top 100 Genes Associated With The Pseudotime",
            Colv = FALSE, labCol = FALSE,  key.xlab = "log(Count+1)", key.ylab = "Frequency",
            ColSideColors = cell_colors_clust[heatclus])
  legend("left", inset = -0.1,
         box.lty = 0, cex = 0.8,
         title = "Time", xpd = TRUE,
         legend=names(cell_colors_clust),
         col=cell_colors_clust,
         pch=15)
  dev.off()
  
  
}
