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
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
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
      unique_trace_base <- intersect(levels(obj@colData[,trace_base]), unique(obj@colData[,trace_base]))
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
      closest_base <- names(dist_result)[which(dist_result == min(dist_result, na.rm = TRUE))]
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
  
  
  # ******************************************************************************************
  # Pathway Analysis with clusterProfiler package
  # Input: geneList     = a vector of gene Entrez IDs for pathway analysis [numeric or character]
  #        org          = organism that will be used in the analysis ["human" or "mouse"]
  #                       should be either "human" or "mouse"
  #        database     = pathway analysis database (KEGG or GO) ["KEGG" or "GO"]
  #        title        = title of the pathway figure [character]
  #        pv_threshold = pathway analysis p-value threshold (not DE analysis threshold) [numeric]
  #        displayNum   = the number of pathways that will be displayed [numeric]
  #                       (If there are many significant pathways show the few top pathways)
  #        imgPrint     = print a plot of pathway analysis [TRUE/FALSE]
  #        dir          = file directory path of the output pathway figure [character]
  #
  # Output: Pathway analysis results in figure - using KEGG and GO pathways
  #         The x-axis represents the number of DE genes in the pathway
  #         The y-axis represents pathway names
  #         The color of a bar indicates adjusted p-value from the pathway analysis
  #         For Pathview Result, all colored genes are found DE genes in the pathway,
  #         and the color indicates log2(fold change) of the DE gene from DE analysis
  # ******************************************************************************************
  pathwayAnalysis_CP <- function(geneList,
                                 org,
                                 database,
                                 title="Pathway_Results",
                                 pv_threshold=0.05,
                                 displayNum=Inf,
                                 imgPrint=TRUE,
                                 dir="./") {
    
    ### load library
    if(!require(clusterProfiler, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("clusterProfiler")
      require(clusterProfiler, quietly = TRUE)
    }
    if(!require(ggplot2)) {
      install.packages("ggplot2")
      library(ggplot2)
    }
    
    
    ### collect gene list (Entrez IDs)
    geneList <- geneList[which(!is.na(geneList))]
    
    if(!is.null(geneList)) {
      ### make an empty list
      p <- list()
      
      if(database == "KEGG") {
        ### KEGG Pathway
        kegg_enrich <- enrichKEGG(gene = geneList, organism = org, pvalueCutoff = pv_threshold)
        
        if(is.null(kegg_enrich)) {
          writeLines("KEGG Result does not exist")
          return(NULL)
        } else {
          kegg_enrich@result <- kegg_enrich@result[which(kegg_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(kegg_enrich@result) <= displayNum)) {
              result <- kegg_enrich@result
              description <- kegg_enrich@result$Description
            } else {
              result <- kegg_enrich@result[1:displayNum,]
              description <- kegg_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(kegg_enrich) > 0) {
              p[[1]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 30) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("KEGG ", title)) +
                theme(axis.text = element_text(size = 30))
              
              ggsave(file = paste0(dir, "KEGG_", title, ".png"), plot = p[[1]], width = 35, height = 15, dpi = 350)
            } else {
              writeLines("KEGG Result does not exist")
            }
          }
          
          return(kegg_enrich@result)
        }
      } else if(database == "GO") {
        ### GO Pathway
        if(org == "human") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Hs.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else if(org == "mouse") {
          go_enrich <- enrichGO(gene = geneList, OrgDb = 'org.Mm.eg.db', readable = T, ont = "BP", pvalueCutoff = pv_threshold)
        } else {
          go_enrich <- NULL
          writeLines(paste("Unknown org variable:", org))
        }
        
        if(is.null(go_enrich)) {
          writeLines("GO Result does not exist")
          return(NULL)
        } else {
          go_enrich@result <- go_enrich@result[which(go_enrich@result$p.adjust < pv_threshold),]
          
          if(imgPrint == TRUE) {
            if((displayNum == Inf) || (nrow(go_enrich@result) <= displayNum)) {
              result <- go_enrich@result
              description <- go_enrich@result$Description
            } else {
              result <- go_enrich@result[1:displayNum,]
              description <- go_enrich@result$Description[1:displayNum]
            }
            
            if(nrow(go_enrich) > 0) {
              p[[2]] <- ggplot(result, aes(x=Description, y=Count)) + labs(x="", y="Gene Counts") + 
                theme_classic(base_size = 30) + geom_bar(aes(fill = p.adjust), stat="identity") + coord_flip() +
                scale_x_discrete(limits = rev(description)) +
                guides(fill = guide_colorbar(ticks=FALSE, title="P.Val", barheight=10)) +
                ggtitle(paste0("GO ", title)) +
                theme(axis.text = element_text(size = 30))
              
              ggsave(file = paste0(dir, "GO_", title, ".png"), plot = p[[2]], width = 35, height = 15, dpi = 350)
            } else {
              writeLines("GO Result does not exist")
            }
          }
          
          return(go_enrich@result)
        }
      } else {
        stop("database prameter should be \"GO\" or \"KEGG\"")
      }
    } else {
      writeLines("geneList = NULL")
    }
  }
  
  
  ### a function that finds genes that change over specific trajectory
  ### and generates a feature plot on UMAP and a dotplot across the trajectory
  ### use simple wilcoxn test between many comparisons to find those genes
  ### and perform pathway analysis on the genes as well
  ###
  ### obj: seurat object
  ### trace_base: a column in the meta.data that the user is interested in
  ### base_order: a lineage order that the user wants to dig in
  ### fdr_threshold: differential expression threshold
  ### top_gene_num: the number of top genes that will be used in the dotplot for visualization
  ### organism: human or mouse
  ### result_path: the path of the result; must include file title as well (no extensions)
  ###              e.g., /usr/hkim29/file_title
  ###
  plot_trajectory_gexp <- function(obj,
                                   trace_base,
                                   base_order,
                                   fdr_threshold = 0.05,
                                   top_gene_num = 20,
                                   organism = c("human", "mouse"),
                                   result_path = "./") {
    
    ### library
    if(!require(poolr, quietly = TRUE)) {
      install.packages("poolr")
      require(poolr, quietly = TRUE)
    }
    if(!require(xlsx, quietly = TRUE)) {
      install.packages("xlsx")
      require(xlsx, quietly = TRUE)
    }
    
    ### check obj
    if(as.character(class(obj)) != "Seurat") {
      stop("ERROR: obj should be a seurat object.")
    }
    
    ### set unique trace base
    if(is.null(levels(obj@meta.data[,trace_base]))) {
      unique_trace_base <- unique(as.character(obj@meta.data[,trace_base]))
    } else {
      unique_trace_base <- intersect(levels(obj@meta.data[,trace_base]), unique(obj@meta.data[,trace_base]))
    }
    
    ### check base_order
    if(sum(as.numeric(base_order %in% unique_trace_base)) != length(base_order)) {
      stop("ERROR: each base_order should be one of trace_base.")
    }
    
    ### set ident with the trace_base
    obj <- SetIdent(object = obj,
                    cells = rownames(obj@meta.data),
                    value = obj@meta.data[,trace_base])
    
    ### find genes of across the trajectory
    de_genes <- list()
    shared_genes <- NULL
    cnt <- 1
    for(i in 1:(length(base_order)-1)) {
      de_genes[[cnt]] <- FindMarkers(obj,
                                     ident.1 = base_order[i],
                                     ident.2 = base_order[i+1],
                                     min.pct = 0.1,
                                     logfc.threshold = 0,
                                     test.use = "wilcox",
                                     verbose = TRUE)
      
      if(is.null(shared_genes)) {
        shared_genes <- rownames(de_genes[[i]])
      } else {
        shared_genes <- intersect(shared_genes, rownames(de_genes[[i]]))
      }
      
      names(de_genes)[cnt] <- paste(base_order[i], base_order[i+1], sep = "_vs_")
      cnt <- cnt + 1
    }
    
    ### find significantly increasing gexp and decreasing gexp
    ### and if any one de was significant among all the comparisons across the trajectory
    
    ### logFC data frame
    logFC_df <- sapply(de_genes, function(x) {
      return(x[shared_genes, "avg_log2FC"])
    }, USE.NAMES = TRUE)
    rownames(logFC_df) <- shared_genes
    logFC_df <- data.frame(logFC_df,
                           stringsAsFactors = FALSE, check.names = FALSE)
    
    ### sum of the signs
    sign_df <- apply(logFC_df, 1, function(x) {
      return(sum(sign(x)))
    })
    
    ### only sum = 3 or -3 
    sign_df <- sign_df[union(which(sign_df == 3),
                             which(sign_df == -3))]
    consistent_genes <- names(sign_df)
    
    ### fdr data frame - only for the consistent genes across the comparisons
    fdr_df <- sapply(de_genes, function(x) {
      return(x[consistent_genes, "p_val_adj"])
    }, USE.NAMES = TRUE)
    rownames(fdr_df) <- consistent_genes
    fdr_df <- data.frame(fdr_df,
                         stringsAsFactors = FALSE, check.names = FALSE)
    
    ### select genes that show < fdr_threshold for every comparison
    filtered_df <- apply(fdr_df, 1, function(x) {
      return(sum(as.numeric(x < fdr_threshold)))
    })
    filtered_df <- filtered_df[which(filtered_df == 3)]
    target_genes <- names(filtered_df)
    
    ### make a result data frame
    result_df <- data.frame(Gene=target_genes,
                            logFC_df[target_genes,],
                            fdr_df[target_genes,],
                            stringsAsFactors = FALSE, check.names = FALSE)
    colnames(result_df)[2:(length(de_genes)+1)] <- paste0("log2FC_", colnames(result_df)[2:(length(de_genes)+1)])
    colnames(result_df)[(length(de_genes)+2):ncol(result_df)] <- paste0("FDR_", colnames(result_df)[2:(length(de_genes)+1)])
    
    ### put Stouffer FDR
    result_df$Stouffer_FDR <- apply(result_df[,(length(de_genes)+2):ncol(result_df)], 1, function(x) {
      return(stouffer(x)$p)
    })
    
    ### put total sign
    result_df$Sign <- sign(sign_df)[result_df$Gene]
    result_df$Sign[which(result_df$Sign == 1)] <- "+"
    result_df$Sign[which(result_df$Sign == -1)] <- "-"
    
    ### order based on Stouffer FDR
    result_df <- result_df[order(result_df$Stouffer_FDR,
                                 result_df$Sign),]
    
    ### write out as an Excel file
    write.xlsx(result_df,
               file = paste0(result_path, "_Excel_Table.xlsx"),
               sheetName = "DE_Table",
               row.names = FALSE)
    
    ### subset only for the given base
    temp_obj <- subset(obj,
                       idents = base_order)
    
    ### set idents and the level
    temp_obj$new_anno <- factor(temp_obj$new_anno,
                                levels = base_order)
    temp_obj <- SetIdent(object = temp_obj,
                         cells = rownames(temp_obj@meta.data),
                         value = temp_obj@meta.data[,trace_base])
    
    ### print out a dotplot
    p <- DotPlot(temp_obj,
                 features = result_df$Gene[1:top_gene_num],
                 group.by = trace_base) +
      scale_size(range = c(0, 20)) +
      xlab("") +
      ylab("") +
      coord_flip() +
      guides(color = guide_colorbar(title = "Scaled Expression")) +
      scale_color_gradientn(colours = c("lightgrey", "#a50026")) +
      theme_classic(base_size = 35) +
      theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
            axis.title = element_text(size = 35, hjust = 0.5, color = "black", face = "bold"),
            axis.text.x = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
            axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
            axis.text = element_text(color = "black", face = "bold"),
            legend.title = element_text(size = 30, color = "black", face = "bold"),
            legend.text = element_text(size = 25, color = "black", face = "bold"),
            legend.key.size = unit(0.7, 'cm'))
    ggsave(file = paste0(result_path, "_Dotplot.pdf"),
           plot = p, width = 20, height = 15, dpi = 350)
    
    ### get entrez ids for the genes
    if(organism[1] == "human") {
      if(!require(org.Hs.eg.db, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install("org.Hs.eg.db")
        require(org.Hs.eg.db, quietly = TRUE)
      }
      de_entrez_ids <- mapIds(org.Hs.eg.db,
                              result_df$Gene[which(result_df$Stouffer_FDR < fdr_threshold)],
                              "ENTREZID", "SYMBOL")
      de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
    } else if(organism[1] == "mouse") {
      if(!require(org.Mm.eg.db, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install("org.Mm.eg.db")
        require(org.Mm.eg.db, quietly = TRUE)
      }
      de_entrez_ids <- mapIds(org.Mm.eg.db,
                              result_df$Gene[which(result_df$Stouffer_FDR < fdr_threshold)],
                              "ENTREZID", "SYMBOL")
      de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
    } else {
      stop("ERRROR: orginism should be either human or mouse.")
    }
    
    ### pathway analysis
    pathway_result <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                         org = organism[1],
                                         database = "GO",
                                         title = paste0(basename(result_path), "_Pathway_Results"),
                                         pv_threshold = 0.05,
                                         displayNum = 30,
                                         imgPrint = TRUE,
                                         dir = paste0(dirname(result_path), "/"))
    
    write.xlsx(pathway_result, file = paste0(result_path, "_Pathway_Table.xlsx"),
               row.names = FALSE, sheetName = paste0("GO_Results"))
    
    ### garbage collection
    gc()
    
  }
  
  ### Mo-4 -> Mo-5 -> cDC2-1 -> cDC2-5
  plot_trajectory_gexp(obj = seurat_obj,
                       trace_base = "new_anno",
                       base_order = c("Mo-4", "Mo-5", "cDC2-1", "cDC2-5"),
                       fdr_threshold = 0.01,
                       top_gene_num = 20,
                       organism = "mouse",
                       result_path = paste0(outputDir, "/Mo-4_To_cDC2-5"))
  
  
  ### Mo-6 vs cDC2-1
  ### DE amd pathway analysis
  
  ### set idents
  seurat_obj <- SetIdent(object = seurat_obj,
                         cells = rownames(seurat_obj@meta.data),
                         value = seurat_obj$new_anno)
  
  ### Mo-6 vs cDC2-1
  de_result <- FindMarkers(seurat_obj,
                           ident.1 = "Mo-6",
                           ident.2 = "cDC2-1",
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(de_result),
                         de_result,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/PID4604_MO-6_vs_cDC2-1.xlsx"),
              sheetName = "DE_MO-6_vs_cDC2-1", row.names = FALSE)
  
  ### Convert symbols to Entrez IDs
  de_entrez_ids <- mapIds(org.Mm.eg.db,
                          rownames(de_result)[intersect(which(de_result$p_val_adj < 1E-100),
                                                        which(abs(de_result$avg_log2FC) > 0.6))],
                          "ENTREZID", "SYMBOL")
  de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
  
  
  ### pathway analysis
  pathway_result <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                       org = "mouse",
                                       database = "GO",
                                       title = paste0("Pathway_MO-6_vs_cDC2-1"),
                                       pv_threshold = 0.05,
                                       displayNum = 30,
                                       imgPrint = TRUE,
                                       dir = paste0(outputDir))
  
  write.xlsx(pathway_result, file = paste0(outputDir, "Pathway_Table_MO-6_vs_cDC2-1.xlsx"),
             row.names = FALSE, sheetName = paste0("GO_Results"))
  
  
  ### compare two trajectories
  ### Mo-3 -> Mo-8 -> Mo-6 -> IM
  ### Mo-3 -> Mo-4 -> Mo-5 -> cDC2-1
  plot_trajectory_gexp(obj = seurat_obj,
                       trace_base = "new_anno",
                       base_order = c("Mo-3", "Mo-8", "Mo-6", "IM"),
                       fdr_threshold = 0.01,
                       top_gene_num = 20,
                       organism = "mouse",
                       result_path = paste0(outputDir, "/Mo-3_To_IM"))
  plot_trajectory_gexp(obj = seurat_obj,
                       trace_base = "new_anno",
                       base_order = c("Mo-3", "Mo-4", "Mo-5", "cDC2-1"),
                       fdr_threshold = 0.01,
                       top_gene_num = 20,
                       organism = "mouse",
                       result_path = paste0(outputDir, "/Mo-3_To_cDC2-1"))
  
  
  
  
  
  
}
