###
#   File name : PID5137_Lineage_Analysis.R
#   Author    : Hyunjin Kim
#   Date      : May 31, 2022
#   Email     : hyunjin.kim2@regeneron.com
#   Purpose   : Run monocle3 on the given single cell dataset.
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
  
  
  ### tree visualization function for Monocle3
  plot_complex_cell_trajectory <- function(cds,
                                           dim = "UMAP",
                                           x=1, 
                                           y=2, 
                                           root_states = NULL,
                                           color_by="State", 
                                           show_tree=TRUE, 
                                           show_backbone=TRUE, 
                                           backbone_color="black", 
                                           markers=NULL, 
                                           show_cell_names=FALSE, 
                                           cell_size=1.5,
                                           cell_link_size=0.75,
                                           cell_name_size=2,
                                           show_branch_points=TRUE, 
                                           ...){
    gene_short_name <- NA
    sample_name <- NA
    data_dim_1 <- NA
    data_dim_2 <- NA
    
    # need to validate cds as ready for this plot (need mst, pseudotime, etc)
    lib_info_with_pseudo <- pData(cds)
    
    if (is.null(cds@reduce_dim_aux)){
      stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
    }
    
    if(dim == "UMAP") {
      reduced_dim_coords <- cds@reduce_dim_aux$UMAP$model$umap_model$embedding
    } else if(dim == "PCA") {
      reduced_dim_coords <- cds@reduce_dim_aux$PCA$model$umap_model$embedding
    } else {
      stop("Error: dim should be either UMAP or PCA.")
    }
    
    if (is.null(reduced_dim_coords)){
      stop("You must first call reduceDimension() before using this function")
    }
    
    dp_mst <- minSpanningTree(cds)
    
    
    if(is.null(root_states)) {
      if(is.null(lib_info_with_pseudo$Pseudotime)){
        root_cell <- row.names(lib_info_with_pseudo)[degree(dp_mst) == 1][1]
      }
      else
        root_cell <- row.names(subset(lib_info_with_pseudo, Pseudotime == 0))
      
      if(cds@dim_reduce_type != "ICA")
        root_cell <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_cell, ]] 
      
    }
    else {
      candidate_root_cells <- row.names(subset(pData(cds), State %in% root_states))
      if(cds@dim_reduce_type == "ICA") {
        root_cell <- candidate_root_cells[which(degree(dp_mst, candidate_root_cells) == 1)]
      } else {
        Y_candidate_root_cells <- V(dp_mst)$name[cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[candidate_root_cells, ]] 
        root_cell <- Y_candidate_root_cells[which(degree(dp_mst, Y_candidate_root_cells) == 1)]
      }
      
    }
    
    tree_coords <- layout_as_tree(dp_mst, root=root_cell)
    
    #ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(x,y),]))
    ica_space_df <- data.frame(tree_coords)
    row.names(ica_space_df) <- colnames(reduced_dim_coords)
    colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")
    
    ica_space_df$sample_name <- row.names(ica_space_df)
    #ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
    #print(ica_space_with_state_df)
    
    
    if (is.null(dp_mst)){
      stop("You must first call orderCells() before using this function")
    }
    
    edge_list <- as.data.frame(get.edgelist(dp_mst))
    colnames(edge_list) <- c("source", "target")
    
    edge_df <- merge(ica_space_df, edge_list, by.x="sample_name", by.y="source", all=TRUE)
    #edge_df <- ica_space_df
    edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="source_prin_graph_dim_1", "prin_graph_dim_2"="source_prin_graph_dim_2"))
    edge_df <- merge(edge_df, ica_space_df[,c("sample_name", "prin_graph_dim_1", "prin_graph_dim_2")], by.x="target", by.y="sample_name", all=TRUE)
    edge_df <- plyr::rename(edge_df, c("prin_graph_dim_1"="target_prin_graph_dim_1", "prin_graph_dim_2"="target_prin_graph_dim_2"))
    
    #S_matrix <- reducedDimS(cds)
    #data_df <- data.frame(t(S_matrix[c(x,y),]))
    
    if(cds@dim_reduce_type == "ICA"){
      S_matrix <- tree_coords[,] #colnames(cds)
      
    } else if(cds@dim_reduce_type %in% c("DDRTree", "SimplePPT", "SGL-tree")){
      S_matrix <- tree_coords[closest_vertex,]
      closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    }
    
    data_df <- data.frame(S_matrix)
    row.names(data_df) <- colnames(reducedDimS(cds))
    colnames(data_df) <- c("data_dim_1", "data_dim_2")
    data_df$sample_name <- row.names(data_df)
    data_df <- merge(data_df, lib_info_with_pseudo, by.x="sample_name", by.y="row.names")
    
    markers_exprs <- NULL
    if (is.null(markers) == FALSE){
      markers_fData <- subset(fData(cds), gene_short_name %in% markers)
      if (nrow(markers_fData) >= 1){
        markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
        #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
        markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
        markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
      }
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
      data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
      #print (head(edge_df))
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, I(cell_size))) + facet_wrap(~feature_label)
    }else{
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) 
    }
    if (show_tree){
      g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=cell_link_size, linetype="solid", na.rm=TRUE, data=edge_df)
    }
    
    # FIXME: setting size here overrides the marker expression funtionality. 
    # Don't do it!
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
      if(class(data_df[, color_by]) == 'numeric') {
        g <- g + geom_jitter(aes_string(color = paste0("log10(", color_by, " + 0.1)")), size=I(cell_size), na.rm = TRUE, height=5) + 
          scale_color_viridis(name = paste0("log10(", color_by, ")"), ...)
      } else {
        g <- g + geom_jitter(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE, height=5) 
      }
    }else {
      if(class(data_df[, color_by]) == 'numeric') {
        g <- g + geom_jitter(aes_string(color = paste0("log10(", color_by, " + 0.1)")), size=I(cell_size), na.rm = TRUE, height=5) + 
          scale_color_viridis(name = paste0("log10(", color_by, " + 0.1)"), ...)
      } else {
        g <- g + geom_jitter(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE, height=5)
      }
    }
    
    if (show_branch_points && cds@dim_reduce_type == 'DDRTree'){
      mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
      branch_point_df <- subset(edge_df, sample_name %in% mst_branch_nodes)[,c("sample_name", "source_prin_graph_dim_1", "source_prin_graph_dim_2")]
      branch_point_df$branch_point_idx <- match(branch_point_df$sample_name, mst_branch_nodes)
      branch_point_df <- branch_point_df[!duplicated(branch_point_df$branch_point_idx), ]
      
      g <- g + geom_point(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2"), 
                          size=2 * cell_size, na.rm=TRUE, data=branch_point_df) +
        geom_text(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", label="branch_point_idx"), 
                  size=1.5 * cell_size, color="white", na.rm=TRUE, data=branch_point_df)
    }
    if (show_cell_names){
      g <- g +geom_text(aes(label=sample_name), size=cell_name_size)
    }
    g <- g + 
      #scale_color_brewer(palette="Set1") +
      theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
      theme(panel.border = element_blank()) +
      # theme(axis.line.x = element_line(size=0.25, color="black")) +
      # theme(axis.line.y = element_line(size=0.25, color="black")) +
      theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
      theme(panel.background = element_rect(fill='white')) +
      theme(legend.key=element_blank()) + 
      xlab('') + 
      ylab('') +
      theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
      #guides(color = guide_legend(label.position = "top")) +
      theme(legend.key = element_blank()) +
      theme(panel.background = element_rect(fill='white')) + 
      theme(line = element_blank(), 
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()) 
    g
  }
  
  ### tree representation of the overall trajectory
  ### root cluster = 5
  
  
  
  
  
  
  
  
  
}