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
  if(!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    require(ggplot2, quietly = TRUE)
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
  
  plot(density(seurat_obj@assays$RNA@data["Msr1",]),
       main = "Normalized Expression of Msr1",
       xlab = "Normalized Expression")
  plot(density(seurat_obj@assays$RNA@scale.data["Msr1",]),
       main = "Scaled Expression of Msr1",
       xlab = "Scaled Expression")
  
  seurat_obj$Msr1 <- "Msr1_Neg"
  seurat_obj@meta.data[colnames(seurat_obj@assays$RNA@scale.data)[which(seurat_obj@assays$RNA@scale.data["Msr1",] > 0)],"Msr1"] <- "Msr1_Pos"
  
  # > max(seurat_obj@assays$RNA@data["Msr1",which(seurat_obj$Msr1 == "Msr1_Neg")])
  # [1] 0.1812776
  
  ### test the umap plot
  DimPlot(seurat_obj,
          group.by = "meta_tissue",
          raster = TRUE)
  
  ### save the new seurat object
  saveRDS(seurat_obj,
          file = "/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_gex_new.rds")
  
  ### load the most updated seurat object
  seurat_obj <- readRDS("/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_gex_new.rds")
  
  ### check some things
  print(identical(rownames(seurat_obj@assays$RNA@counts), rownames(seurat_obj@assays$RNA@data)))
  print(identical(rownames(seurat_obj@assays$RNA@counts), rownames(seurat_obj@assays$RNA@scale.data)))
  print(identical(colnames(seurat_obj@assays$RNA@counts), colnames(seurat_obj@assays$RNA@data)))
  print(identical(colnames(seurat_obj@assays$RNA@counts), colnames(seurat_obj@assays$RNA@scale.data)))
  print(identical(rownames(seurat_obj@meta.data), colnames(seurat_obj@assays$RNA@counts)))
  
  ### compare the changes in frequencies of MSR1+ (representing the macrophages) and MSR- peritoneal macrophages
  ### in response to PBS/LPS, Dexamethasone/LPS, MSR1-ncADC/LPS and Isotype Ctrl-ncADC/LPS treatments
  ### bar plot with Msr1+ cell #
  ### multiple box/violin plots with p-values (Msr1 expression)
  
  ### annotate more
  seurat_obj$new_group1 <- paste(seurat_obj$treatment, seurat_obj$time, sep = ".")
  seurat_obj$new_group2 <- paste(seurat_obj$treatment, seurat_obj$time, seurat_obj$Msr1, sep = ".")
  
  ### frequency table
  plot_df <- data.frame(Group=unique(seurat_obj$new_group2),
                        Cell_Num=0,
                        stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 1:nrow(plot_df)) {
    plot_df$Cell_Num[i] <- length(which(seurat_obj$new_group2 == plot_df$Group[i]))
  }
  plot_df$Group1 <- sapply(plot_df$Group, function(x) {
    return(strsplit(x, split = ".", fixed = TRUE)[[1]][1])
  })
  plot_df$Group2 <- sapply(plot_df$Group, function(x) {
    return(strsplit(x, split = ".", fixed = TRUE)[[1]][2])
  })
  plot_df$Group3 <- sapply(plot_df$Group, function(x) {
    return(strsplit(x, split = ".", fixed = TRUE)[[1]][3])
  })
  plot_df$Group4 <- paste(plot_df$Group1, plot_df$Group2, sep = ".")
  
  ### set some levels
  plot_df$Group4 <- factor(plot_df$Group4, levels = rev(c("PBS.2HR_PRIOR_TO_LPS",
                                                          "Steroid.2HR_PRIOR_TO_LPS",
                                                          "Steroid.24HR_PRIOR_TO_LPS",
                                                          "MSR1_ncADC.24HR_PRIOR_TO_LPS",
                                                          "Ctrl_ncADC.24HR_PRIOR_TO_LPS",
                                                          "PBS.NO_LPS")))
  plot_df$Group3 <- factor(plot_df$Group3, levels = c("Msr1_Neg",
                                                      "Msr1_Pos"))
  
  ### bar plot
  p <- ggplot(data=plot_df, aes_string(x="Group4", y="Cell_Num", fill="Group3", label="Cell_Num")) +
    geom_bar(position = position_dodge(), stat = "identity") +
    ggtitle("") +
    xlab("") + ylab("The Number of Cells") +
    coord_flip() +
    geom_text(aes_string(x="Group4", y="Cell_Num", label = "Cell_Num"),
              position = position_dodge(1), hjust = 1.1,
              size = 10, color = "black") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_fill_manual(values = c("#f46d43", "#74add1"),
                      breaks = c("Msr1_Pos", "Msr1_Neg")) +
    theme_classic(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.title = element_text(size = 30, color = "black", face = "bold"),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 30, color = "black", face = "bold"),
          legend.key.size = unit(1, 'cm'),
          axis.ticks = element_blank())
  ggsave(file = paste0(outputDir, "Barplot_PID5202_Msr1_Frequencies.pdf"), plot = p,
         width = 20, height = 10, dpi = 330)
  
  ### set levels of new_group2
  seurat_obj$new_group2 <- factor(seurat_obj$new_group2,
                                  levels = c("PBS.2HR_PRIOR_TO_LPS.Msr1_Pos",
                                             "PBS.2HR_PRIOR_TO_LPS.Msr1_Neg",
                                             "Steroid.2HR_PRIOR_TO_LPS.Msr1_Pos",
                                             "Steroid.2HR_PRIOR_TO_LPS.Msr1_Neg",
                                             "Steroid.24HR_PRIOR_TO_LPS.Msr1_Pos",
                                             "Steroid.24HR_PRIOR_TO_LPS.Msr1_Neg",
                                             "MSR1_ncADC.24HR_PRIOR_TO_LPS.Msr1_Pos",
                                             "MSR1_ncADC.24HR_PRIOR_TO_LPS.Msr1_Neg",
                                             "Ctrl_ncADC.24HR_PRIOR_TO_LPS.Msr1_Pos",
                                             "Ctrl_ncADC.24HR_PRIOR_TO_LPS.Msr1_Neg",
                                             "PBS.NO_LPS.Msr1_Pos",
                                             "PBS.NO_LPS.Msr1_Neg"))
  
  ### color scale
  plot_colors <- colorRampPalette(c("#d73027", "#fee090", "#4575b4"))(length(levels(seurat_obj$new_group2)))
  names(plot_colors) <- levels(seurat_obj$new_group2)
  
  ### violin plot of Msr1 EXP among the groups
  seurat_obj <- SetIdent(object = seurat_obj,
                         cells = rownames(seurat_obj@meta.data),
                         value = seurat_obj$new_group2)
  p <- VlnPlot(seurat_obj, features = "rna_Msr1",
               raster = TRUE, pt.size = 0,
               cols = plot_colors)
  p[[1]] <- p[[1]] + geom_boxplot(width=0.1) +
    # stat_compare_means(size = 8) +
    xlab("") +
    ylab("Normalized Expression") +
    ggtitle("Msr1 Expression") +
    stat_summary(fun=mean, geom="point", size=3, color="black") +
    theme_classic(base_size = 40) +
    theme(legend.key.size = unit(3, 'cm'),
          legend.position = "none",
          legend.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          legend.text = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 35, color = "black", face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, vjust = 0.5, size = 25, color = "black", face = "bold"),
          axis.title = element_text(angle = 0, size = 30, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = -70, size = 15, vjust = 0.5, hjust = 0, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 25, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"))
  
  ### save the violin plot
  ggsave(file = paste0(outputDir, "Vlnplot_PID5202_Msr1_Frequencies.pdf"), plot = p, width = 20, height = 15, dpi = 350)
  
  
  ### define and compare the transcriptional signature
  ### in MSR+ and MSR- peritoneal cells with Dexamethasone/LPS treatment (in comparison to PBS/LPS);
  ###
  ### define and compare the transcriptional signature in MSR+ peritoneal cells with PBS/LPS, Dexamethasone/LPS,
  ### MSR1-ncADC/LPS and Isotype Ctrl-ncADC/LPS treatments; compare the transcriptional signature
  ### in MSR1+ and MSR- peritoneal cells with MSR1-ncADC/LPS and Isotype Ctrl-ncADC/LPS treatments, if any.
  
  ###
  
  
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
