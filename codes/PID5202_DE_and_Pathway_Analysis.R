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
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(org.Mm.eg.db, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("org.Mm.eg.db")
    require(org.Mm.eg.db, quietly = TRUE)
  }
  if(!require(msigdbr, quietly = TRUE)) {
    install.packages("msigdbr")
    library(msigdbr, quietly = TRUE)
  }
  if(!require(scales, quietly = TRUE)) {
    install.packages("scales")
    require(scales, quietly = TRUE)
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
  seurat_obj@meta.data[colnames(seurat_obj@assays$RNA@data)[which(seurat_obj@assays$RNA@data["Msr1",] > 0)],"Msr1"] <- "Msr1_Pos"
  
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
  
  ### annotate more
  seurat_obj$new_group1 <- paste(seurat_obj$treatment, seurat_obj$time, sep = ".")
  seurat_obj$new_group2 <- paste(seurat_obj$treatment, seurat_obj$time, seurat_obj$Msr1, sep = ".")
  
  ### UMAP plot with some annotations
  umap_labels <- c("time", "treatment", "Msr1", "new_group1", "new_group2")
  for(label_col in umap_labels) {
    p <- DimPlot(object = seurat_obj, reduction = "umap",
                 group.by = label_col,
                 pt.size = 1) +
      ggtitle("") +
      labs(color="") +
      scale_color_viridis_d() +
      theme_classic(base_size = 40) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48, color = "black", face = "bold"),
            axis.text = element_text(size = 48, color = "black", face = "bold"),
            axis.title = element_text(size = 48, color = "black", face = "bold"),
            legend.title = element_text(size = 30, color = "black", face = "bold"),
            legend.text = element_text(size = 24, color = "black", face = "bold")) +
      guides(colour = guide_legend(override.aes = list(size=10)))
    p[[1]]$layers[[1]]$aes_params$alpha <- 1
    ggsave(paste0(outputDir, "UMAP_PID5202_", label_col, ".pdf"), plot = p, width = 20, height = 15, dpi = 350)
    
    if(length(unique(seurat_obj@meta.data[,label_col])) > 8) {
      ncol <- 3
      label_size <- 12
    } else {
      ncol <- 2
      label_size <- 20
    }
    p <- DimPlot(object = seurat_obj, reduction = "umap",
                 group.by = label_col,
                 split.by = label_col,
                 pt.size = 1,
                 ncol = ncol) +
      ggtitle("") +
      labs(color="") +
      scale_color_viridis_d() +
      theme_classic(base_size = label_size) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 48, color = "black", face = "bold"),
            axis.text = element_text(size = 30, color = "black", face = "bold"),
            axis.title = element_text(size = 40, color = "black", face = "bold"),
            legend.title = element_text(size = 30, color = "black", face = "bold"),
            legend.text = element_text(size = 24, color = "black", face = "bold")) +
      guides(colour = guide_legend(override.aes = list(size=10)))
    p[[1]]$layers[[1]]$aes_params$alpha <- 1
    ggsave(paste0(outputDir, "UMAP_PID5202_", label_col, "_SPLIT.pdf"), plot = p, width = 20, height = 15, dpi = 350)
  }
  
  ### information table
  ### number of cells
  ### avg expressed gene # per cell
  ### top 5 marker genes
  ### based on new_group2
  
  ### find all markers
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
  seurat_obj <- SetIdent(object = seurat_obj,
                         cells = rownames(seurat_obj@meta.data),
                         value = seurat_obj$new_group2)
  all_markers <- FindAllMarkers(seurat_obj,
                                min.pct = 0.2,
                                logfc.threshold = 0.5,
                                test.use = "wilcox")
  
  ### write out the DE result
  write.xlsx2(data.frame(Gene=rownames(all_markers),
                         all_markers,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "/Markers_PID5202.xlsx"),
              sheetName = "All_Markers", row.names = FALSE)
  
  ### make an info table
  info_table <- data.frame(Group = levels(seurat_obj$new_group2),
                           Cell_Num = sapply(levels(seurat_obj$new_group2), function(x) {
                             return(length(which(seurat_obj$new_group2 == x)))
                           }),
                           Avg_Gene_Num = sapply(levels(seurat_obj$new_group2), function(x) {
                             target_mat <- seurat_obj@assays$RNA@counts[,which(seurat_obj$new_group2 == x)]
                             expressed_gene_num <- apply(target_mat, 2, function(y) {
                               return(length(which(y > 0)))
                             })
                             return(mean(expressed_gene_num))
                           }),
                           Markers = sapply(levels(seurat_obj$new_group2), function(x) {
                             top_genes <- paste(all_markers$gene[which(all_markers$cluster == x)][1:5], collapse = ";")
                             return(top_genes)
                           }),
                           stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out the info table
  write.xlsx2(info_table, file = paste0(outputDir, "/Info_Table_PID5202.xlsx"),
              sheetName = "Info", row.names = FALSE)
  
  
  ### compare the changes in frequencies of MSR1+ (representing the macrophages) and MSR- peritoneal macrophages
  ### in response to PBS/LPS, Dexamethasone/LPS, MSR1-ncADC/LPS and Isotype Ctrl-ncADC/LPS treatments
  ### bar plot with Msr1+ cell #
  ### multiple box/violin plots with p-values (Msr1 expression)
  
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
  
  #'****************************************************************************************
  #' Gene Set Enrichment Analysis function
  #' 
  #' It receives gene list (character vector) and signature profiles (named numeric vector)
  #' as inputs, performs GSEA and returns a table of GSEA result table and draws
  #' a GSEA plot. It is basically a statistical significance test to check how the
  #' given given genes are biased on the signature profiles.
  #' 
  #' Whether there are multiple gene sets or multiple signatures,
  #' multiple testing (FDR computation) is performed.
  #' But if the input gene set and the input signature are both lists with multiple
  #' items (The length of the two are both more than 1) then we return an error message.
  #' 
  #' The plot file names will be determined by names(gene_list) or names(signature)
  #' If length(gene_list) > 1, then names(gene_list) will be used and
  #' if length(signature) > 1, then names(signature) will be used as file names.
  #' If there is no list names, then file names will be "GSEA_Plot_i.png".
  #' Here, i indicates that the plot is from the i-th row of the GSEA result table.
  #' 
  #' * Some plot drawing codes were from Rtoolbox/R/ReplotGSEA.R written by Thomas Kuilman. 
  #'****************************************************************************************
  #' @title	run_gsea
  #' 
  #' @param gene_list   A list of character vectors containing gene names to be tested
  #' @param signature   A list of named numeric vectors of signature values for GSEA. The gene_list
  #'                    should be included in the names(signature)
  #' @param printPlot   If TRUE, it also generates GSEA plot of the results
  #'                    (Default = FALSE)
  #' @param fdr_cutoff  When printing GSEA plots, print them with the FDR < fdr_cutoff only
  #'                    (Default = 0.05)
  #' @param heatmap_color_type  when 'relative', the heatmap of the GSEA colors the bottom half of the
  #'                            absolute range of the signature as blue and the upper half as red
  #'                            when 'absolute', the heatmap of GSEA colors the negative signature as blue
  #'                            and the positives as red
  #' @param printPath   When printing GSEA plots, print them in the designated path
  #'                    (Default = "./")
  #' @param width       The width of the plot file
  #'                    (Default = 2000)
  #' @param height      The height of the plot file
  #'                    (Default = 1200)
  #' @param res         The resolution of the plot file
  #'                    (Default = 130)
  #' 
  #' @return 	          It tests bias of the "gene_list" on the "signature" range and
  #'                    returns a table including p-values and FDRs (adjusted p-values)
  #'                    If fdr_cutoff == TRUE, it also generates a GSEA plot with the result
  #' 
  run_gsea <- function(gene_list,
                       signature,
                       printPlot = FALSE,
                       fdr_cutoff = 0.05,
                       heatmap_color_type = c("relative", "absolute"),
                       width = 2000,
                       height = 1200,
                       res = 350,
                       printPath = "./") {
    
    ### load required libraries
    if(!require("fgsea", quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("fgsea")
      require("fgsea", quietly = TRUE)
    }
    if(!require("limma", quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("limma")
      require("limma", quietly = TRUE)
    } 
    if(!require("checkmate", quietly = TRUE)) {
      install.packages("checkmate")
      require("checkmate", quietly = TRUE)
    }
    
    ### argument checking
    assertList(gene_list)
    assertList(signature)
    assertLogical(printPlot)
    assertNumeric(fdr_cutoff)
    assertIntegerish(width)
    assertIntegerish(height)
    assertIntegerish(res)
    assertString(printPath)
    if(length(gene_list) > 1 && length(signature) > 1) {
      stop("ERROR: \"gene_list\" and \"signature\" cannot be both \"list\"")
    }
    
    ### set random seed
    set.seed(1234)
    
    ### run GSEA
    ### if there are more than one signatures
    if(length(signature) > 1) {
      ### combine GSEA results of every signature inputs
      for(i in 1:length(signature)) {
        temp <- data.frame(fgseaMultilevel(pathways = gene_list, stats = signature[[i]]))
        if(i == 1) {
          gsea_result <- temp
        } else {
          gsea_result <- rbind(gsea_result, temp)
        }
      }
      
      ### compute FDRs
      corrected_gsea_result <- gsea_result[order(gsea_result$pval),]
      corrected_gsea_result$padj <- p.adjust(corrected_gsea_result$pval, method = "BH")
      gsea_result <- corrected_gsea_result[rownames(gsea_result),]
    } else {
      ### if there are more than one gene sets
      gsea_result <- data.frame(fgseaMultilevel(pathways = gene_list, stats = signature[[1]], minSize = -Inf, maxSize = Inf))
    }
    
    ### order GSEA by FDR
    gsea_result <- gsea_result[order(gsea_result$padj, gsea_result$pval),]
    
    ### print GSEA plot
    sIdx <- which(gsea_result$padj < fdr_cutoff)
    if(length(sIdx) > 100) {
      sIdx <- sIdx[1:100]
    }
    if(printPlot && length(sIdx) > 0) {
      for(i in sIdx) {
        ### get required values ready
        if(length(signature) > 1) {
          geneset <- gene_list[[i]]
          stats <- signature[[i]]
          stats <- stats[order(-stats)]
          fileName <- names(signature)[i]
        } else {
          geneset <- gene_list[[gsea_result$pathway[i]]]
          stats <- signature[[1]]
          stats <- stats[order(-stats)]
          fileName <- gsea_result$pathway[i]
        }
        if(is.null(fileName)) {
          fileName <- paste0("GSEA_Plot_", i)
        }
        stats <- stats[!is.na(stats)]
        gsea.hit.indices <- which(names(stats) %in% geneset)
        es.temp <- calcGseaStat(stats, gsea.hit.indices, returnAllExtremes = TRUE)
        if(es.temp$res >= 0) {
          gsea.es.profile <- es.temp$tops
        } else {
          gsea.es.profile <- es.temp$bottoms
        }
        enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
        metric.range <- c(min(stats), max(stats))
        gsea.p.value <- round(gsea_result$pval[i] ,5)
        gsea.fdr <- round(gsea_result$padj[i] ,5)
        gsea.enrichment.score <- round(gsea_result$ES[i], 5)
        gsea.normalized.enrichment.score <- round(gsea_result$NES[i], 5)
        
        ### print GSEA result plot
        png(paste0(printPath, fileName, ".png"), width = width, height = height, res = res)
        
        ### set layout
        layout.show(layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2)))
        
        ### draw the GSEA plot
        par(mar = c(0, 5, 2, 2))
        plot(c(1, gsea.hit.indices, length(stats)),
             c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
             xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
             ylim = enrichment.score.range,
             main = list(fileName, font = 1, cex = 1),
             panel.first = {
               abline(h = seq(round(enrichment.score.range[1], digits = 1),
                              enrichment.score.range[2], 0.1),
                      col = "gray95", lty = 2)
               abline(h = 0, col = "gray50", lty = 2)
             }
        )
        
        ### add informative text to the GSEA plot
        plot.coordinates <- par("usr")
        if(es.temp$res < 0) {
          text(length(stats) * 0.01, plot.coordinates[3] * 0.98,
               paste("P-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
                     gsea.enrichment.score, "\nNormalized ES:",
                     gsea.normalized.enrichment.score), adj = c(0, 0))
        } else {
          text(length(stats) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
               paste("P-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
                     gsea.enrichment.score, "\nNormalized ES:",
                     gsea.normalized.enrichment.score), adj = c(1, 1))
        }
        
        ### draw hit indices
        par(mar = c(0, 5, 0, 2))
        plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
             ylab = "", xlim = c(1, length(stats)))
        abline(v = gsea.hit.indices, lwd = 0.75)
        
        ### create color palette for the heatmap
        par(mar = c(0, 5, 0, 2))
        if(heatmap_color_type[1] == "relative") {
          rank.colors <- stats - metric.range[1]
          rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
          rank.colors <- ceiling(rank.colors * 511 + 1)
          rank.colors <- colorRampPalette(c("blue", "white", "red"))(512)[rank.colors]
        } else {
          rank.colors1 <- stats[which(stats >= 0)]
          rank.colors1 <- rank.colors1 - min(rank.colors1)
          rank.colors1 <- rank.colors1 / (max(rank.colors1) - min(rank.colors1))
          rank.colors1 <- ceiling(rank.colors1 * 255 + 1)
          rank.colors1 <- colorRampPalette(c("white", "red"))(256)[rank.colors1]
          rank.colors2 <- stats[which(stats < 0)]
          rank.colors2 <- rank.colors2 - min(rank.colors2)
          rank.colors2 <- rank.colors2 / (max(rank.colors2) - min(rank.colors2))
          rank.colors2 <- ceiling(rank.colors2 * 255 + 1)
          rank.colors2 <- colorRampPalette(c("blue", "white"))(256)[rank.colors2]
          rank.colors <- c(rank.colors1, rank.colors2)
        }
        
        ### draw the heatmap
        rank.colors <- rle(rank.colors)
        barplot(matrix(rank.colors$lengths), col = rank.colors$values,
                border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(stats)))
        box()
        text(length(stats) / 2, 0.7,
             labels = "Signature")
        text(length(stats) * 0.01, 0.7, "Largest", adj = c(0, NA))
        text(length(stats) * 0.99, 0.7, "Smallest", adj = c(1, NA))
        
        ### draw signature values
        par(mar = c(5, 5, 0, 2))
        rank.metric <- rle(round(stats, digits = 2))
        plot(stats, type = "n", xaxs = "i",
             xlab = "Rank in ordered gene list", xlim = c(0, length(stats)),
             ylim = metric.range, yaxs = "i",
             ylab = "Signature values",
             panel.first = abline(h = seq(metric.range[1] / 2,
                                          metric.range[2] - metric.range[1] / 4,
                                          metric.range[2] / 2), col = "gray95", lty = 2))
        
        barplot(rank.metric$values, col = "lightgrey", lwd = 0.1,
                xlim = c(0, length(stats)), ylim = c(-1, 1),
                width = rank.metric$lengths, border = NA,
                space = 0, add = TRUE, xaxt = "n")
        box()
        
        ### print out the file
        dev.off()
      }
    }
    
    return(gsea_result)
    
  }
  
  ###
  ### define and compare the transcriptional signature
  ### in MSR+ and MSR- peritoneal cells with Dexamethasone/LPS treatment (in comparison to PBS/LPS);
  ### with MSR1-ncADC/LPS and Isotype Ctrl-ncADC/LPS treatments
  ### 1. [MSR1_ncADC.24HR_PRIOR_TO_LPS vs Ctrl_ncADC.24HR_PRIOR_TO_LPS] in Mrs1+ cells
  ### 2. [MSR1_ncADC.24HR_PRIOR_TO_LPS vs Ctrl_ncADC.24HR_PRIOR_TO_LPS] in Mrs1- cells
  ### 3. [Steroid.2HR_PRIOR_TO_LPS vs PBS.2HR_PRIOR_TO_LPS] in Mrs1+ cells
  ### 4. [Steroid.2HR_PRIOR_TO_LPS vs PBS.2HR_PRIOR_TO_LPS] in Mrs1- cells
  ###
  ### compare the transcriptional signature in MSR1+ and MSR- peritoneal cells
  ### with MSR1-ncADC/LPS and Isotype Ctrl-ncADC/LPS treatments, if any.
  ### 5. [Mrs1+ vs Mrs1-] in MSR1_ncADC.24HR_PRIOR_TO_LPS cells
  ### 6. [Mrs1+ vs Mrs1-] in Ctrl_ncADC.24HR_PRIOR_TO_LPS cells
  ### 7. [Mrs1+ vs Mrs1-] in Steroid.2HR_PRIOR_TO_LPS cells
  ### 8. [Mrs1+ vs Mrs1-] in Steroid.24HR_PRIOR_TO_LPS cells
  ### 9. [Mrs1+ vs Mrs1-] in PBS.2HR_PRIOR_TO_LPS cells
  ### 10. [Mrs1+ vs Mrs1-] in PBS.NO_LPS cells
  ###
  
  ### comparison set-up
  comparisons <- vector("list", length = 10)
  comparisons[[1]] <- c("MSR1_ncADC.24HR_PRIOR_TO_LPS.Msr1_Pos",
                        "Ctrl_ncADC.24HR_PRIOR_TO_LPS.Msr1_Pos")
  comparisons[[2]] <- c("MSR1_ncADC.24HR_PRIOR_TO_LPS.Msr1_Neg",
                        "Ctrl_ncADC.24HR_PRIOR_TO_LPS.Msr1_Neg")
  comparisons[[3]] <- c("Steroid.2HR_PRIOR_TO_LPS.Msr1_Pos",
                        "PBS.2HR_PRIOR_TO_LPS.Msr1_Pos")
  comparisons[[4]] <- c("Steroid.2HR_PRIOR_TO_LPS.Msr1_Neg",
                        "PBS.2HR_PRIOR_TO_LPS.Msr1_Neg")
  comparisons[[5]] <- c("MSR1_ncADC.24HR_PRIOR_TO_LPS.Msr1_Pos",
                        "MSR1_ncADC.24HR_PRIOR_TO_LPS.Msr1_Neg")
  comparisons[[6]] <- c("Ctrl_ncADC.24HR_PRIOR_TO_LPS.Msr1_Pos",
                        "Ctrl_ncADC.24HR_PRIOR_TO_LPS.Msr1_Neg")
  comparisons[[7]] <- c("Steroid.2HR_PRIOR_TO_LPS.Msr1_Pos",
                        "Steroid.2HR_PRIOR_TO_LPS.Msr1_Neg")
  comparisons[[8]] <- c("Steroid.24HR_PRIOR_TO_LPS.Msr1_Pos",
                        "Steroid.24HR_PRIOR_TO_LPS.Msr1_Neg")
  comparisons[[9]] <- c("PBS.2HR_PRIOR_TO_LPS.Msr1_Pos",
                        "PBS.2HR_PRIOR_TO_LPS.Msr1_Neg")
  comparisons[[10]] <- c("PBS.NO_LPS.Msr1_Pos",
                        "PBS.NO_LPS.Msr1_Neg")
  
  ### set empy list of the results
  de_results <- vector("list", length = 10)
  pathway_results <- vector("list", length = 10)
  gsea_results <- vector("list", length = 10)
  
  ### GSEA db preparation
  ### MSIGDB
  m_df <- msigdbr(species = "Mus musculus")
  m_list <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  rm(m_df)
  gc()
  
  ### set idents
  seurat_obj <- SetIdent(object = seurat_obj,
                         cells = rownames(seurat_obj@meta.data),
                         value = seurat_obj$new_group2)
  
  ### automated run for every comparison
  for(i in 1:length(comparisons)) {
    
    ### make the output directory for each comparison
    outputDir2 <- paste0(outputDir, "/", comparisons[[i]][1], "_vs_", comparisons[[i]][2], "/")
    dir.create(paste0(outputDir2, "GSEA"), recursive = TRUE)
    
    ### DE analysis
    de_results[[i]] <- FindMarkers(seurat_obj,
                                   ident.1 = comparisons[[i]][1],
                                   ident.2 = comparisons[[i]][2],
                                   min.pct = 0.1,
                                   logfc.threshold = 0.1,
                                   test.use = "wilcox")
    
    ### write out the DE result
    write.xlsx2(data.frame(Gene=rownames(de_results[[i]]),
                           de_results[[i]],
                           stringsAsFactors = FALSE, check.names = FALSE),
                file = paste0(outputDir2, "/DE_Results_",
                              comparisons[[i]][1], "_vs_",
                              comparisons[[i]][2], ".xlsx"),
                sheetName = paste0(comparisons[[i]][1], "_vs_",
                                   comparisons[[i]][2]), row.names = FALSE)
    
    ### dotplot to check the DE genes (top up-regulated and down-regulated)
    temp_seurat_obj <- subset(seurat_obj,
                              idents = c(comparisons[[i]][1], comparisons[[i]][2]))
    temp_seurat_obj$new_group2 <- factor(as.character(temp_seurat_obj$new_group2),
                                         levels = c(comparisons[[i]][1], comparisons[[i]][2]))
    top_and_bottom_genes <- rownames(de_results[[i]])[c(which(de_results[[i]]$avg_log2FC < 0)[1:10],
                                                        which(de_results[[i]]$avg_log2FC > 0)[1:10])]
    
    ### print out the dotplot
    p <- DotPlot(temp_seurat_obj,
                 features = rev(top_and_bottom_genes),
                 group.by = "new_group2") +
      scale_size(range = c(0, 20)) +
      xlab("") +
      ylab("") +
      coord_flip() +
      guides(color = guide_colorbar(title = "Scaled Expression")) +
      scale_color_gradientn(colours = c("lightgrey", "#a50026")) +
      theme_classic(base_size = 35) +
      theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
            axis.title = element_text(size = 35, hjust = 0.5, color = "black", face = "bold"),
            axis.text.x = element_text(angle = -60, size = 30, vjust = 0.5, hjust = 0, color = "black", face = "bold"),
            axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
            axis.text = element_text(color = "black", face = "bold"),
            legend.title = element_text(size = 30, color = "black", face = "bold"),
            legend.text = element_text(size = 25, color = "black", face = "bold"),
            legend.key.size = unit(0.7, 'cm'))
    ggsave(file = paste0(outputDir2, "/Dotplot_",
                         comparisons[[i]][1], "_vs_",
                         comparisons[[i]][2], ".pdf"),
           plot = p, width = 20, height = 30, dpi = 350)
    
    ### get entrez ids for the genes
    ### up-regulated in treatment only
    de_entrez_ids <- mapIds(org.Mm.eg.db,
                            rownames(de_results[[i]])[intersect(which(de_results[[i]]$p_val_adj < 0.01),
                                                                which(de_results[[i]]$avg_log2FC > 0))],
                            "ENTREZID", "SYMBOL")
    de_entrez_ids <- de_entrez_ids[!is.na(de_entrez_ids)]
    
    ### pathway analysis
    pathway_results[[i]] <- pathwayAnalysis_CP(geneList = de_entrez_ids,
                                               org = "mouse",
                                               database = "GO",
                                               title = paste0("Pathway_",
                                                              comparisons[[i]][1], "_vs_",
                                                              comparisons[[i]][2]),
                                               pv_threshold = 0.05,
                                               displayNum = 30,
                                               imgPrint = TRUE,
                                               dir = paste0(outputDir2, "/"))
    write.xlsx(pathway_results[[i]], file = paste0(outputDir2, "/Pathway_Table_",
                                                   comparisons[[i]][1], "_vs_",
                                                   comparisons[[i]][2], ".xlsx"),
               row.names = FALSE, sheetName = paste0("GO_Results"))
    
    ### GSEA
    signat <- de_results[[i]]$avg_log2FC
    names(signat) <- rownames(de_results[[i]])
    gsea_results[[i]] <- run_gsea(gene_list = m_list, signature = list(signat),
                                  fdr_cutoff = 0.05,
                                  printPlot = TRUE, printPath = paste0(outputDir2, "/GSEA/"))
    gsea_results[[i]] <- gsea_results[[i]][order(gsea_results[[i]]$padj),]
    
    ### write out the result file
    write.xlsx2(gsea_results[[i]], file = paste0(outputDir2, "/GSEA_Table_",
                                                comparisons[[i]][1], "_vs_",
                                                comparisons[[i]][2], ".xlsx"),
                sheetName = "GSEA_Results", row.names = FALSE)
    
    ### garbage collection
    gc()
    
  }
  
  
  
  
  
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
