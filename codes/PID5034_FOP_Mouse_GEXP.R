###
#   File name : PID5034_FOP_Mouse_GEXP.R
#   Author    : Hyunjin Kim
#   Date      : Jun 13, 2022
#   Email     : hyunjin.kim2@regeneron.com
#   Purpose   : See gene expression differences among tissues in mouse endothelial cells
#
#   Instruction
#               1. Source("PID5034_FOP_Mouse_GEXP.R")
#               2. Run the function "gexp_analysis" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_PID5034_FOP_Mouse_GEXP.R/PID5034_FOP_Mouse_GEXP.R")
#               > gexp_analysis(gexp_mat_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/PID5034/exprMatrix.tsv",
#                               meta_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/PID5034/meta.tsv",
#                               outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5034/")
###

gexp_analysis <- function(gexp_mat_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/PID5034/exprMatrix.tsv",
                          meta_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/PID5034/meta.tsv",
                          outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5034/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(ggpubr, quietly = TRUE)) {
    install.packages("ggpubr")
    require(ggpubr, quietly = TRUE)
  }
  
  ### set genes of interest
  interesting_genes <- c("Inhba", "Inhbb",
                         "Acvr1", "Acvr1b", "Acvr1c", "Acvrl1", "Bmpr1a", "Bmpr1b", "Tgfbr1",
                         "Acvr2a", "Acvr2b", "Amhr2", "Bmpr2", "Tgfbr2")
  
  ### load files
  gexp <- read.table(file = gexp_mat_path,
                     header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)
  gene_names <- strsplit(gexp$gene, split = "|", fixed = TRUE)
  gene_names <- sapply(gene_names, function(x) x[1])
  rownames(gexp) <- gene_names
  gexp <- gexp[,-1]
  
  meta <- read.table(file = meta_path,
                     header = TRUE, sep = "\t", row.names = 1,
                     stringsAsFactors = FALSE, check.names = FALSE)
  
  ### create seurat object
  seurat_obj <- CreateSeuratObject(gexp, project = "Mouse_Endo_Atlas", assay = "RNA",
                                   min.cells = 0, min.features = 0, names.field = 1,
                                   names.delim = "_", meta.data = meta)
  rm(gexp)
  rm(meta)
  gc()
  
  ### check whether the orders are the same
  print(identical(rownames(seurat_obj@meta.data), colnames(seurat_obj@assays$RNA@counts)))
  print(identical(names(Idents(object = seurat_obj)), rownames(seurat_obj@meta.data)))
  
  ### save the obj as RDS file
  # saveRDS(seurat_obj, file = "./data/PID5034/Mouse_Endo_Atlas.RDS")
  
  ### load the seurat object
  seurat_obj <- readRDS(file = "./data/PID5034/Mouse_Endo_Atlas.RDS")
  
  ### factorize the tissue column
  seurat_obj$Tissue <- factor(seurat_obj$Tissue,
                              levels = levels(as.factor(seurat_obj$Tissue)))
  
  ### set tissue as Ident
  seurat_obj <- SetIdent(object = seurat_obj,
                         cells = rownames(seurat_obj@meta.data),
                         value = seurat_obj@meta.data$Tissue)
  
  ### color setting
  tissue_colors <- colorRampPalette(colors=c("#f46d43", "#fdae61", "#fee090", "#abd9e9", "#74add1", "#4575b4"))(length(unique(seurat_obj$Tissue)))
  names(tissue_colors) <- levels(as.factor(seurat_obj$Tissue))
  
  ### violin plot
  p <- VlnPlot(seurat_obj, features = interesting_genes,
               pt.size = 0, ncol = 5)
  for(i in 1:length(interesting_genes)) {
    p[[i]] <- p[[i]] + geom_boxplot(width=0.1) +
      # stat_compare_means(size = 3) +
      stat_summary(fun=mean, geom="point", size=0.5, color="red") +
      theme_classic(base_size = 12) +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, vjust = 0.5, color = "black", face = "bold"),
            axis.title = element_text(angle = 0, vjust = 0.5, hjust = 0.5, color = "black", face = "bold"),
            axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0, color = "black", face = "bold"))
  }
  
  ### save the violin plot
  ggsave(file = paste0(outputDir, "ViolinPlot_FOP_Mouse_Endo_GEXP.pdf"), plot = p, width = 15, height = 10, dpi = 350)
  
  
  ### factorize the tissue column
  seurat_obj$Tissue <- factor(seurat_obj$Tissue,
                              levels = rev(levels(as.factor(seurat_obj$Tissue))))
  
  ### dot plot
  p <- DotPlot(seurat_obj,
               features = interesting_genes,
               group.by = "Tissue") +
    scale_size(range = c(5, 20)) +
    xlab("") +
    ylab("") +
    scale_color_gradientn(colours = c("#313695", "#ffffbf", "#a50026")) +
    theme_classic(base_size = 35) +
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
          axis.title = element_text(size = 35, hjust = 0.5, color = "black", face = "bold"),
          axis.text.x = element_text(angle = 90, size = 30, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.text.y = element_text(angle = 0, size = 35, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
          axis.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(size = 30, color = "black", face = "bold"),
          legend.text = element_text(size = 25, color = "black", face = "bold"),
          legend.key.size = unit(0.7, 'cm'))
  ggsave(file = paste0(outputDir, "DotPlot_FOP_Mouse_Endo_GEXP.pdf"),
         plot = p, width = 20, height = 15, dpi = 350)
  
  
}
