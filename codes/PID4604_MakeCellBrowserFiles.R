###
#   File name : MakeCellBrowserFiles.R
#   Author    : Hyunjin Kim
#   Date      : Jun 14, 2022
#   Email     : hyunjin.kim2@regeneron.com
#   Purpose   : Load Seurat object and make CellBrowser files.
#
#   Instruction
#               1. Source("MakeCellBrowserFiles.R")
#               2. Run the function "make_files" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_MakeCellBrowserFiles.R/MakeCellBrowserFiles.R")
#               > make_files(Seurat_RObj_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/tx_fil_dim20_res250.Robj",
#                            outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID4604/CellBrowser/")
###

make_files <- function(Seurat_RObj_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/tx_fil_dim20_res250.Robj",
                       outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID4604/CellBrowser/") {
  
  ### load libraries
  if(!require(Seurat, quietly = TRUE)) {
    install.packages("Seurat")
    require(Seurat, quietly = TRUE)
  }
  if(!require(remotes, quietly = TRUE)) {
    install.packages("remotes")
    require(remotes, quietly = TRUE)
  }
  if(!require(R.utils, quietly = TRUE)) {
    install.packages("R.utils")
    require(R.utils, quietly = TRUE)
  }
  if(!require(SeuratWrappers, quietly = TRUE)) {
    remotes::install_github('satijalab/seurat-wrappers')
    require(SeuratWrappers, quietly = TRUE)
  }
  ### install specific version of seurat if needed
  remotes::install_version(package = 'Seurat', version = package_version('3.1.5'))
  
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
  
  ### remove cells with new annotation as "-"
  seurat_obj <- subset(seurat_obj,
                       cells = rownames(seurat_obj@meta.data)[which(seurat_obj$new_anno != "-")])
  
  
  ### update the old seurat object to a new seurat object with the new version
  ### this has to be run because ExportToCellbrowser() needs the same seurat version
  ### between the installed seurat and the object seurat
  # seurat_obj <- UpdateSeuratObject(seurat_obj)
  # gc()
  
  
  ### make CellBrowser files
  ExportToCellbrowser(object = seurat_obj,
                      dir = outputDir,
                      markers.n = 100,
                      skip.markers = FALSE,
                      skip.expr.matrix = FALSE,
                      skip.metadata = FALSE,
                      skip.reductions = FALSE)
  
  ### Error in -partition$avg_logFC : invalid argument to unary operator
  ### while running embedded FindAllMarkers() in ExportToCellbrowser()
  ### so, I'm running it manually
  
  ### set ident with the seurat cluster info
  seurat_obj <- SetIdent(object = seurat_obj,
                         cells = rownames(seurat_obj@meta.data),
                         value = seurat_obj$seurat_clusters) 
  
  ### DE analysis
  de_result <- FindAllMarkers(seurat_obj,
                              min.pct = 0.2,
                              logfc.threshold = 0.3,
                              max.cells.per.ident = 10000,
                              test.use = "wilcox")
  
  ### top 30 DE genes only for each cluster
  remove_idx <- NULL
  for(clstr in as.character(unique(de_result$cluster))) {
    clstr_idx <- which(de_result$cluster == clstr)
    if(length(clstr_idx) > 30) {
      remove_idx <- c(remove_idx, clstr_idx[31:length(clstr_idx)])
    }
  }
  if(length(remove_idx) > 0) {
    de_result <- de_result[-remove_idx,]
  }
  
  ### save as markers.tsv
  ### there's something weird happen if I generate this file manually
  ### it's because of the header: just copy and paste other header to the new markers.tsv
  write.table(de_result,
              file = paste0(outputDir, "markers.tsv"),
              sep = "\t",
              col.names = NA,
              quote = FALSE)
  
  
  ### since we wants the new annotation as the basic idents, generate new ones based on new annotation; not based on seurat clusters
  
  ### set new_anno as Ident
  seurat_obj <- SetIdent(object = seurat_obj,
                         cells = rownames(seurat_obj@meta.data),
                         value = seurat_obj@meta.data$new_anno)
  
  ### create new output dir
  outputDir2 <- paste0(outputDir, "/new_anno/")
  dir.create(outputDir2, showWarnings = FALSE, recursive = TRUE)
  
  ### make CellBrowser files
  ExportToCellbrowser(object = seurat_obj,
                      dir = outputDir2,
                      markers.n = 100,
                      skip.markers = FALSE,
                      skip.expr.matrix = FALSE,
                      skip.metadata = FALSE,
                      skip.reductions = FALSE)
  
  ### there's something weird happened with markers.tsv,
  ### so run it manually
  de_result <- FindAllMarkers(seurat_obj,
                              min.pct = 0.2,
                              logfc.threshold = 0.3,
                              max.cells.per.ident = 10000,
                              test.use = "wilcox")
  
  ### top 30 DE genes only for each cluster
  remove_idx <- NULL
  for(clstr in as.character(unique(de_result$cluster))) {
    clstr_idx <- which(de_result$cluster == clstr)
    if(length(clstr_idx) > 30) {
      remove_idx <- c(remove_idx, clstr_idx[31:length(clstr_idx)])
    }
  }
  if(length(remove_idx) > 0) {
    de_result <- de_result[-remove_idx,]
  }
  
  ### save as markers.tsv
  ### there's something weird happen if I generate this file manually
  ### it's because of the header: just copy and paste other header to the new markers.tsv
  write.table(de_result,
              file = paste0(outputDir2, "markers.tsv"),
              sep = "\t",
              col.names = NA,
              quote = FALSE)
  
  ### since the CellBrowser didn't work with the output files
  ### I assume that the exp mat file is too large
  ### so I'm going to try with down-sampled version and see if it works
  
  ### down-sampling
  ### make sure every seurat_cluster has at least 500 cells
  ### if max. number is smaller than that, just keep all the cells
  # set.seed(1234)
  # fixed_cell_num <- 500
  # downsampling_idx <- NULL
  # for(clstr in levels(seurat_obj$seurat_clusters)) {
  #   target_idx <- which(seurat_obj$seurat_clusters == clstr)
  #   if(length(target_idx) > fixed_cell_num) {
  #     target_idx <- sample(target_idx, fixed_cell_num)
  #   }
  #   
  #   downsampling_idx <- c(downsampling_idx, target_idx)
  # }
  # 
  # ### subset based on the down-sampling
  # sub_seurat_obj <- subset(seurat_obj,
  #                          cells = rownames(seurat_obj@meta.data)[downsampling_idx])
  # 
  # ### new output dir
  # new_outputDir <- paste0(outputDir, "/downsampled/")
  # dir.create(new_outputDir, recursive = TRUE)
  # 
  # ### make CellBrowser files
  # ExportToCellbrowser(object = sub_seurat_obj,
  #                     dir = new_outputDir,
  #                     markers.n = 100,
  #                     skip.markers = FALSE,
  #                     skip.expr.matrix = FALSE,
  #                     skip.metadata = FALSE,
  #                     skip.reductions = FALSE)
  
  
  ### test it internally
  # cd /Users/hyunjin.kim2/Documents/SimpleTasks/results/PID4604/CellBrowser/
  # cbBuild -o ./ -p 8888
  ### if you see 'Point your internet browser to http://127.0.0.1:8888 (or the IP address of this server)'
  ### then go to http://127.0.0.1:8888
  
  ### recover the seurat version
  remotes::install_version(package = 'Seurat', version = package_version('4.1.1'))
  
}