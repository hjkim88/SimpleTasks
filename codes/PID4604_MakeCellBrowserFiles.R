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
  
  ### update the old seurat object to a new seurat object with the new version
  ### this has to be run because ExportToCellbrowser() needs the same seurat version
  ### between the installed seurat and the object seurat
  # seurat_obj <- UpdateSeuratObject(seurat_obj)
  # gc()
  
  
  ### make CellBrowser files
  ExportToCellbrowser(object = seurat_obj,
                      dir = outputDir)
  
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
                              max.cells.per.ident = 100,
                              test.use = "wilcox")
  
  ### top 100 DE genes only for each cluster
  remove_idx <- NULL
  for(clstr in as.character(unique(de_result$cluster))) {
    clstr_idx <- which(de_result$cluster == clstr)
    if(length(clstr_idx) > 100) {
      remove_idx <- c(remove_idx, clstr_idx[101:length(clstr_idx)])
    }
  }
  if(length(remove_idx) > 0) {
    de_result <- de_result[-remove_idx,]
  }
  
  ### save as markers.tsv
  write.table(de_result,
              file = paste0(outputDir, "markers.tsv"),
              sep = "\t",
              col.names = NA,
              quote = FALSE)
  
  ### test it internally
  # cd /Users/hyunjin.kim2/Documents/SimpleTasks/results/PID4604/CellBrowser/
  # cbBuild -o ./ -p 8888
  ### if you see 'Point your internet browser to http://127.0.0.1:8888 (or the IP address of this server)'
  ### then go to http://127.0.0.1:8888
  
  ### recover the seurat version
  remotes::install_version(package = 'Seurat', version = package_version('4.1.1'))
  
}