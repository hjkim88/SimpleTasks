###
#   File name : PID5094_MakeCellBrowserFiles.R
#   Author    : Hyunjin Kim
#   Date      : May 30, 2022
#   Email     : hyunjin.kim2@regeneron.com
#   Purpose   : Load Seurat object and make CellBrowser files.
#
#   Instruction
#               1. Source("PID5094_MakeCellBrowserFiles.R")
#               2. Run the function "make_files" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_PID5094_MakeCellBrowserFiles.R/PID5094_MakeCellBrowserFiles.R")
#               > make_files(Seurat_RObj_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/tx_rmct_dim10_res03_label.Robj",
#                            outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5094/CellBrowser/")
###

make_files <- function(Seurat_RObj_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/tx_rmct_dim10_res03_label.Robj",
                       outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5094/CellBrowser/") {
  
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
  
  ### test it internally
  # cd /Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5094/CellBrowser/
  # cbBuild -o ./ -p 8888
  ### if you see 'Point your internet browser to http://127.0.0.1:8888 (or the IP address of this server)'
  ### then go to http://127.0.0.1:8888
  
  ### recover the seurat version
  remotes::install_version(package = 'Seurat', version = package_version('4.1.1'))
  
}