###
#   File name : MakeDesignCSV_For_Cooper.R
#   Author    : Hyunjin Kim
#   Date      : Jul 5, 2022
#   Email     : hyunjin.kim2@regeneron.com
#   Purpose   : Make a design csv file for Cooper's single cell preprocessing
#
#   Instruction
#               1. Source("MakeDesignCSV_For_Cooper.R")
#               2. Run the function "make_design_csv" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_MakeDesignCSV_For_Cooper.R/MakeDesignCSV_For_Cooper.R")
#               > make_design_csv(meta_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/PID5202/qcSummary_p5202e16639.xlsx",
#                                 outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/")
###

make_design_csv <- function(meta_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/PID5202/qcSummary_p5202e16639.xlsx",
                            outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/") {
  
  ### load libraries
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### load the data
  meta_data <- read.xlsx(file = meta_path,
                         sheetIndex = 1,
                         row.names = 1,
                         stringsAsFactors = FALSE, check.names = FALSE)
  meta_data <- data.frame(t(meta_data),
                          stringsAsFactors = FALSE, check.names = FALSE)
  
  ### prepare for the columns
  sample_ID <- rownames(meta_data)
  rna_data <- paste0("/mnt/MP", meta_data$`path to the results`, "/outs/filtered_feature_bc_matrix/")
  tcr_data <- ""
  bcr_data <- ""
  outdir <- "/data/home/hyunjin.kim2/projects/PID5202/data/"
  
  species <- meta_data$species
  cell_type <- meta_data$`Cell Type`
  tissue_detail <- meta_data$tissue_detail
  
  result_table <- data.frame(sample_ID = sample_ID,
                             rna_data = rna_data,
                             tcr_data = tcr_data,
                             bcr_data = bcr_data,
                             outdir = outdir,
                             species = species,
                             cell_type = cell_type,
                             tissue = tissue_detail,
                             stringsAsFactors = FALSE, check.names = FALSE)
  
  ### write out as a csv file
  write.csv(result_table,
            file = paste0(outputDir, "PID5202_prep_table.csv"),
            row.names = FALSE)
  
}
