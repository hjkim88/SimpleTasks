###
#   File name : DESeq2_GPP.R
#   Author    : Hyunjin Kim
#   Date      : Jul 5, 2022
#   Email     : hyunjin.kim2@regeneron.com
#   Purpose   : Run the same DESeq2 (as Wei Keat's example) pipeline to GPP (bulk) data
#
#   Instruction
#               1. Source("DESeq2_GPP.R")
#               2. Run the function "deseq2_analysis" - specify the input file path and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_DESeq2_GPP.R/DESeq2_GPP.R")
#               > deseq2_analysis(gexp_mat_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/GPP/counts_matrix_GPP_allGenes.txt",
#                                 meta_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/GPP/E-MTAB-11144.sdrf.txt",
#                                 outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/GPP/")
###

deseq2_analysis <- function(gexp_mat_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/GPP/counts_matrix_GPP_allGenes.txt",
                            meta_path="/Users/hyunjin.kim2/Documents/SimpleTasks/data/GPP/E-MTAB-11144.sdrf.txt",
                            outputDir="/Users/hyunjin.kim2/Documents/SimpleTasks/results/GPP/") {
  
  ### load libraries
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(biomaRt, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("biomaRt")
    require(biomaRt, quietly = TRUE)
  }
  if(!require(dplyr, quietly = TRUE)) {
    install.packages("dplyr")
    require(dplyr, quietly = TRUE)
  }
  
  ### load the data
  gexp_data <- read.table(file = gexp_mat_path,
                          header = TRUE,
                          sep = "\t",
                          row.names = 1,
                          stringsAsFactors = FALSE, check.names = FALSE)
  meta_data <- read.table(file = meta_path,
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE)
  
  ### refine some columns in the meta data
  meta_data$`Factor Value[organism part]`[which(meta_data$`Factor Value[organism part]` %in% c("", "  ", "        "))] <- NA
  meta_data$`Factor Value[sampling site]`[which(meta_data$`Factor Value[sampling site]` %in% c("", "  ", "        "))] <- NA
  meta_data$`Factor Value[compound]`[which(meta_data$`Factor Value[compound]` %in% c("", "  ", "        "))] <- NA
  meta_data$`Factor Value[time]`[which(meta_data$`Factor Value[time]` %in% c("", "  ", "        "))] <- NA
  meta_data$`Unit[time unit]`[which(meta_data$`Unit[time unit]` %in% c("", "  ", "        "))] <- NA
  
  ### make a comparison group column
  meta_data$comp_group <- paste(meta_data$`Factor Value[organism part]`,
                                meta_data$`Factor Value[sampling site]`,
                                meta_data$`Factor Value[compound]`,
                                meta_data$`Factor Value[time]`,
                                meta_data$`Unit[time unit]`,
                                sep = "_")
  meta_data$comp_group <- make.names(meta_data$comp_group)
  
  
  ### An IDS filtering function
  ### Only include genes with at least 75% of counts greater than 10 in any one group
  ###
  ### count_mat: count matrix - row: genes, column: samples
  ### group_info: a character vector of group info that has the same length as the ncol(count_mat)
  ### groups_in_use: a vector of two groups for DE analysis
  ###                must be one of the group_info element
  ### count_threshold: the raw count cut-off
  ### sample_percentage_threshold: the cut-off for the percentage of the samples
  ###
  ids_bulk_rnaseq_filtering <- function(count_mat,
                                        group_info,
                                        groups_in_use,
                                        count_threshold = 10,
                                        sample_percentage_threshold = 0.75) {
    
    ### check the params
    if(ncol(count_mat) != length(group_info)) {
      stop("ERROR: ncol(count_mat) must be the same as length(group_info).")
    }
    if(length(which(groups_in_use %in% unique(group_info))) != length(groups_in_use)) {
      stop("ERROR: groups_in_use must be one of the group_info element.")
    }
    if(length(groups_in_use) != 2) {
      stop("ERROR: length(groups_in_use) must be 2.")
    }
    
    ### get group indicies
    group_idx <- vector("list", length = length(groups_in_use))
    names(group_idx) <- groups_in_use
    
    for(grp in names(group_idx)) {
      group_idx[[grp]] <- which(group_info == grp)
    }
    
    ### filtering
    count_mat <- count_mat[which((((rowSums(count_mat[,group_idx[[groups_in_use[1]]]]>=count_threshold))/length(group_idx[[groups_in_use[1]]]))>=sample_percentage_threshold) | (((rowSums(count_mat[,group_idx[[groups_in_use[2]]]]>=count_threshold))/length(group_idx[[groups_in_use[2]]]))>=sample_percentage_threshold)),]
    
    return(count_mat)
    
  }
  
  
  #####################################################
  ### A function to perform DE analysis with DESeq2 ###
  #####################################################
  #' @title deseqWithComparisons
  #' @param rCnt raw count matrix (rows: genes, columns: samples)
  #' @param grp a character vector of class info of the samples
  #'            or a list with two character vectors representing
  #'            two different independent classes of the samples
  #' @param exp_class a string of the experiment group's name
  #'                  or a vector with two strings that are from
  #'                  the grp[[1]] for a comparison
  #' @param ctrl_class a string of the control group's name
  #'                  or a vector with two strings that are from
  #'                  the grp[[2]] for a comparison
  #'                  
  #' * If exp_class = "A" and ctrl_class = "B", the direction of
  #'   differential expression is A - B. If exp_class = c("A", "B") and
  #'   and ctrl_class = c("C", "D"), then the direction of
  #'   differential expression is (AC - AD) - (BC - BD).
  #'   
  #' @param bat_eff a character vector of batch effect info of the samples
  #' @param filtering filtering method for the raw count
  #' @param adj_p_thresh numeric. Filters out from the results genes with
  #' 		                 adjusted p-value larger than this value
  #' @return data.frame
  #' @export
  #' @author Hyunjin Kim
  ####################################################
  deseqWithComparisons <- function(rCnt,
                                   grp,
                                   exp_class,
                                   ctrl_class,
                                   bat_eff=NULL,
                                   filtering=c("default", "ids_method"),
                                   adj_p_thresh = 1) {
    
    ### load library
    if(!require(DESeq2, quietly = TRUE)) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("DESeq2")
      require(DESeq2, quietly = TRUE)
    }
    if(!require("checkmate", quietly = TRUE)) {
      install.packages("checkmate")
      require("checkmate", quietly = TRUE)
    }
    
    ### argument checking
    assertDataFrame(rCnt)
    assert(checkVector(grp), checkList(grp))
    if(testList(grp)) {
      assert(length(grp) == 2)
      assert(length(exp_class) == 2)
      assert(length(ctrl_class) == 2)
      assert(length(which(grp[[1]] == exp_class[1])) > 0)
      assert(length(which(grp[[1]] == exp_class[2])) > 0)
      assert(length(which(grp[[2]] == ctrl_class[1])) > 0)
      assert(length(which(grp[[2]] == ctrl_class[2])) > 0)
    } else {
      assert(length(exp_class) == 1)
      assert(length(ctrl_class) == 1)
      assert(length(which(grp == exp_class)) > 0, length(which(grp == ctrl_class)) > 0)
    }
    assertCharacter(exp_class)
    assertCharacter(ctrl_class)
    assert(checkNull(bat_eff), checkCharacter(bat_eff), checkFactor(bat_eff))
    assertNumeric(adj_p_thresh)
    
    ### sometimes, there are some variables which can not be transformed into R variable names
    ### so, just to be safe, change all the variables to R-usuable ones
    exp_class <- make.names(exp_class)
    ctrl_class <- make.names(ctrl_class)
    
    ### there are two different grp options
    if(testList(grp)) {
      sampleType1 <- make.names(as.character(grp[[1]]))
      sampleType2 <- make.names(as.character(grp[[2]]))
      
      ### there are two options: considering batch effect or not
      if(is.null(bat_eff)) {
        ### make a data frame for design matrix
        Coldata <- data.frame(sampleType1, sampleType2)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rCnt)
        Coldata$sampleType1 <- relevel(as.factor(Coldata$sampleType1), ref = exp_class[2])
        Coldata$sampleType2 <- relevel(as.factor(Coldata$sampleType2), ref = ctrl_class[2])
        
        ### data preparation for DE analysis
        deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType1*sampleType2)
      } else {
        ### make a data frame for design matrix
        batch_eff <- make.names(as.character(bat_eff))
        Coldata <- data.frame(sampleType1, sampleType2, batch_eff)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rCnt)
        Coldata$sampleType1 <- relevel(as.factor(Coldata$sampleType1), ref = exp_class[2])
        Coldata$sampleType2 <- relevel(as.factor(Coldata$sampleType2), ref = ctrl_class[2])
        Coldata$batch_eff <- as.factor(Coldata$batch_eff)
        
        ### data preparation for DE analysis
        deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType1*sampleType2+batch_eff)
      }
      
      ### filtering some useless genes out
      deSeqData <- deSeqData[rowSums(counts(deSeqData))>1,]
      
      ### run DE analysis
      dea <- DESeq(deSeqData)
      
      ### get DE genes using the contrast
      deresults <- results(dea, name = paste0("sampleType1", exp_class[1], ".sampleType2", ctrl_class[1]))
    } else {
      sampleType <- make.names(as.character(grp))
      
      ### there are two options: considering batch effect or not
      if(is.null(bat_eff)) {
        ### make a data frame for design matrix
        Coldata <- data.frame(sampleType)
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rCnt)
        Coldata$sampleType <- relevel(as.factor(Coldata$sampleType), ref = ctrl_class)
        
        ### filtering some useless genes out
        if(filtering[1] == "ids_method") {
          rCnt <- ids_bulk_rnaseq_filtering(count_mat = rCnt,
                                            group_info = sampleType,
                                            groups_in_use = c(exp_class, ctrl_class),
                                            count_threshold = 10,
                                            sample_percentage_threshold = 0.75)
        } else {
          rCnt <- rCnt[rowSums(counts(rCnt))>1,]
        }
        
        ### data preparation for DE analysis
        deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType)
      } else {
        ### make a data frame for design matrix
        batch_eff <- make.names(as.character(bat_eff))
        Coldata <- data.frame(sampleType, batch_eff)
        
        ### filtering some useless genes out
        if(filtering[1] == "ids_method") {
          rCnt <- ids_bulk_rnaseq_filtering(count_mat = rCnt,
                                            group_info = sampleType,
                                            groups_in_use = c(exp_class, ctrl_class),
                                            count_threshold = 10,
                                            sample_percentage_threshold = 0.75)
        } else {
          rCnt <- rCnt[rowSums(counts(rCnt))>1,]
        }
        
        ### set names and the refererence group
        rownames(Coldata) <- colnames(rCnt)
        Coldata$sampleType <- relevel(as.factor(Coldata$sampleType), ref = ctrl_class)
        Coldata$batch_eff <- as.factor(Coldata$batch_eff)
        
        ### data preparation for DE analysis
        deSeqData <- DESeqDataSetFromMatrix(countData=rCnt, colData=Coldata, design= ~sampleType+batch_eff)
      }
      
      ### run DE analysis
      dea <- DESeq(deSeqData)
      
      ### get DE genes using the contrast
      deresults <- results(dea, contrast = c("sampleType", exp_class, ctrl_class))
    }
    
    ### change p-values that have NA values to 1
    if(length(which(is.na(deresults$pvalue))) > 0) {
      deresults$pvalue[which(is.na(deresults$pvalue))] <- 1
    }
    if(length(which(is.na(deresults$padj))) > 0) {
      deresults$padj[which(is.na(deresults$padj))] <- 1
    }
    
    ### ordering and filtering the DE result
    deresults <- deresults[order(deresults$padj, deresults$pvalue, na.last = TRUE), ,drop = FALSE]
    deresults <- deresults[deresults$padj <= adj_p_thresh, ,drop = FALSE]
    
    ### there are two different grp options
    if(testList(grp)) {
      ### add baseMean for each group
      nCnt <- counts(dea, normalized=TRUE, replaced=TRUE)
      exp1_ctrl1_rowMeans <- apply(nCnt[,intersect(which(Coldata$sampleType1 == exp_class[1]),
                                                   which(Coldata$sampleType2 == ctrl_class[1])),
                                        drop=FALSE], 1, mean)
      exp1_ctrl2_rowMeans <- apply(nCnt[,intersect(which(Coldata$sampleType1 == exp_class[1]),
                                                   which(Coldata$sampleType2 == ctrl_class[2])),
                                        drop=FALSE], 1, mean)
      exp2_ctrl1_rowMeans <- apply(nCnt[,intersect(which(Coldata$sampleType1 == exp_class[2]),
                                                   which(Coldata$sampleType2 == ctrl_class[1])),
                                        drop=FALSE], 1, mean)
      exp2_ctrl2_rowMeans <- apply(nCnt[,intersect(which(Coldata$sampleType1 == exp_class[2]),
                                                   which(Coldata$sampleType2 == ctrl_class[2])),
                                        drop=FALSE], 1, mean)
      deresults <- data.frame(baseMean=deresults[,1],
                              V1=exp1_ctrl1_rowMeans[rownames(deresults)],
                              V2=exp1_ctrl2_rowMeans[rownames(deresults)],
                              V3=exp2_ctrl1_rowMeans[rownames(deresults)],
                              V4=exp2_ctrl2_rowMeans[rownames(deresults)],
                              deresults[,2:6],
                              stringsAsFactors = FALSE, check.names = FALSE)
      colnames(deresults)[2:5] <- c(paste0("baseMean_", exp_class[1], "_", ctrl_class[1]),
                                    paste0("baseMean_", exp_class[1], "_", ctrl_class[2]),
                                    paste0("baseMean_", exp_class[2], "_", ctrl_class[1]),
                                    paste0("baseMean_", exp_class[2], "_", ctrl_class[2]))
    } else {
      ### add baseMean for each group
      nCnt <- counts(dea, normalized=TRUE, replaced=TRUE)
      exp_rowMeans <- apply(nCnt[,which(Coldata$sampleType == exp_class), drop=FALSE], 1, mean)
      ctrl_rowMeans <- apply(nCnt[,which(Coldata$sampleType == ctrl_class), drop=FALSE], 1, mean)
      deresults <- data.frame(baseMean=deresults[,1],
                              V1=exp_rowMeans[rownames(deresults)],
                              V2=ctrl_rowMeans[rownames(deresults)],
                              deresults[,2:6],
                              stringsAsFactors = FALSE, check.names = FALSE)
      colnames(deresults)[2:3] <- c(paste0("baseMean_", exp_class), paste0("baseMean_", ctrl_class))  
    }
    
    return(deresults)
    
  }
  
  
  ### biomart objects
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  ensembl_genes <-  rownames(gexp_data)
  ensembl_symbol_map <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),
                              values = ensembl_genes, mart= mart)
  rownames(ensembl_symbol_map) <- ensembl_symbol_map$ensembl_gene_id
  
  ### automatically run DE analysis for the specific comparisons
  ###
  ### Skin: lesion/none vs non-lesional/none
  ### Skin: lesion/spesolimab vs lesional/none
  ### Blood: spesolimab 1 week vs none
  ### Blood: spesolimab 2 week vs none
  ### Blood: spesolimab 4 week vs none
  ###
  comparisons <- list()
  comparisons[[1]] <- c("skin_lesion_none_NA_NA", "skin_non.lesional.site_none_NA_NA")
  comparisons[[2]] <- c("skin_lesion_spesolimab_1_week", "skin_lesion_none_NA_NA")
  comparisons[[3]] <- c("blood_NA_spesolimab_1_week", "blood_NA_none_NA_NA")
  comparisons[[4]] <- c("blood_NA_spesolimab_2_week", "blood_NA_none_NA_NA")
  comparisons[[5]] <- c("blood_NA_spesolimab_4_week", "blood_NA_none_NA_NA")
  
  de_results <- vector("list", length(comparisons))
  names(de_results) <- sapply(comparisons, function(x) paste0(x[1], "_vs_", x[2]))
  
  for(comp in comparisons) {
    
    ### print comparison name
    comp_name <- paste0(comp[1], "_vs_", comp[2])
    print(comp_name)
    
    ## run the DE analysis for [Skin: lesion/none vs non-lesional/none]
    de_results[[comp_name]] <- deseqWithComparisons(rCnt = gexp_data,
                                                    grp = as.character(meta_data$comp_group),
                                                    exp_class = as.character(comp[1]),
                                                    ctrl_class = as.character(comp[2]),
                                                    bat_eff = as.character(meta_data$`Characteristics[individual]`),
                                                    filtering = "ids_method",
                                                    adj_p_thresh = 1)
    
    ### annotate gene symbols for the Ensembl IDs
    de_results[[comp_name]]$EID <- rownames(de_results[[comp_name]])
    de_results[[comp_name]] <- left_join(de_results[[comp_name]], ensembl_symbol_map, by = c("EID"="ensembl_gene_id"))
    de_results[[comp_name]]$hgnc_symbol[which(de_results[[comp_name]]$hgnc_symbol %in% c(NA, ""))] <- ""
    de_results[[comp_name]] <- data.frame(Ensembl_ID=de_results[[comp_name]]$EID,
                                          Gene_Symbol=de_results[[comp_name]]$hgnc_symbol,
                                          de_results[[comp_name]][,-which(colnames(de_results[[comp_name]]) %in% c("EID", "hgnc_symbol"))],
                                          stringsAsFactors = FALSE, check.names = FALSE)
    
    ### write out the result as an Excel file
    write.xlsx(de_results[[comp_name]],
               file = paste0(outputDir, "/", "DESeq2_", comp_name, ".xlsx"),
               sheetName = comp_name,
               row.names = FALSE)
    
  }
    
}
