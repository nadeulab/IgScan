#' Recalculate IgScan IDs for Single Cell Data
#'
#' @description This function recalculates and standardizes IgScan clonotype-related IDs
#' in single-cell sequencing datasets. Using this function may be of interest
#' when a IgScan-annotated single cell object has been filtered and some cells
#' have been removed. It supports input objects from both the `SingleCellExperiment`
#' and `Seurat` classes.
#'
#' @param single_cell_object A single-cell object combined with IgScan annotation
#' in either `SingleCellExperiment` or `Seurat` format.
#'
#' @return A single cell object with the IgScan clonotype-related IDs updated in the metadata.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' recalculated_object <- recalculate_IDs_single_cell(sc_object)
#' }
#'
recalculate_IDs_single_cell <- function(single_cell_object){

  if(class(single_cell_object)[1] == "SingleCellExperiment"){
    data_frame <- single_cell_object@colData
  } else if(class(single_cell_object)[1] == "Seurat"){
    data_frame <- single_cell_object@meta.data
  }

  ## To correct ClonotypeID
  clt_dict <- aggregate(x = data_frame$ClonotypeID[!is.na(data_frame$ClonotypeID)], by = list(data_frame$ClonotypeID[!is.na(data_frame$ClonotypeID)]), FUN = length)
  clt_dict$chain <- sapply(clt_dict$Group.1, function(x) substr(x, 1, 3))
  clt_dict$cloneID_no_x <- sapply(clt_dict$Group.1, function(x) strsplit(x, "x")[[1]][1])
  clt_dict <- clt_dict[order(clt_dict$x, decreasing = T),]

  chain_count <- list(IGH = 1, IGK = 1, IGL = 1)
  corrected_clones <- c()
  for(row in 1:nrow(clt_dict)){
    if(clt_dict$cloneID_no_x[row] %in% corrected_clones){
      clt_dict$NewName[row] <- clt_dict$NewName[clt_dict$cloneID_no_x == clt_dict$cloneID_no_x[row]][1]
    } else{
      clt_dict$NewName[row] <- paste0(clt_dict$chain[row], ".C", chain_count[[clt_dict$chain[row]]])
      corrected_clones <- c(corrected_clones, clt_dict$cloneID_no_x[row])
      chain_count[[clt_dict$chain[row]]] <- chain_count[[clt_dict$chain[row]]]+1
    }
  }

  if(any(duplicated(clt_dict$NewName))){
    for(sim_c in clt_dict$NewName[duplicated(clt_dict$NewName)]){
      clt_dict$NewName[clt_dict$NewName == sim_c] <- paste0(sim_c, "x", 1:length(clt_dict$NewName[clt_dict$NewName == sim_c]))
    }
  }
  data_frame$ClonotypeID <- clt_dict$NewName[match(data_frame$ClonotypeID, clt_dict$Group.1)]

  ## To correct ClonotypeVariantID
  cv_dict_all <- data.frame()
  for(c in unique(data_frame$ClonotypeID)){
    cv_dict <- aggregate(x = data_frame$ClonotypeVariantID[data_frame$ClonotypeID == c], by = list(data_frame$ClonotypeVariantID[data_frame$ClonotypeID == c]), FUN = length)
    cv_dict <- cv_dict[order(cv_dict$x, decreasing = T),]
    cv_dict$NewName <- paste0(c, ".CV", 1:nrow(cv_dict))
    cv_dict_all <- rbind(cv_dict_all, cv_dict)
  }
  data_frame$ClonotypeVariantID <- cv_dict_all$NewName[match(data_frame$ClonotypeVariantID, cv_dict_all$Group.1)]

  ## To correct SubcloneID
  sbc_dict_all <- data.frame()
  for(cv in unique(data_frame$ClonotypeVariantID)){
    sbc_dict <- aggregate(x = data_frame$SubcloneID[data_frame$ClonotypeVariantID == cv], by = list(data_frame$SubcloneID[data_frame$ClonotypeVariantID == cv]), FUN = length)
    sbc_dict <- sbc_dict[order(sbc_dict$x, decreasing = T),]
    sbc_dict$NewName <- paste0(cv, ".S", 1:nrow(sbc_dict))
    sbc_dict_all <- rbind(sbc_dict_all, sbc_dict)
  }
  data_frame$SubcloneID <- sbc_dict_all$NewName[match(data_frame$SubcloneID, sbc_dict_all$Group.1)]

  ## Update merged IDs
  for(cell in unique(data_frame$igSubcloneID)){
    old_cltID_list <- strsplit(unique(data_frame$igClonotypeID[data_frame$igSubcloneID == cell]), "-")[[1]]
    data_frame$igClonotypeID[data_frame$igSubcloneID == cell] <- paste(clt_dict$NewName[match(old_cltID_list, clt_dict$Group.1)], collapse = "-")

    old_cvID_list <- strsplit(unique(data_frame$igClonotypeVariantID[data_frame$igSubcloneID == cell]), "-")[[1]]
    data_frame$igClonotypeVariantID[data_frame$igSubcloneID == cell] <- paste(cv_dict_all$NewName[match(old_cvID_list, cv_dict_all$Group.1)], collapse = "-")

    old_scID_list <- strsplit(cell, "-")[[1]]
    data_frame$igSubcloneID[data_frame$igSubcloneID == cell] <- paste(sbc_dict_all$NewName[match(old_scID_list, sbc_dict_all$Group.1)], collapse = "-")
  }

  ## Correct numeric IDs based on new merged IDs
  per_cell_clt_dict <- aggregate(x = data_frame$igClonotypeID_num[data_frame$igClonotypeID_num != ""], by = list(data_frame$igClonotypeID_num[data_frame$igClonotypeID_num != ""]), FUN = length)
  per_cell_clt_dict <- per_cell_clt_dict[order(per_cell_clt_dict$x, decreasing = T),]
  per_cell_clt_dict$NewName <- paste0("C", 1:nrow(per_cell_clt_dict))
  data_frame$igClonotypeID_num <- per_cell_clt_dict$NewName[match(data_frame$igClonotypeID_num, per_cell_clt_dict$Group.1)]

  per_cell_cv_dict_all <- data.frame()
  per_cell_sbc_in_clt_dict_all <- data.frame()
  for(c in unique(data_frame$igClonotypeID_num[!is.na(data_frame$igClonotypeID_num)])){
    cv_dict <- aggregate(x = data_frame$igClonotypeVariantID_num[data_frame$igClonotypeID_num == c], by = list(data_frame$igClonotypeVariantID_num[data_frame$igClonotypeID_num == c]), FUN = length)
    cv_dict <- cv_dict[order(cv_dict$x, decreasing = T),]
    cv_dict$NewName <- paste0(c, ".CV", 1:nrow(cv_dict))
    per_cell_cv_dict_all <- rbind(per_cell_cv_dict_all, cv_dict)

    sbc_in_clt_dict <- aggregate(x = data_frame$igSubcloneID_in_Clonotype_num[data_frame$igClonotypeID_num == c], by = list(data_frame$igSubcloneID_in_Clonotype_num[data_frame$igClonotypeID_num == c]), FUN = length)
    sbc_in_clt_dict <- sbc_in_clt_dict[order(sbc_in_clt_dict$x, decreasing = T),]
    sbc_in_clt_dict$NewName <- paste0(c, ".", 1:nrow(sbc_in_clt_dict))
    per_cell_sbc_in_clt_dict_all <- rbind(per_cell_sbc_in_clt_dict_all, sbc_in_clt_dict)
  }
  data_frame$igClonotypeVariantID_num <- per_cell_cv_dict_all$NewName[match(data_frame$igClonotypeVariantID_num, per_cell_cv_dict_all$Group.1)]
  data_frame$igSubcloneID_in_Clonotype_num <- per_cell_sbc_in_clt_dict_all$NewName[match(data_frame$igSubcloneID_in_Clonotype_num, per_cell_sbc_in_clt_dict_all$Group.1)]

  per_cell_sbc_in_cv_dict_all <- data.frame()
  for(cv in unique(data_frame$igClonotypeVariantID_num[!is.na(data_frame$igClonotypeVariantID_num)])){
    sbc_in_cv_dict <- aggregate(x = data_frame$igSubcloneID_in_ClonotypeVariant_num[data_frame$igClonotypeVariantID_num == cv], by = list(data_frame$igSubcloneID_in_ClonotypeVariant_num[data_frame$igClonotypeVariantID_num == cv]), FUN = length)
    sbc_in_cv_dict <- sbc_in_cv_dict[order(sbc_in_cv_dict$x, decreasing = T),]
    sbc_in_cv_dict$NewName <- paste0(cv, ".S", 1:nrow(sbc_in_cv_dict))
    per_cell_sbc_in_cv_dict_all <- rbind(per_cell_sbc_in_cv_dict_all, sbc_in_cv_dict)
  }
  data_frame$igSubcloneID_in_ClonotypeVariant_num <- per_cell_sbc_in_cv_dict_all$NewName[match(data_frame$igSubcloneID_in_ClonotypeVariant_num, per_cell_sbc_in_cv_dict_all$Group.1)]

  if(class(single_cell_object)[1] == "SingleCellExperiment"){
    single_cell_object@colData <- data_frame
  } else if(class(single_cell_object)[1] == "Seurat"){
    single_cell_object@meta.data <- data_frame
  }

  return(single_cell_object)
}
