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
#' @param group_col A vector with the name of the column or columns containing the grouping variable/s.
#' Default is 'orig.ident', thus recalculating BCR IDs based in each sample independently.
#' @param threads The number of threads to perform BCR ID recalculation. Default is 1.
#'
#' @return A single cell object with the IgScan clonotype-related IDs updated in the metadata.
#'
#' @export
#'
#' @import dplyr
#' @importFrom qs qread qsave
#' @import Seurat
#' @import SeuratObject
#' @import SingleCellExperiment
#' @importFrom parallel mclapply
#'
#' @examples
#' \dontrun{
#' recalculated_object <- recalculate_IDs_single_cell(sc_object, group_col = "orig.ident", threads = 4)
#' }
#'
recalculate_IDs_single_cell <- function(single_cell_object, group_col = "orig.ident", threads = 1){

  if(class(single_cell_object)[1] == "SingleCellExperiment"){
    meta_data <- as.data.frame(single_cell_object@colData)
  } else if(class(single_cell_object)[1] == "Seurat"){
    meta_data <- as.data.frame(single_cell_object@meta.data)
  }

  if(!all(group_col %in% colnames(meta_data))){stop(paste0("\nUnknown column (", group_col[!group_col %in% colnames(meta_data)], ") selected for BCR ID recalculation! Please, set a valid column name."))}

  meta_data$tmp_col <- apply(meta_data[,group_col, drop = FALSE], 1, function(row) paste(row, collapse = "_"))

  recalc_df_list <- mclapply(unique(meta_data$tmp_col), function(col_v){
    tmp_df <- .extract_IgScanDf_from_SingleCell(meta_data[meta_data$tmp_col == col_v,])

    if(nrow(tmp_df) == 0){return(NULL)}

    ## To correct ClonotypeID
    clt_dict <- aggregate(x = tmp_df$ClonotypeID[!is.na(tmp_df$ClonotypeID)], by = list(tmp_df$ClonotypeID[!is.na(tmp_df$ClonotypeID)]), FUN = length)
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
    tmp_df$ClonotypeID <- clt_dict$NewName[match(tmp_df$ClonotypeID, clt_dict$Group.1)]

    ## To correct ClonotypeVariantID
    cv_dict_all <- data.frame()
    for(c in unique(tmp_df$ClonotypeID)){
      cv_dict <- aggregate(x = tmp_df$ClonotypeVariantID[tmp_df$ClonotypeID == c], by = list(tmp_df$ClonotypeVariantID[tmp_df$ClonotypeID == c]), FUN = length)
      cv_dict <- cv_dict[order(cv_dict$x, decreasing = T),]
      cv_dict$NewName <- paste0(c, ".CV", 1:nrow(cv_dict))
      cv_dict_all <- rbind(cv_dict_all, cv_dict)
    }
    tmp_df$ClonotypeVariantID <- cv_dict_all$NewName[match(tmp_df$ClonotypeVariantID, cv_dict_all$Group.1)]

    ## To correct SubcloneID
    sbc_dict_all <- data.frame()
    for(cv in unique(tmp_df$ClonotypeVariantID)){
      sbc_dict <- aggregate(x = tmp_df$SubcloneID[tmp_df$ClonotypeVariantID == cv], by = list(tmp_df$SubcloneID[tmp_df$ClonotypeVariantID == cv]), FUN = length)
      sbc_dict <- sbc_dict[order(sbc_dict$x, decreasing = T),]
      sbc_dict$NewName <- paste0(cv, ".S", 1:nrow(sbc_dict))
      sbc_dict_all <- rbind(sbc_dict_all, sbc_dict)
    }
    tmp_df$SubcloneID <- sbc_dict_all$NewName[match(tmp_df$SubcloneID, sbc_dict_all$Group.1)]

    ## Update merged IDs
    tmp_df$igClonotypeID_tmp <- NA
    tmp_df$igClonotypeVariantID_tmp <- NA
    tmp_df$igSubcloneID_tmp <- NA
    for(cell in unique(tmp_df$igSubcloneID)){
      old_cltID_list <- strsplit(unique(tmp_df$igClonotypeID[tmp_df$igSubcloneID == cell]), "-")[[1]]
      tmp_df$igClonotypeID_tmp[tmp_df$igSubcloneID == cell] <- paste(clt_dict$NewName[match(old_cltID_list, clt_dict$Group.1)], collapse = "-")

      old_cvID_list <- strsplit(unique(tmp_df$igClonotypeVariantID[tmp_df$igSubcloneID == cell]), "-")[[1]]
      tmp_df$igClonotypeVariantID_tmp[tmp_df$igSubcloneID == cell] <- paste(cv_dict_all$NewName[match(old_cvID_list, cv_dict_all$Group.1)], collapse = "-")

      old_scID_list <- strsplit(cell, "-")[[1]]
      tmp_df$igSubcloneID_tmp[tmp_df$igSubcloneID == cell] <- paste(sbc_dict_all$NewName[match(old_scID_list, sbc_dict_all$Group.1)], collapse = "-")
    }
    tmp_df$igClonotypeID <- tmp_df$igClonotypeID_tmp
    tmp_df$igClonotypeVariantID <- tmp_df$igClonotypeVariantID_tmp
    tmp_df$igSubcloneID <- tmp_df$igSubcloneID_tmp
    tmp_df <- tmp_df[,-((ncol(tmp_df)-2):ncol(tmp_df))]

    if(!all(is.na(tmp_df$igClonotypeID_num[tmp_df$igClonotypeID_num != ""]))){
      ## Correct numeric IDs based on new merged IDs
      per_cell_clt_dict <- aggregate(x = tmp_df$igClonotypeID_num[tmp_df$igClonotypeID_num != "", drop = FALSE], by = list(tmp_df$igClonotypeID_num[tmp_df$igClonotypeID_num != ""]), FUN = length)
      per_cell_clt_dict <- per_cell_clt_dict[order(per_cell_clt_dict$x, decreasing = T),]
      per_cell_clt_dict$NewName <- paste0("C", 1:nrow(per_cell_clt_dict))
      tmp_df$igClonotypeID_num <- per_cell_clt_dict$NewName[match(tmp_df$igClonotypeID_num, per_cell_clt_dict$Group.1)]

      per_cell_cv_dict_all <- data.frame()
      per_cell_sbc_in_clt_dict_all <- data.frame()
      for(c in unique(tmp_df$igClonotypeID_num[!is.na(tmp_df$igClonotypeID_num)])){
        cv_dict <- aggregate(x = tmp_df$igClonotypeVariantID_num[tmp_df$igClonotypeID_num == c, drop = FALSE], by = list(tmp_df$igClonotypeVariantID_num[tmp_df$igClonotypeID_num == c]), FUN = length)
        cv_dict <- cv_dict[order(cv_dict$x, decreasing = T),]
        cv_dict$NewName <- paste0(c, ".CV", 1:nrow(cv_dict))
        per_cell_cv_dict_all <- rbind(per_cell_cv_dict_all, cv_dict)

        sbc_in_clt_dict <- aggregate(x = tmp_df$igSubcloneID_in_Clonotype_num[tmp_df$igClonotypeID_num == c, drop = FALSE], by = list(tmp_df$igSubcloneID_in_Clonotype_num[tmp_df$igClonotypeID_num == c]), FUN = length)
        sbc_in_clt_dict <- sbc_in_clt_dict[order(sbc_in_clt_dict$x, decreasing = T),]
        sbc_in_clt_dict$NewName <- paste0(c, ".", 1:nrow(sbc_in_clt_dict))
        per_cell_sbc_in_clt_dict_all <- rbind(per_cell_sbc_in_clt_dict_all, sbc_in_clt_dict)
      }
      tmp_df$igClonotypeVariantID_num <- per_cell_cv_dict_all$NewName[match(tmp_df$igClonotypeVariantID_num, per_cell_cv_dict_all$Group.1)]
      tmp_df$igSubcloneID_in_Clonotype_num <- per_cell_sbc_in_clt_dict_all$NewName[match(tmp_df$igSubcloneID_in_Clonotype_num, per_cell_sbc_in_clt_dict_all$Group.1)]

      per_cell_sbc_in_cv_dict_all <- data.frame()
      for(cv in unique(tmp_df$igClonotypeVariantID_num[!is.na(tmp_df$igClonotypeVariantID_num)])){
        sbc_in_cv_dict <- aggregate(x = tmp_df$igSubcloneID_in_ClonotypeVariant_num[tmp_df$igClonotypeVariantID_num == cv, drop = FALSE], by = list(tmp_df$igSubcloneID_in_ClonotypeVariant_num[tmp_df$igClonotypeVariantID_num == cv]), FUN = length)
        sbc_in_cv_dict <- sbc_in_cv_dict[order(sbc_in_cv_dict$x, decreasing = T),]
        sbc_in_cv_dict$NewName <- paste0(cv, ".S", 1:nrow(sbc_in_cv_dict))
        per_cell_sbc_in_cv_dict_all <- rbind(per_cell_sbc_in_cv_dict_all, sbc_in_cv_dict)
      }
      tmp_df$igSubcloneID_in_ClonotypeVariant_num <- per_cell_sbc_in_cv_dict_all$NewName[match(tmp_df$igSubcloneID_in_ClonotypeVariant_num, per_cell_sbc_in_cv_dict_all$Group.1)]
    }
    return(tmp_df)
  }, mc.cores = threads)

  recalc_df_list <- Filter(Negate(is.null), recalc_df_list)
  recalc_df <- do.call(rbind, recalc_df_list)

  if(class(single_cell_object)[1] == "SingleCellExperiment"){
    single_cell_object <- combine_IgScan_SingleCellExperiment(igscan_out = recalc_df, sce = single_cell_object)
  } else if(class(single_cell_object)[1] == "Seurat"){
    single_cell_object <- combine_IgScan_Seurat(igscan_out = recalc_df, seurat_object = single_cell_object)
  }
  return(single_cell_object)
}

.extract_IgScanDf_from_SingleCell <- function(meta_data) {

  chain_dictionary <- c("IGH1", "IGH2", "IGL1", "IGL2")

  igscan_out <- data.frame()

  for(sbc_all in unique(meta_data$igSubcloneID_all)){
    subclone_IDs <- strsplit(sbc_all, split = "-")[[1]]
    writen_contigs <- 1

    for(i in 1:4){
      sc <- subclone_IDs[i]
      if(sc == "NA" | is.na(sc)){next}

      chain <- chain_dictionary[i]

      barcode = rownames( meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all),] )

      append_df <- data.frame(contig_id = paste0(barcode, "_", writen_contigs), barcode = barcode)

      append_df$Raw_sequence <- NA
      append_df$Raw_VDJ_sequence <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_Raw_VDJ_sequence")])
      append_df$VDJ_sequence <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_VDJ_sequence")])
      append_df$IgBlast_Germline_alignment <- NA
      append_df$VDJ_sequence_correctedCDR3 <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_VDJ_sequence_correctedCDR3")])
      append_df$VDJ_sequence_correctedCDR3_aa <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_VDJ_sequence_correctedCDR3_aa")])
      append_df$Consensus_Germline <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_Consensus_Germline")])
      append_df$Consensus_Germline_aa <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_Consensus_Germline_aa")])
      append_df$VDJ_genes <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_VDJ_genes")])
      append_df$C_gene <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_C_gene")])
      append_df$Functionality <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_Functionality")])
      append_df$Junction_aa <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_Junction_aa")])
      append_df$Junction_lenght <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_Junction_lenght")])
      append_df$V_identity <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_V_identity")])
      append_df$VDJ_positions <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_VDJ_positions")])
      append_df$V_length <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_V_length")])
      append_df$InDels <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_InDels")])
      if(paste0(chain, "_CLL_Stereotype_Subsets") %in% colnames(meta_data)){append_df$CLL_Stereotype_Subsets <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_CLL_Stereotype_Subsets")])}

      append_df$ClonotypeID <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_ClonotypeID")])
      append_df$Clonotype_Consensus_CDR3aa <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_Clonotype_Consensus_CDR3aa")])

      if(paste0(chain, "_Clonotype_CLL_Stereotype_Subsets") %in% colnames(meta_data)){Clonotype_CLL_Stereotype_Subsets <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_Clonotype_CLL_Stereotype_Subsets")])}

      append_df$ClonotypeVariantID <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_ClonotypeVariantID")])
      append_df$SubcloneID <- unique(meta_data[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all), paste0(chain, "_SubcloneID")])

      append_df$completeBCR <- unique(meta_data$completeBCR[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igClonotypeID_num <- unique(meta_data$igClonotypeID_num[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igClonotypeID <- unique(meta_data$igClonotypeID[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igClonotypeVariantID_num <- unique(meta_data$igClonotypeVariantID_num[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igClonotypeVariantID <- unique(meta_data$igClonotypeVariantID[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igSubcloneID_in_ClonotypeVariant_num <- unique(meta_data$igSubcloneID_in_ClonotypeVariant_num[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igSubcloneID_in_Clonotype_num <- unique(meta_data$igSubcloneID_in_Clonotype_num[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igSubcloneID <- unique(meta_data$igSubcloneID[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igRaw_VDJ_sequence <- unique(meta_data$igRaw_VDJ_sequence[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igVDJ_sequence <- unique(meta_data$igVDJ_sequence[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igVDJ_sequence_aa <- unique(meta_data$igVDJ_sequence_aa[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igClonotype_Consensus_Germline <- unique(meta_data$igClonotype_Consensus_Germline[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igClonotype_Consensus_Germline_aa <- unique(meta_data$igClonotype_Consensus_Germline_aa[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igVDJ_positions <- unique(meta_data$igVDJ_positions[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igInDels <- unique(meta_data$igInDels[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      append_df$igClonotype_Consensus_CDR3aa <- unique(meta_data$igClonotype_Consensus_CDR3aa[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])
      if("igCLL_Stereotype_Subsets" %in% colnames(meta_data)){igCLL_Stereotype_Subsets <- unique(meta_data$igCLL_Stereotype_Subsets[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)])}
      append_df$SampleID <- as.character(unique(meta_data$orig.ident[meta_data$igSubcloneID_all == sbc_all & !is.na(meta_data$igSubcloneID_all)]))

      igscan_out <- rbind(igscan_out, append_df)
      writen_contigs <- writen_contigs + 1
    }
  }
  return(igscan_out)
}
