#' Combine IgScan Output with Seurat Object
#'
#' @description This function integrates immunogenetic data from IgScan with a Seurat
#' object. It maps clonotype-related information and sequence data to the metadata of
#' the Seurat object based on matching barcodes.
#'
#' @param igscan_out A data frame containing the output from IgScan for a single cell dataset.
#' @param seurat_object A Seurat object where metadata will be updated with the corresponding immunogenetic information from `igscan_out`.
#'
#' @return A Seurat object with updated metadata, including the IgScan annotation organized by cell barcode.
#'
#' @export
#'
#' @importFrom qs qread qsave
#' @importFrom stringr str_count
#' @import Seurat
#' @import SeuratObject
#'
#' @examples
#' \dontrun{
#' updated_seurat_object <- combine_IgScan_Seurat(igscan_out, seurat_object)
#' }
#'
combine_IgScan_Seurat <- function(igscan_out, seurat_object){

  if(!"barcode" %in% colnames(igscan_out)){igscan_out$barcode <- sapply(igscan_out$contig_id, function(x) strsplit(x, "_")[[1]][1]}

  seurat_object@meta.data$completeBCR <- igscan_out$completeBCR[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igClonotypeID <- igscan_out$igClonotypeID[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igClonotypeID_num <- igscan_out$igClonotypeID_num[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igClonotypeVariantID <- igscan_out$igClonotypeVariantID[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igClonotypeVariantID_num <- igscan_out$igClonotypeVariantID_num[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igSubcloneID <- igscan_out$igSubcloneID[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igSubcloneID_in_ClonotypeVariant_num <- igscan_out$igSubcloneID_in_ClonotypeVariant_num[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igSubcloneID_in_Clonotype_num <- igscan_out$igSubcloneID_in_Clonotype_num[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igRaw_VDJ_sequence <- igscan_out$igRaw_VDJ_sequence[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igVDJ_sequence <- igscan_out$igVDJ_sequence[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igVDJ_sequence_aa <- igscan_out$igVDJ_sequence_aa[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igClonotype_Consensus_Germline <- igscan_out$igClonotype_Consensus_Germline[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igClonotype_Consensus_Germline_aa <- igscan_out$igClonotype_Consensus_Germline_aa[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igVDJ_positions <- igscan_out$igVDJ_positions[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igInDels <- igscan_out$igInDels[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  seurat_object@meta.data$igClonotype_Consensus_CDR3aa <- igscan_out$igClonotype_Consensus_CDR3aa[match(rownames(seurat_object@meta.data), igscan_out$barcode)]
  if("igCLL_Stereotype_Subsets" %in% colnames(igscan_out)){seurat_object@meta.data$igCLL_Stereotype_Subsets <- igscan_out$igCLL_Stereotype_Subsets[match(rownames(seurat_object@meta.data), igscan_out$barcode)]}

  seurat_object@meta.data$igSubcloneID_all <- sapply(seurat_object@meta.data$igSubcloneID, function(x){
    if(is.na(x)){
      NA
    }else{
      ighCo <- 2 - str_count(x, "IGH")
      iglCo <- 2 - str_count(x, "IGK|IGL")
      rearr <- strsplit(x, "-")[[1]]
      paste0(ifelse(ighCo != 2, paste0(paste0(rearr[1:2-ighCo], collapse="-"),"-"), ""),
             ifelse(ighCo != 0, paste0(paste0(rep("NA", ighCo), collapse="-"), "-"), ""),
             ifelse(iglCo != 2, paste0(rearr[(3-ighCo):length(rearr)], collapse="-"), "NA-NA"),
             ifelse(iglCo != 0, rep("-NA", iglCo), ""))
    }
  })

  split_pos <- 1
  for(n in c("IGH1", "IGH2", "IGL1", "IGL2")){

    subclone <- sapply(seurat_object@meta.data$igSubcloneID_all, function(x) strsplit(x, split = "-")[[1]][split_pos] )

    ## Immunogenetic data
    seurat_object@meta.data[[paste0(n, "_VDJ_genes")]] <- igscan_out$VDJ_genes[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_C_gene")]] <- igscan_out$C_gene[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_Junction_aa")]] <- igscan_out$Junction_aa[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_Junction_lenght")]] <- igscan_out$Junction_lenght[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_Functionality")]] <- igscan_out$Functionality[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_V_length")]] <- igscan_out$V_length[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_V_identity")]] <- igscan_out$V_identity[match(subclone, igscan_out$SubcloneID)]
    if("CLL_Stereotype_Subsets" %in% colnames(igscan_out)){seurat_object@meta.data[[paste0(n, "_CLL_Stereotype_Subsets")]] <- igscan_out$CLL_Stereotype_Subsets[match(subclone, igscan_out$SubcloneID)]}

    ## Clonotype level data
    seurat_object@meta.data[[paste0(n, "_ClonotypeID")]] <- igscan_out$ClonotypeID[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_Clonotype_Consensus_CDR3aa")]] <- igscan_out$Clonotype_Consensus_CDR3aa[match(subclone, igscan_out$SubcloneID)]
    if("Clonotype_CLL_Stereotype_Subsets" %in% colnames(igscan_out)){seurat_object@meta.data[[paste0(n, "_Clonotype_CLL_Stereotype_Subsets")]] <- igscan_out$Clonotype_CLL_Stereotype_Subsets[match(subclone, igscan_out$SubcloneID)]}
    seurat_object@meta.data[[paste0(n, "_ClonotypeVariantID")]] <- igscan_out$ClonotypeVariantID[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_SubcloneID")]] <- igscan_out$SubcloneID[match(subclone, igscan_out$SubcloneID)]

    ## Sequence data
    seurat_object@meta.data[[paste0(n, "_Raw_VDJ_sequence")]] <- igscan_out$Raw_VDJ_sequence[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_VDJ_sequence")]] <- igscan_out$VDJ_sequence[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_VDJ_sequence_correctedCDR3")]] <- igscan_out$VDJ_sequence_correctedCDR3[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_VDJ_sequence_correctedCDR3_aa")]] <- igscan_out$VDJ_sequence_correctedCDR3_aa[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_Consensus_Germline")]] <- igscan_out$Consensus_Germline[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_Consensus_Germline_aa")]] <- igscan_out$Consensus_Germline_aa[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_VDJ_positions")]] <- igscan_out$VDJ_positions[match(subclone, igscan_out$SubcloneID)]
    seurat_object@meta.data[[paste0(n, "_InDels")]] <- igscan_out$InDels[match(subclone, igscan_out$SubcloneID)]

    split_pos <- split_pos + 1
  }

  seurat_object@meta.data[seurat_object@meta.data == ""] <- NA

  return(seurat_object)
}
