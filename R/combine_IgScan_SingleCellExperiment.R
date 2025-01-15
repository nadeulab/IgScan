#' Combine IgScan Output with Single Cell Experiment
#'
#' @description This function integrates immunogenetic data from IgScan with a Single
#' Cell Experiment (SCE) object. It maps clonotype-related information and sequence data to
#' the colData of the SCE object based on matching barcodes.
#'
#' @param igscan_out A data frame containing the output from IgScan for a single cell dataset.
#' @param seurat_object A SCE object where colData will be updated with the corresponding immunogenetic information from `igscan_out`.
#'
#' @return A SCE object with updated colData, including the IgScan annotation organized by cell barcode.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' updated_sce_object <- combine_IgScan_Seurat(igscan_out, sce_object)
#' }
#'
combine_IgScan_SingleCellExperiment <- function(igscan_out, sce){

  igscan_out$barcode <- sapply(igscan_out$contig_id, function(x) strsplit(x, "_")[[1]][1])

  sce@colData$completeBCR <- igscan_out$completeBCR[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igClonotypeID <- igscan_out$igClonotypeID[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igClonotypeID_num <- igscan_out$igClonotypeID_num[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igClonotypeVariantID <- igscan_out$igClonotypeVariantID[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igClonotypeVariantID_num <- igscan_out$igClonotypeVariantID_num[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igSubcloneID <- igscan_out$igSubcloneID[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igSubcloneID_in_ClonotypeVariant_num <- igscan_out$igSubcloneID_in_ClonotypeVariant_num[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igSubcloneID_in_Clonotype_num <- igscan_out$igSubcloneID_in_Clonotype_num[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igVDJ_sequence <- igscan_out$igVDJ_sequence[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igVDJ_sequence_aa <- igscan_out$igVDJ_sequence_aa[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igClonotype_Consensus_Germline <- igscan_out$igClonotype_Consensus_Germline[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igClonotype_Consensus_Germline_aa <- igscan_out$igClonotype_Consensus_Germline_aa[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igVDJ_positions <- igscan_out$igVDJ_positions[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igInDels <- igscan_out$igInDels[match(rownames(sce@colData), igscan_out$barcode)]
  sce@colData$igClonotype_Consensus_CDR3aa <- igscan_out$igClonotype_Consensus_CDR3aa[match(rownames(sce@colData), igscan_out$barcode)]
  if("ig_CLL_Stereotype_Subsets" %in% colnames(igscan_out)){sce@colData$ig_CLL_Stereotype_Subsets <- igscan_out$ig_CLL_Stereotype_Subsets[match(rownames(sce@colData), igscan_out$barcode)]}

  sce@colData$igSubcloneID_all <- sapply(sce@colData$igSubcloneID, function(x){
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

    subclone <- sapply(sce@colData$igSubcloneID_all, function(x) strsplit(x, split = "-")[[1]][split_pos] )

    ## Immunogenetic data
    sce@colData[[paste0(n, "_VDJ_genes")]] <- igscan_out$VDJ_genes[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_C_gene")]] <- igscan_out$C_gene[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_Junction_aa")]] <- igscan_out$Junction_aa[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_Junction_lenght")]] <- igscan_out$Junction_lenght[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_Functionality")]] <- igscan_out$Functionality[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_V_length")]] <- igscan_out$V_length[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_V_identity")]] <- igscan_out$V_identity[match(subclone, igscan_out$SubcloneID)]
    if("CLL_Stereotype_Subsets" %in% colnames(igscan_out)){sce@colData[[paste0(n, "_CLL_Stereotype_Subsets")]] <- igscan_out$CLL_Stereotype_Subsets[match(subclone, igscan_out$SubcloneID)]}

    ## Clonotype level data
    sce@colData[[paste0(n, "_ClonotypeID")]] <- igscan_out$ClonotypeID[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_Clonotype_Consensus_CDR3aa")]] <- igscan_out$Clonotype_Consensus_CDR3aa[match(subclone, igscan_out$SubcloneID)]
    if("Clonotype_CLL_Stereotype_Subsets" %in% colnames(igscan_out)){sce@colData[[paste0(n, "_Clonotype_CLL_Stereotype_Subsets")]] <- igscan_out$Clonotype_CLL_Stereotype_Subsets[match(subclone, igscan_out$SubcloneID)]}
    sce@colData[[paste0(n, "_ClonotypeVariantID")]] <- igscan_out$ClonotypeVariantID[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_SubcloneID")]] <- igscan_out$SubcloneID[match(subclone, igscan_out$SubcloneID)]

    ## Sequence data
    sce@colData[[paste0(n, "_VDJ_sequence")]] <- igscan_out$VDJ_sequence[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_IgBlast_Germline_alignment")]] <- igscan_out$IgBlast_Germline_alignment[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_VDJ_sequence_correctedCDR3")]] <- igscan_out$VDJ_sequence_correctedCDR3[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_VDJ_sequence_correctedCDR3_aa")]] <- igscan_out$VDJ_sequence_correctedCDR3_aa[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_Consensus_Germline")]] <- igscan_out$Consensus_Germline[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_Consensus_Germline_aa")]] <- igscan_out$Consensus_Germline_aa[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_VDJ_positions")]] <- igscan_out$VDJ_positions[match(subclone, igscan_out$SubcloneID)]
    sce@colData[[paste0(n, "_InDels")]] <- igscan_out$InDels[match(subclone, igscan_out$SubcloneID)]

    split_pos <- split_pos + 1
  }

  col_data <- as.data.frame(colData(sce))
  col_data[col_data == ""] <- NA
  colData(sce) <- DataFrame(col_data)

  return(sce)
}
