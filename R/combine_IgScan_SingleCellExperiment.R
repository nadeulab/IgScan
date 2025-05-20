#' Combine IgScan Output with a Single Cell Experiment Object
#'
#' @description This function integrates immunogenetic data from IgScan with a Single
#' Cell Experiment (SCE) object. It maps clonotype-related information and sequence data to
#' the colData of the SCE object based on matching barcodes.
#'
#' @param igscan_out A data frame containing the output from IgScan for a single cell dataset.
#' @param sce A SCE object where colData will be updated with the corresponding immunogenetic information from `igscan_out`.
#' @param sce_sample_col A vector with the name of the column or columns containing the sample identifiers
#' in the metadata of the SCE object. There is no default value for this parameter.
#' @param igscan_sample_col A vector with the name of the column or columns containing the sample identifiers
#' in the IgScan annotation data frame. Default is 'SampleID'.
#' @param threads The number of threads to be used. Default is 1.

#' @return A SCE object with updated colData, including the IgScan annotation organized by cell barcode.
#'
#' @export
#'
#' @importFrom qs qread qsave
#' @importFrom stringr str_count
#' @import SingleCellExperiment
#' @importFrom singleCellTK combineSCE
#'
#' @examples
#' \dontrun{
#' updated_sce_object <- combine_IgScan_Seurat(igscan_out, sce_object, sce_sample_col = "Sample", igscan_sample_col = "SampleID", threads = 4)
#' }
#'
combine_IgScan_SingleCellExperiment <- function(igscan_out, sce, sce_sample_col, igscan_sample_col = "SampleID", threads = 1){

  if(!all(sce_sample_col %in% colnames(colData(sce)))){stop(paste0("\nUnknown column (", sce_sample_col[!sce_sample_col %in% colnames(colData(sce))], ")! Please, set a valid column name."))}
  if(!all(igscan_sample_col %in% colnames(igscan_out))){stop(paste0("\nUnknown column (", igscan_sample_col[!igscan_sample_col %in% colnames(igscan_out)], ")! Please, set a valid column name."))}

  colData(sce)$tmp_col <- apply(colData(sce)[,sce_sample_col, drop = FALSE], 1, function(row) paste(row, collapse = "_"))
  igscan_out$tmp_col <- apply(igscan_out[,igscan_sample_col, drop = FALSE], 1, function(row) paste(row, collapse = "_"))

  if(!"barcode" %in% colnames(igscan_out)){igscan_out$barcode <- sapply(igscan_out$contig_id, function(x) strsplit(x, "_")[[1]][1])}

  tmp_object_list <- mclapply(unique(colData(sce)$tmp_col), function(sample_id){

    tmp_sce <- sce[, sce$tmp_col == sample_id]
    tmp_igscan <- igscan_out[igscan_out$tmp_col == sample_id, ]

    if(nrow(tmp_igscan) == 0 | ncol(tmp_sce) == 0){return(NULL)}

    colData(tmp_sce)$completeBCR <- tmp_igscan$completeBCR[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igClonotypeID <- tmp_igscan$igClonotypeID[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igClonotypeID_num <- tmp_igscan$igClonotypeID_num[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igClonotypeVariantID <- tmp_igscan$igClonotypeVariantID[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igClonotypeVariantID_num <- tmp_igscan$igClonotypeVariantID_num[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igSubcloneID <- tmp_igscan$igSubcloneID[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igSubcloneID_in_ClonotypeVariant_num <- tmp_igscan$igSubcloneID_in_ClonotypeVariant_num[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igSubcloneID_in_Clonotype_num <- tmp_igscan$igSubcloneID_in_Clonotype_num[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igRaw_VDJ_sequence <- tmp_igscan$igRaw_VDJ_sequence[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igVDJ_sequence <- tmp_igscan$igVDJ_sequence[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igVDJ_sequence_aa <- tmp_igscan$igVDJ_sequence_aa[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igClonotype_Consensus_Germline <- tmp_igscan$igClonotype_Consensus_Germline[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igClonotype_Consensus_Germline_aa <- tmp_igscan$igClonotype_Consensus_Germline_aa[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igVDJ_positions <- tmp_igscan$igVDJ_positions[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igInDels <- tmp_igscan$igInDels[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    colData(tmp_sce)$igClonotype_Consensus_CDR3aa <- tmp_igscan$igClonotype_Consensus_CDR3aa[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]
    if("igCLL_Stereotype_Subsets" %in% colnames(tmp_igscan)){colData(tmp_sce)$igCLL_Stereotype_Subsets <- tmp_igscan$igCLL_Stereotype_Subsets[match(rownames(colData(tmp_sce)), tmp_igscan$barcode)]}

    colData(tmp_sce)$igSubcloneID_all <- sapply(colData(tmp_sce)$igSubcloneID, function(x){
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

      subclone <- sapply(colData(tmp_sce)$igSubcloneID_all, function(x) strsplit(x, split = "-")[[1]][split_pos])

      ## Immunogenetic data
      colData(tmp_sce)[[paste0(n, "_VDJ_genes")]] <- tmp_igscan$VDJ_genes[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_C_gene")]] <- tmp_igscan$C_gene[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_Junction_aa")]] <- tmp_igscan$Junction_aa[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_Junction_lenght")]] <- tmp_igscan$Junction_lenght[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_Functionality")]] <- tmp_igscan$Functionality[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_V_length")]] <- tmp_igscan$V_length[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_V_identity")]] <- tmp_igscan$V_identity[match(subclone, tmp_igscan$SubcloneID)]
      if("CLL_Stereotype_Subsets" %in% colnames(tmp_igscan)){colData(tmp_sce)[[paste0(n, "_CLL_Stereotype_Subsets")]] <- tmp_igscan$CLL_Stereotype_Subsets[match(subclone, tmp_igscan$SubcloneID)]}

      ## Clonotype level data
      colData(tmp_sce)[[paste0(n, "_ClonotypeID")]] <- tmp_igscan$ClonotypeID[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_Clonotype_Consensus_CDR3aa")]] <- tmp_igscan$Clonotype_Consensus_CDR3aa[match(subclone, tmp_igscan$SubcloneID)]
      if("Clonotype_CLL_Stereotype_Subsets" %in% colnames(tmp_igscan)){colData(tmp_sce)[[paste0(n, "_Clonotype_CLL_Stereotype_Subsets")]] <- tmp_igscan$Clonotype_CLL_Stereotype_Subsets[match(subclone, tmp_igscan$SubcloneID)]}
      colData(tmp_sce)[[paste0(n, "_ClonotypeVariantID")]] <- tmp_igscan$ClonotypeVariantID[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_SubcloneID")]] <- tmp_igscan$SubcloneID[match(subclone, tmp_igscan$SubcloneID)]

      ## Sequence data
      colData(tmp_sce)[[paste0(n, "_Raw_VDJ_sequence")]] <- tmp_igscan$Raw_VDJ_sequence[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_VDJ_sequence")]] <- tmp_igscan$VDJ_sequence[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_VDJ_sequence_correctedCDR3")]] <- tmp_igscan$VDJ_sequence_correctedCDR3[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_VDJ_sequence_correctedCDR3_aa")]] <- tmp_igscan$VDJ_sequence_correctedCDR3_aa[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_Consensus_Germline")]] <- tmp_igscan$Consensus_Germline[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_Consensus_Germline_aa")]] <- tmp_igscan$Consensus_Germline_aa[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_VDJ_positions")]] <- tmp_igscan$VDJ_positions[match(subclone, tmp_igscan$SubcloneID)]
      colData(tmp_sce)[[paste0(n, "_InDels")]] <- tmp_igscan$InDels[match(subclone, tmp_igscan$SubcloneID)]

      split_pos <- split_pos + 1
    }

    tmp_col_data <- as.data.frame(colData(tmp_sce))
    tmp_col_data[tmp_col_data == ""] <- NA
    colData(tmp_sce) <- DataFrame(tmp_col_data)

    return(tmp_sce)

  }, mc.cores = threads)

  tmp_object_list <- Filter(Negate(is.null), tmp_object_list)

  if(length(tmp_object_list) == 0){stop("No cells found in none of the sample identifiers specified. Please, ensure that sample identifiers are valid.\n")}

  combined_sce_object <- combineSCE(tmp_object_list, combined = TRUE)

  combined_col_data <- as.data.frame(colData(combined_sce_object))
  combined_col_data <- combined_col_data[, colnames(combined_col_data) != "tmp_col"]
  colData(combined_sce_object) <- DataFrame(combined_col_data)

  return(combined_sce_object)
}
