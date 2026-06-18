#' Combine IgScan Output with a Seurat Object
#'
#' @description This function integrates immunogenetic data from IgScan with a Seurat
#' object. It maps clonotype-related information and sequence data to the metadata of
#' the Seurat object based on matching barcodes.
#'
#' @param igscan_out A data frame containing the output from IgScan for a single cell dataset.
#' @param seurat_object A Seurat object where metadata will be updated with the corresponding immunogenetic information from `igscan_out`.
#' @param seurat_sample_col A vector with the name of the column or columns containing the sample identifiers
#' in the metadata of the Seurat object. Default is 'orig.ident'.
#' @param igscan_sample_col A vector with the name of the column or columns containing the sample identifiers
#' in the IgScan annotation data frame. Default is 'SampleID'.
#' @param threads The number of threads to be used. Default is 1.
#'
#' @note
#' When several sample columns are passed to the function, a temporary ID is generated (e.g. Variable1_Variable2).
#' Importantly, sample identifiers must be identical and coincide between the Seurat object and IgScan annotation file.
#'
#' @return A Seurat object with updated metadata, including the IgScan annotation organized by cell barcode.
#'
#' @export
#'
#' @importFrom stringr str_count
#' @import Seurat
#' @import SeuratObject
#' @importFrom parallel mclapply
#'
#' @examples
#' \dontrun{
#' updated_seurat_object <- combine_IgScan_Seurat(igscan_out, seurat_object)
#' }
#'
combine_IgScan_Seurat <- function(igscan_out, seurat_object, seurat_sample_col = "orig.ident", igscan_sample_col = "SampleID", threads = 1){

  if(!all(seurat_sample_col %in% colnames(seurat_object@meta.data))){stop(paste0("\nUnknown column (", seurat_sample_col[!seurat_sample_col %in% colnames(seurat_object@meta.data)], ")! Please, set a valid column name."))}
  if(!all(igscan_sample_col %in% colnames(igscan_out))){stop(paste0("\nUnknown column (", igscan_sample_col[!igscan_sample_col %in% colnames(igscan_out)], ")! Please, set a valid column name."))}

  seurat_object@meta.data$tmp_col <- apply(seurat_object@meta.data[,seurat_sample_col, drop = FALSE], 1, function(row) paste(row, collapse = "_"))
  igscan_out$tmp_col <- apply(igscan_out[,igscan_sample_col, drop = FALSE], 1, function(row) paste(row, collapse = "_"))

  if(!"barcode" %in% colnames(igscan_out)){igscan_out$barcode <- sapply(igscan_out$contig_id, function(x) strsplit(x, "_")[[1]][1])}

  tmp_object_list <- mclapply(unique(seurat_object@meta.data$tmp_col), function(sample_id){

    tmp_seurat <- subset(seurat_object, subset = tmp_col == sample_id)
    tmp_igscan <- igscan_out[igscan_out$tmp_col == sample_id, ]

    if(ncol(tmp_seurat) == 0){
      return(NULL)
    } else if(nrow(tmp_igscan) == 0){
      return(tmp_seurat@meta.data)
    }

    tmp_seurat@meta.data$completeBCR <- tmp_igscan$completeBCR[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igClonotypeID <- tmp_igscan$igClonotypeID[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igClonotypeID_num <- tmp_igscan$igClonotypeID_num[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igClonotypeVariantID <- tmp_igscan$igClonotypeVariantID[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igClonotypeVariantID_num <- tmp_igscan$igClonotypeVariantID_num[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igSubcloneID <- tmp_igscan$igSubcloneID[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igSubcloneID_in_ClonotypeVariant_num <- tmp_igscan$igSubcloneID_in_ClonotypeVariant_num[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igSubcloneID_in_Clonotype_num <- tmp_igscan$igSubcloneID_in_Clonotype_num[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igRaw_VDJ_sequence <- tmp_igscan$igRaw_VDJ_sequence[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igVDJ_sequence <- tmp_igscan$igVDJ_sequence[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igVDJ_sequence_aa <- tmp_igscan$igVDJ_sequence_aa[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igClonotype_Consensus_Germline <- tmp_igscan$igClonotype_Consensus_Germline[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igClonotype_Consensus_Germline_aa <- tmp_igscan$igClonotype_Consensus_Germline_aa[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igVDJ_positions <- tmp_igscan$igVDJ_positions[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igInDels <- tmp_igscan$igInDels[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    tmp_seurat@meta.data$igClonotype_Consensus_CDR3aa <- tmp_igscan$igClonotype_Consensus_CDR3aa[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]
    if("igCLL_Stereotype_Subsets" %in% colnames(tmp_igscan)){tmp_seurat@meta.data$igCLL_Stereotype_Subsets <- tmp_igscan$igCLL_Stereotype_Subsets[match(rownames(tmp_seurat@meta.data), tmp_igscan$barcode)]}

    tmp_seurat@meta.data$igSubcloneID_all <- sapply(tmp_seurat@meta.data$igSubcloneID, function(x){
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

      subclone <- sapply(tmp_seurat@meta.data$igSubcloneID_all, function(x) strsplit(x, split = "-")[[1]][split_pos])

      ## Immunogenetic data
      tmp_seurat@meta.data[[paste0(n, "_VDJ_genes")]] <- tmp_igscan$VDJ_genes[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_C_gene")]] <- tmp_igscan$C_gene[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_Junction_aa")]] <- tmp_igscan$Junction_aa[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_Junction_length")]] <- tmp_igscan$Junction_length[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_Functionality")]] <- tmp_igscan$Functionality[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_V_length")]] <- tmp_igscan$V_length[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_V_identity")]] <- tmp_igscan$V_identity[match(subclone, tmp_igscan$SubcloneID)]
      if("CLL_Stereotype_Subsets" %in% colnames(tmp_igscan)){tmp_seurat@meta.data[[paste0(n, "_CLL_Stereotype_Subsets")]] <- tmp_igscan$CLL_Stereotype_Subsets[match(subclone, tmp_igscan$SubcloneID)]}

      ## Clonotype level data
      tmp_seurat@meta.data[[paste0(n, "_ClonotypeID")]] <- tmp_igscan$ClonotypeID[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_Clonotype_Consensus_CDR3aa")]] <- tmp_igscan$Clonotype_Consensus_CDR3aa[match(subclone, tmp_igscan$SubcloneID)]
      if("Clonotype_CLL_Stereotype_Subsets" %in% colnames(tmp_igscan)){tmp_seurat@meta.data[[paste0(n, "_Clonotype_CLL_Stereotype_Subsets")]] <- tmp_igscan$Clonotype_CLL_Stereotype_Subsets[match(subclone, tmp_igscan$SubcloneID)]}
      tmp_seurat@meta.data[[paste0(n, "_ClonotypeVariantID")]] <- tmp_igscan$ClonotypeVariantID[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_SubcloneID")]] <- tmp_igscan$SubcloneID[match(subclone, tmp_igscan$SubcloneID)]

      ## Sequence data
      tmp_seurat@meta.data[[paste0(n, "_Raw_sequence")]] <- tmp_igscan$Raw_sequence[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_Raw_VDJ_sequence")]] <- tmp_igscan$Raw_VDJ_sequence[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_VDJ_sequence")]] <- tmp_igscan$VDJ_sequence[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_IgBlast_Germline_alignment")]] <- tmp_igscan$IgBlast_Germline_alignment[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_VDJ_sequence_correctedCDR3")]] <- tmp_igscan$VDJ_sequence_correctedCDR3[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_VDJ_sequence_correctedCDR3_aa")]] <- tmp_igscan$VDJ_sequence_correctedCDR3_aa[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_Consensus_Germline")]] <- tmp_igscan$Consensus_Germline[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_Consensus_Germline_aa")]] <- tmp_igscan$Consensus_Germline_aa[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_VDJ_positions")]] <- tmp_igscan$VDJ_positions[match(subclone, tmp_igscan$SubcloneID)]
      tmp_seurat@meta.data[[paste0(n, "_InDels")]] <- tmp_igscan$InDels[match(subclone, tmp_igscan$SubcloneID)]

      split_pos <- split_pos + 1
    }

    tmp_seurat@meta.data[tmp_seurat@meta.data == ""] <- NA

    return(tmp_seurat@meta.data)

  }, mc.cores = threads)

  tmp_object_list <- Filter(Negate(is.null), tmp_object_list)

  if(length(tmp_object_list) == 0){stop("No cells found in none of the sample identifiers specified. Please, ensure that sample identifiers are valid.\n")}

  combined_seurat_object_metadata <- do.call(rbind, tmp_object_list)
  seurat_object@meta.data <- combined_seurat_object_metadata[Cells(seurat_object),]
  seurat_object@meta.data <- seurat_object@meta.data[, colnames(seurat_object@meta.data) != "tmp_col"]

  return(seurat_object)
}
