#' Rescue cells with incomplete BCR (single chains)
#'
#' This function identifies and rescues cells with incomplete BCR information (labeled
#' as "Single_chain_1" or "Single_chain_2") by attempting to infer and assign clonotype
#' and subclonotype IDs based on related cells that have complete BCR information. It
#' supports both `SingleCellExperiment` and `Seurat` objects. After rescuing single-chain
#' entries, the function also recalculates clonotype-related identifiers.
#'
#' @param single_cell_object A single-cell object combined with IgScan annotation
#' in either `SingleCellExperiment` or `Seurat` format. This function also supports
#' IgScan data.frames coming from the analysis of single cell data.
#' @param group_col A vector with the name of the column or columns containing the grouping variable/s.
#' Default is 'orig.ident', thus rescuing and recalculating BCR IDs based in each sample independently.
#' @param threads The number of threads to perform single_chain cells rescuing. Default is 1.
#' @param relaxed_rescue Logical value indicating whether to apply a relaxed clonotype rescue mode.
#' The default (recommended) is `FALSE`, requiring exact V(D)J nucleotide sequence matches to rescue
#' single-chain cells. If set to `TRUE`, which we recommend only for the analysis of Mission Bio
#' single-cell V(D)J data, rescue is performed at the clonotype level rather than by exact nucleotide
#' sequence identity.
#'
#' @return The input single cell object (or IgScan data.frame) with updated
#' metadata, where the `completeBCR`field in some "single chain" cells may now be marked as "Yes_rescue".
#'
#' @export
#'
#' @import Seurat
#' @import SeuratObject
#' @import SingleCellExperiment
#' @importFrom parallel mclapply
#'
#' @examples
#' \dontrun{
#' sce <- rescue_single_chain_cells(sce, group_col = "sample", threads = 4)
#' seurat <- rescue_single_chain_cells(seurat, group_col = "orig.ident", threads = 2)
#' igscan_df <- rescue_single_chain_cells(igscan_df, group_col = "SampleID", threads = 2)
#' }
#'
rescue_single_chain_cells <- function(single_cell_object, group_col = "orig.ident", threads = 1, relaxed_rescue = FALSE){

  if(class(single_cell_object)[1] == "SingleCellExperiment"){
    single_cell_object <- .rescue_single_chain_sce(single_cell_object, group_col, threads, relaxed_rescue)
    single_cell_object <- recalculate_IDs_single_cell(single_cell_object, group_col, threads)

  } else if(class(single_cell_object)[1] == "Seurat"){
    single_cell_object <- .rescue_single_chain_seurat(single_cell_object, group_col, threads, relaxed_rescue)
    single_cell_object <- recalculate_IDs_single_cell(single_cell_object, group_col, threads)

  } else if(class(single_cell_object)[1] == "data.frame"){
    single_cell_object <- .rescue_single_chain_igscandf(single_cell_object, group_col, threads, relaxed_rescue)
    single_cell_object <- recalculate_IDs_single_cell(single_cell_object, group_col, threads)
  }
  return(single_cell_object)
}

.rescue_single_chain_sce <- function(single_cell_object, group_col, threads, relaxed_rescue){

  if(!all(group_col %in% colnames(colData(single_cell_object)))) {stop(paste0("\nUnknown column (", group_col[!group_col %in% colnames(colData(single_cell_object))], ") selected for BCR ID recalculation! Please, set a valid column name."))}
  colData(single_cell_object)$tmp_col <- apply(as.data.frame(colData(single_cell_object)[, group_col, drop = FALSE]), 1, function(row) paste(row, collapse = "_"))

  rescue_df_list <- parallel::mclapply(unique(tmp_col), function(col_v){

    col_ids <- which(colData(single_cell_object)$tmp_col == col_v)
    tmp <- single_cell_object[, col_ids]

    if(ncol(tmp) == 0){return(NULL)}

    cols_dict <- c("completeBCR", "igClonotypeID", "igClonotypeID_num", "igClonotypeVariantID",
                   "igClonotypeVariantID_num", "igSubcloneID", "igSubcloneID_in_ClonotypeVariant_num",
                   "igSubcloneID_in_Clonotype_num")

    meta <- as.data.frame(colData(tmp))

    dict_yes <- meta[meta$completeBCR %in% "Yes", cols_dict, drop = FALSE]
    dict_yes <- dict_yes[!duplicated(dict_yes),]

    dict_single <- meta[meta$completeBCR %in% c("Single_chain_1", "Single_chain_2"), cols_dict, drop = FALSE]
    dict_single <- dict_single[!duplicated(dict_single),]

    if(relaxed_rescue){
      for(row in 1:nrow(dict_single)){
        clt_single <- dict_single$igClonotypeID[row]
        resc_df <- dict_yes[grepl(paste0("^", clt_single, "-"), dict_yes$igClonotypeID) | grepl(paste0("-", clt_single, "$"), dict_yes$igClonotypeID) | grepl(paste0("-", clt_single, "-"), dict_yes$igClonotypeID),]

        if(length(unique(resc_df$igClonotypeID_num)) > 1 | nrow(resc_df) == 0){ next }

        selection <- meta$igClonotypeID == clt_single & meta$completeBCR %in% c("Single_chain_1", "Single_chain_2")
        meta$igClonotypeID_num[selection] <- unique(resc_df$igClonotypeID_num)
        meta$igClonotypeVariantID_num[selection] <- names(which.min(sapply(unique(resc_df$igClonotypeVariantID_num), function(x) as.numeric(strsplit(x, "\\.CV")[[1]][2]))))
        meta$igSubcloneID_in_ClonotypeVariant_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_ClonotypeVariant_num), function(x) as.numeric(strsplit(x, "\\.S")[[1]][2]))))
        meta$igSubcloneID_in_Clonotype_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_Clonotype_num), function(x) as.numeric(strsplit(x, "\\.")[[1]][2]))))
        meta$completeBCR[selection] <- "Yes_rescue"
      }
    } else{
      for(row in 1:nrow(dict_single)){
        sbc_single <- dict_single$igSubcloneID[row]
        resc_df <- dict_yes[grepl(paste0("^", sbc_single, "-"), dict_yes$igSubcloneID) | grepl(paste0("-", sbc_single, "$"), dict_yes$igSubcloneID) | grepl(paste0("-", sbc_single, "-"), dict_yes$igSubcloneID),]

        if(length(unique(resc_df$igClonotypeID_num)) > 1 | nrow(resc_df) == 0){ next }

        selection <- meta$igSubcloneID == sbc_single & meta$completeBCR %in% c("Single_chain_1", "Single_chain_2")
        meta$igClonotypeID_num[selection] <- unique(resc_df$igClonotypeID_num)
        meta$igClonotypeVariantID_num[selection] <- names(which.min(sapply(unique(resc_df$igClonotypeVariantID_num), function(x) as.numeric(strsplit(x, "\\.CV")[[1]][2]))))
        meta$igSubcloneID_in_ClonotypeVariant_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_ClonotypeVariant_num), function(x) as.numeric(strsplit(x, "\\.S")[[1]][2]))))
        meta$igSubcloneID_in_Clonotype_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_Clonotype_num), function(x) as.numeric(strsplit(x, "\\.")[[1]][2]))))
        meta$completeBCR[selection] <- "Yes_rescue"
      }
    }
    meta$cell_id <- colnames(tmp)
    return(meta)
  }, mc.cores = threads)

  rescue_df_list <- Filter(Negate(is.null), rescue_df_list)
  rescue_df <- do.call(rbind, rescue_df_list)
  rownames(rescue_df) <- rescue_df$cell_id
  rescue_df <- rescue_df[,!colnames(rescue_df) %in% c("tmp_col", "cell_id")]

  colData(single_cell_object) <- DataFrame(rescue_df[colnames(single_cell_object),])

  return(single_cell_object)
}

.rescue_single_chain_seurat <- function(single_cell_object, group_col, threads, relaxed_rescue){

  if(!all(group_col %in% colnames(single_cell_object@meta.data))){stop(paste0("\nUnknown column (", group_col[!group_col %in% colnames(single_cell_object@meta.data)], ") selected for BCR ID recalculation! Please, set a valid column name."))}
  single_cell_object@meta.data$tmp_col <- apply(single_cell_object@meta.data[,group_col, drop = FALSE], 1, function(row) paste(row, collapse = "_"))

  rescue_df_list <- mclapply(unique(single_cell_object@meta.data$tmp_col), function(col_v){

    tmp <- subset(single_cell_object, subset = tmp_col == col_v)

    if(nrow(tmp@meta.data) == 0){return(NULL)}

    cols_dict <- c("completeBCR", "igClonotypeID", "igClonotypeID_num", "igClonotypeVariantID",
                   "igClonotypeVariantID_num", "igSubcloneID", "igSubcloneID_in_ClonotypeVariant_num",
                   "igSubcloneID_in_Clonotype_num")

    dict_yes <- tmp@meta.data[tmp@meta.data$completeBCR %in% "Yes", cols_dict, drop = FALSE]
    dict_yes <- dict_yes[!duplicated(dict_yes),]

    dict_single <- tmp@meta.data[tmp@meta.data$completeBCR %in% c("Single_chain_1", "Single_chain_2"), cols_dict, drop = FALSE]
    dict_single <- dict_single[!duplicated(dict_single),]

    if(relaxed_rescue){
      for(row in 1:nrow(dict_single)){
        clt_single <- dict_single$igClonotypeID[row]
        resc_df <- dict_yes[grepl(paste0("^", clt_single, "-"), dict_yes$igClonotypeID) | grepl(paste0("-", clt_single, "$"), dict_yes$igClonotypeID) | grepl(paste0("-", clt_single, "-"), dict_yes$igClonotypeID),]

        if(length(unique(resc_df$igClonotypeID_num)) > 1 | nrow(resc_df) == 0){next}

        selection <- tmp$igClonotypeID == clt_single & tmp$completeBCR %in% c("Single_chain_1", "Single_chain_2")
        tmp$igClonotypeID_num[selection] <- unique(resc_df$igClonotypeID_num)
        tmp$igClonotypeVariantID_num[selection] <- names(which.min(sapply(unique(resc_df$igClonotypeVariantID_num), function(x) as.numeric(strsplit(x, "\\.CV")[[1]][2]))))
        tmp$igSubcloneID_in_ClonotypeVariant_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_ClonotypeVariant_num), function(x) as.numeric(strsplit(x, "\\.S")[[1]][2]))))
        tmp$igSubcloneID_in_Clonotype_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_Clonotype_num), function(x) as.numeric(strsplit(x, "\\.")[[1]][2]))))
        tmp$completeBCR[selection] <- "Yes_rescue"
      }
    } else{
      for(row in 1:nrow(dict_single)){
        sbc_single <- dict_single$igSubcloneID[row]
        resc_df <- dict_yes[grepl(paste0("^", sbc_single, "-"), dict_yes$igSubcloneID) | grepl(paste0("-", sbc_single, "$"), dict_yes$igSubcloneID) | grepl(paste0("-", sbc_single, "-"), dict_yes$igSubcloneID),]

        if(length(unique(resc_df$igClonotypeID_num)) > 1 | nrow(resc_df) == 0){next}

        selection <- tmp$igSubcloneID == sbc_single & tmp$completeBCR %in% c("Single_chain_1", "Single_chain_2")
        tmp$igClonotypeID_num[selection] <- unique(resc_df$igClonotypeID_num)
        tmp$igClonotypeVariantID_num[selection] <- names(which.min(sapply(unique(resc_df$igClonotypeVariantID_num), function(x) as.numeric(strsplit(x, "\\.CV")[[1]][2]))))
        tmp$igSubcloneID_in_ClonotypeVariant_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_ClonotypeVariant_num), function(x) as.numeric(strsplit(x, "\\.S")[[1]][2]))))
        tmp$igSubcloneID_in_Clonotype_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_Clonotype_num), function(x) as.numeric(strsplit(x, "\\.")[[1]][2]))))
        tmp$completeBCR[selection] <- "Yes_rescue"
      }
    }
    return(tmp@meta.data)
  }, mc.cores = threads)

  rescue_df_list <- Filter(Negate(is.null), rescue_df_list)
  rescue_df <- do.call(rbind, rescue_df_list)
  rescue_df <- rescue_df[, colnames(rescue_df) != "tmp_col"]

  single_cell_object@meta.data <- rescue_df[Cells(single_cell_object),]

  return(single_cell_object)
}

.rescue_single_chain_igscandf <- function(igscan_df, group_col, threads, relaxed_rescue){

  if(!all(group_col %in% colnames(igscan_df))){stop(paste0("\nUnknown column (", group_col[!group_col %in% colnames(igscan_df@meta.data)], ") selected for BCR ID recalculation! Please, set a valid column name."))}
  igscan_df$tmp_col <- apply(igscan_df[,group_col, drop = FALSE], 1, function(row) paste(row, collapse = "_"))

  rescue_df_list <- mclapply(unique(igscan_df$tmp_col), function(col_v){

    tmp <- subset(igscan_df, subset = tmp_col == col_v)

    if(nrow(tmp) == 0){return(NULL)}

    cols_dict <- c("completeBCR", "igClonotypeID", "igClonotypeID_num", "igClonotypeVariantID",
                   "igClonotypeVariantID_num", "igSubcloneID", "igSubcloneID_in_ClonotypeVariant_num",
                   "igSubcloneID_in_Clonotype_num")

    dict_yes <- tmp[tmp$completeBCR %in% "Yes", cols_dict, drop = FALSE]
    dict_yes <- dict_yes[!duplicated(dict_yes),]

    dict_single <- tmp[tmp$completeBCR %in% c("Single_chain_1", "Single_chain_2"), cols_dict, drop = FALSE]
    dict_single <- dict_single[!duplicated(dict_single),]

    if(relaxed_rescue){
      for(row in 1:nrow(dict_single)){
        clt_single <- dict_single$igClonotypeID[row]
        resc_df <- dict_yes[grepl(paste0("^", clt_single, "-"), dict_yes$igClonotypeID) | grepl(paste0("-", clt_single, "$"), dict_yes$igClonotypeID) | grepl(paste0("-", clt_single, "-"), dict_yes$igClonotypeID),]

        if(length(unique(resc_df$igClonotypeID_num)) > 1 | nrow(resc_df) == 0){next}

        selection <- tmp$igClonotypeID == clt_single & tmp$completeBCR %in% c("Single_chain_1", "Single_chain_2")
        tmp$igClonotypeID_num[selection] <- unique(resc_df$igClonotypeID_num)
        tmp$igClonotypeVariantID_num[selection] <- names(which.min(sapply(unique(resc_df$igClonotypeVariantID_num), function(x) as.numeric(strsplit(x, "\\.CV")[[1]][2]))))
        tmp$igSubcloneID_in_ClonotypeVariant_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_ClonotypeVariant_num), function(x) as.numeric(strsplit(x, "\\.S")[[1]][2]))))
        tmp$igSubcloneID_in_Clonotype_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_Clonotype_num), function(x) as.numeric(strsplit(x, "\\.")[[1]][2]))))
        tmp$completeBCR[selection] <- "Yes_rescue"
      }
    } else{
      for(row in 1:nrow(dict_single)){
        sbc_single <- dict_single$igSubcloneID[row]
        resc_df <- dict_yes[grepl(paste0("^", sbc_single, "-"), dict_yes$igSubcloneID) | grepl(paste0("-", sbc_single, "$"), dict_yes$igSubcloneID) | grepl(paste0("-", sbc_single, "-"), dict_yes$igSubcloneID),]

        if(length(unique(resc_df$igClonotypeID_num)) > 1 | nrow(resc_df) == 0){next}

        selection <- tmp$igSubcloneID == sbc_single & tmp$completeBCR %in% c("Single_chain_1", "Single_chain_2")
        tmp$igClonotypeID_num[selection] <- unique(resc_df$igClonotypeID_num)
        tmp$igClonotypeVariantID_num[selection] <- names(which.min(sapply(unique(resc_df$igClonotypeVariantID_num), function(x) as.numeric(strsplit(x, "\\.CV")[[1]][2]))))
        tmp$igSubcloneID_in_ClonotypeVariant_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_ClonotypeVariant_num), function(x) as.numeric(strsplit(x, "\\.S")[[1]][2]))))
        tmp$igSubcloneID_in_Clonotype_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_Clonotype_num), function(x) as.numeric(strsplit(x, "\\.")[[1]][2]))))
        tmp$completeBCR[selection] <- "Yes_rescue"
      }
    }
    return(tmp)
  }, mc.cores = threads)

  rescue_df_list <- Filter(Negate(is.null), rescue_df_list)
  rescue_df <- do.call(rbind, rescue_df_list)
  rescue_df <- rescue_df[, colnames(rescue_df) != "tmp_col"]

  return(rescue_df)
}
