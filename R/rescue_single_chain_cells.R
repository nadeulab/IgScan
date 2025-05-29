#' Rescue cells with incomplete BCR (single chains)
#'
#' This function identifies and rescues cells with incomplete BCR information (labeled
#' as "Single_chain_1" or "Single_chain_2") by attempting to infer and assign clonotype
#' and subclonotype IDs based on related cells that have complete BCR information. It
#' supports both `SingleCellExperiment` and `Seurat` objects. After rescuing single-chain
#' entries, the function also recalculates clonotype-related identifiers.
#'
#' @param single_cell_object A single-cell object combined with IgScan annotation
#' in either `SingleCellExperiment` or `Seurat` format.
#' @param group_col A vector with the name of the column or columns containing the grouping variable/s.
#' Default is 'orig.ident', thus rescuing and recalculating BCR IDs based in each sample independently.
#' @param threads The number of threads to perform single_chain cells rescuing. Default is 1.
#'
#' @return A single cell object with updated metadata, where the `completeBCR`field in
#' some "single chain" cells may now be marked as "Yes_rescue".
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
#' }
#'
rescue_single_chain_cells <- function(single_cell_object, group_col = "orig.ident", threads = 1){

  if(class(single_cell_object)[1] == "SingleCellExperiment"){
    single_cell_object <- .rescue_single_chain_sce(single_cell_object, group_col, threads)
  } else if(class(single_cell_object)[1] == "Seurat"){
    single_cell_object <- .rescue_single_chain_seurat(single_cell_object, group_col, threads)
  }

  single_cell_object <- recalculate_IDs_single_cell(single_cell_object, group_col, threads)

  return(single_cell_object)
}

.rescue_single_chain_sce <- function(single_cell_object, group_col, threads){

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

    for(row in 1:nrow(dict_single)){
      sbc_single <- dict_single$igSubcloneID[row]
      resc_df <- dict_yes[grepl(paste0(sbc_single, "-"), dict_yes$igSubcloneID) | grepl(paste0("-", sbc_single), dict_yes$igSubcloneID) | grepl(paste0("-", sbc_single, "-"), dict_yes$igSubcloneID),]

      if(length(unique(resc_df$igClonotypeID_num)) > 1 | nrow(resc_df) == 0){ next }

      selection <- meta$igSubcloneID == sbc_single & meta$completeBCR %in% c("Single_chain_1", "Single_chain_2")
      meta$igClonotypeID_num[selection] <- unique(resc_df$igClonotypeID_num)
      meta$igClonotypeVariantID_num[selection] <- names(which.min(sapply(unique(resc_df$igClonotypeVariantID_num), function(x) as.numeric(strsplit(x, "\\.CV")[[1]][2]))))
      meta$igSubcloneID_in_ClonotypeVariant_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_ClonotypeVariant_num), function(x) as.numeric(strsplit(x, "\\.S")[[1]][2]))))
      meta$igSubcloneID_in_Clonotype_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_Clonotype_num), function(x) as.numeric(strsplit(x, "\\.")[[1]][2]))))
      meta$completeBCR[selection] <- "Yes_rescue"
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

.rescue_single_chain_seurat <- function(single_cell_object, group_col, threads){

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

    for(row in 1:nrow(dict_single)){
      sbc_single <- dict_single$igSubcloneID[row]
      resc_df <- dict_yes[grepl(paste0(sbc_single, "-"), dict_yes$igSubcloneID) | grepl(paste0("-", sbc_single), dict_yes$igSubcloneID) | grepl(paste0("-", sbc_single, "-"), dict_yes$igSubcloneID),]

      if(length(unique(resc_df$igClonotypeID_num)) > 1 | nrow(resc_df) == 0){next}

      selection <- tmp@meta.data$igSubcloneID == sbc_single & tmp@meta.data$completeBCR %in% c("Single_chain_1", "Single_chain_2")
      tmp@meta.data$igClonotypeID_num[selection] <- unique(resc_df$igClonotypeID_num)
      tmp@meta.data$igClonotypeVariantID_num[selection] <- names(which.min(sapply(unique(resc_df$igClonotypeVariantID_num), function(x) as.numeric(strsplit(x, "\\.CV")[[1]][2]))))
      tmp@meta.data$igSubcloneID_in_ClonotypeVariant_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_ClonotypeVariant_num), function(x) as.numeric(strsplit(x, "\\.S")[[1]][2]))))
      tmp@meta.data$igSubcloneID_in_Clonotype_num[selection] <- names(which.min(sapply(unique(resc_df$igSubcloneID_in_Clonotype_num), function(x) as.numeric(strsplit(x, "\\.")[[1]][2]))))
      tmp@meta.data$completeBCR[selection] <- "Yes_rescue"
    }
    return(tmp@meta.data)
  }, mc.cores = threads)

  rescue_df_list <- Filter(Negate(is.null), rescue_df_list)
  rescue_df <- do.call(rbind, rescue_df_list)
  rescue_df <- rescue_df[, colnames(rescue_df) != "tmp_col"]

  single_cell_object@meta.data <- rescue_df[Cells(single_cell_object),]

  return(single_cell_object)
}
