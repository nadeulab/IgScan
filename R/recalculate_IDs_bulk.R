#' Recalculate IgScan IDs for Bulk NGS Data
#'
#' @description This function recalculates and standardizes clonotype-related IDs
#' in IgScan outputs of bulk NGS datasets. Using this function may be of interest
#' when a IgScan-annotated dataframe has been filtered and some sequences
#' have been removed. It also updates the cumulative number of reads and frequencies
#' of the different IDs.
#'
#' @param igscan_data_frame An IgScan output dataframe coming from the annotation
#' of bulk NGS data.
#' @param group_col A vector with the name of the column or columns containing the grouping variable/s.
#' Default is 'SampleID', thus recalculating BCR IDs based in each sample independently.
#' @param threads The number of threads to perform BCR ID recalculation. Default is 1.
#'
#' @return An IgScan output data frame with updated clonotype-related IDs and
#' cumulative number of reads and frequencies.
#'
#' @export
#'
#' @import dplyr
#' @importFrom parallel mclapply
#'
#' @examples
#' \dontrun{
#' updated_igscan_out <- recalculate_IDs_bulk(igscan_output_df group_col = "SampleID", threads = 4)
#' }

recalculate_IDs_bulk <- function(igscan_data_frame, group_col = "SampleID", threads = 1){

  igscan_data_frame <- as.data.frame(igscan_data_frame)

  if(!all(group_col %in% colnames(igscan_data_frame))){stop(paste0("Unknown column (", group_col, ") selected for BCR ID recalculation! Please, set a valid column name."))}

  igscan_data_frame$tmp_col <- apply(igscan_data_frame[,group_col, drop = FALSE], 1, function(row) paste(row, collapse = "_"))

  recalc_df_list <- mclapply(unique(igscan_data_frame$tmp_col), function(col_v){

    tmp_df <- igscan_data_frame[igscan_data_frame$tmp_col == col_v,]

    ## Update Clonotype, ClonotypeVariant and Subclone Frequencies after removing
    tmp_df <- tmp_df %>%
      group_by(ClonotypeID) %>%
      mutate(Clonotype_nReads = sum(Subclone_nReads))

    tmp_df$Clonotype_freq <- (tmp_df$Clonotype_nReads/sum(tmp_df$Clonotype_nReads[!duplicated(tmp_df$ClonotypeID)]))*100

    tmp_df <- tmp_df %>%
      group_by(ClonotypeID) %>%
      mutate(Subclone_freq_in_Clonotype = (Subclone_nReads/sum(Subclone_nReads))*100,
             ClonotypeVariant_freq = (ClonotypeVariant_nReads / sum(ClonotypeVariant_nReads[!duplicated(ClonotypeVariantID)]))*100)

    tmp_df <- tmp_df %>%
      group_by(ClonotypeVariantID) %>%
      mutate(Subclone_freq_in_ClonotypeVariant = (Subclone_nReads/sum(Subclone_nReads))*100)

    # Update Clonotype, ClonotypeVariant and Subclone IDs
    clt_dict <- tmp_df[!duplicated(tmp_df$ClonotypeID), c("ClonotypeID","Clonotype_freq")]
    clt_dict <- clt_dict[order(clt_dict$Clonotype_freq, decreasing = T),]
    clt_dict$NewClonotypeID <- paste0("C",1:nrow(clt_dict))
    tmp_df$ClonotypeID <- clt_dict$NewClonotypeID[match(tmp_df$ClonotypeID, clt_dict$ClonotypeID)]

    ## ClonotypeVariantID
    cv_dict_all <- data.frame()
    sc_in_clt_dict_all <- data.frame()
    for(c in unique(tmp_df$ClonotypeID)){
      cv_dict <- tmp_df[tmp_df$ClonotypeID == c, c("ClonotypeVariantID","ClonotypeVariant_freq")]
      cv_dict <- cv_dict[!duplicated(cv_dict$ClonotypeVariantID),]
      cv_dict <- cv_dict[order(cv_dict$ClonotypeVariant_freq, decreasing = T),]
      cv_dict$New_cvID <- paste0(c,".CV",1:nrow(cv_dict))
      cv_dict_all <- rbind(cv_dict_all, cv_dict)

      sc_in_clt_dict <- tmp_df[tmp_df$ClonotypeID == c, c("SubcloneID_in_Clonotype","Subclone_freq_in_Clonotype")]
      sc_in_clt_dict <- sc_in_clt_dict[order(sc_in_clt_dict$Subclone_freq_in_Clonotype, decreasing = T),]
      sc_in_clt_dict$New_scID <- paste0(c,".",1:nrow(sc_in_clt_dict))
      sc_in_clt_dict_all <- rbind(sc_in_clt_dict_all, sc_in_clt_dict)
    }
    tmp_df$ClonotypeVariantID <- cv_dict_all$New_cvID[match(tmp_df$ClonotypeVariantID, cv_dict_all$ClonotypeVariantID)]
    tmp_df$SubcloneID_in_Clonotype <- sc_in_clt_dict_all$New_scID[match(tmp_df$SubcloneID_in_Clonotype, sc_in_clt_dict_all$SubcloneID_in_Clonotype)]

    ## SubcloneID_in_CV
    sc_dict_all <- data.frame()
    for(cv in unique(tmp_df$ClonotypeVariantID)){
      sc_dict <- tmp_df[tmp_df$ClonotypeVariantID == cv, c("SubcloneID","Subclone_freq_in_ClonotypeVariant")]
      sc_dict <- sc_dict[order(sc_dict$Subclone_freq_in_ClonotypeVariant, decreasing = T),]
      sc_dict$New_scID <- paste0(cv,".S",1:nrow(sc_dict))
      sc_dict_all <- rbind(sc_dict_all, sc_dict)
    }
    tmp_df$SubcloneID <- sc_dict_all$New_scID[match(tmp_df$SubcloneID, sc_dict_all$SubcloneID)]
    return(tmp_df)
  }, mc.cores = threads)

  recalc_df <- do.call(rbind, recalc_df_list)
  recalc_df <- recalc_df[, colnames(recalc_df) != "tmp_col"]

  return(as.data.frame(recalc_df))
}
