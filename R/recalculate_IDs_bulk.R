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
#'
#' @return An IgScan output data frame with updated clonotype-related IDs and
#' cumulative number of reads and frequencies.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' updated_igscan_out <- recalculate_IDs_bulk(igscan_output_df)
#' }

recalculate_IDs_bulk <- function(igscan_data_frame){

  ## Update Clonotype, ClonotypeVariant and Subclone Frequencies after removing
  igscan_data_frame <- igscan_data_frame %>%
    group_by(ClonotypeID) %>%
    mutate(Clonotype_nReads = sum(Subclone_nReads))

  igscan_data_frame$Clonotype_freq <- (igscan_data_frame$Clonotype_nReads/sum(igscan_data_frame$Clonotype_nReads[!duplicated(igscan_data_frame$ClonotypeID)]))*100

  igscan_data_frame <- igscan_data_frame %>%
    group_by(ClonotypeID) %>%
    mutate(Subclone_freq_in_Clonotype = (Subclone_nReads/sum(Subclone_nReads))*100,
           ClonotypeVariant_freq = (clonotypeVariant_nReads / sum(clonotypeVariant_nReads[!duplicated(ClonotypeVariantID)]))*100)

  igscan_data_frame <- igscan_data_frame %>%
    group_by(ClonotypeVariantID) %>%
    mutate(Subclone_freq_in_ClonotypeVariant = (Subclone_nReads/sum(Subclone_nReads))*100)

  # Update Clonotype, ClonotypeVariant and Subclone IDs
  clt_dict <- igscan_data_frame[!duplicated(igscan_data_frame$ClonotypeID), c("ClonotypeID","Clonotype_freq")]
  clt_dict <- clt_dict[order(clt_dict$Clonotype_freq, decreasing = T),]
  clt_dict$NewClonotypeID <- paste0("C",1:nrow(clt_dict))
  igscan_data_frame$ClonotypeID <- clt_dict$NewClonotypeID[match(igscan_data_frame$ClonotypeID, clt_dict$ClonotypeID)]

  ## ClonotypeVariantID
  cv_dict_all <- data.frame()
  sc_in_clt_dict_all <- data.frame()
  for(c in unique(igscan_data_frame$ClonotypeID)){
    cv_dict <- igscan_data_frame[igscan_data_frame$ClonotypeID == c, c("ClonotypeVariantID","ClonotypeVariant_freq")]
    cv_dict <- cv_dict[!duplicated(cv_dict$ClonotypeVariantID),]
    cv_dict <- cv_dict[order(cv_dict$ClonotypeVariant_freq, decreasing = T),]
    cv_dict$New_cvID <- paste0(c,".CV",1:nrow(cv_dict))
    cv_dict_all <- rbind(cv_dict_all, cv_dict)

    sc_in_clt_dict <- igscan_data_frame[igscan_data_frame$ClonotypeID == c, c("SubcloneID_in_Clonotype","Subclone_freq_in_Clonotype")]
    sc_in_clt_dict <- sc_in_clt_dict[order(sc_in_clt_dict$Subclone_freq_in_Clonotype, decreasing = T),]
    sc_in_clt_dict$New_scID <- paste0(c,".",1:nrow(sc_in_clt_dict))
    sc_in_clt_dict_all <- rbind(sc_in_clt_dict_all, sc_in_clt_dict)
  }
  igscan_data_frame$ClonotypeVariantID <- cv_dict_all$New_cvID[match(igscan_data_frame$ClonotypeVariantID, cv_dict_all$ClonotypeVariantID)]
  igscan_data_frame$SubcloneID_in_Clonotype <- sc_in_clt_dict_all$New_scID[match(igscan_data_frame$SubcloneID_in_Clonotype, sc_in_clt_dict_all$SubcloneID_in_Clonotype)]

  ## SubcloneID_in_CV
  sc_dict_all <- data.frame()
  for(cv in unique(igscan_data_frame$ClonotypeVariantID)){
    sc_dict <- igscan_data_frame[igscan_data_frame$ClonotypeVariantID == cv, c("SubcloneID","Subclone_freq_in_ClonotypeVariant")]
    sc_dict <- sc_dict[order(sc_dict$Subclone_freq_in_ClonotypeVariant, decreasing = T),]
    sc_dict$New_scID <- paste0(cv,".S",1:nrow(sc_dict))
    sc_dict_all <- rbind(sc_dict_all, sc_dict)
  }
  igscan_data_frame$SubcloneID <- sc_dict_all$New_scID[match(igscan_data_frame$SubcloneID, sc_dict_all$SubcloneID)]

  return(as.data.frame(igscan_data_frame))
}
