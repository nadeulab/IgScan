#' Filter BCR contamination in bulk IG-NGS data
#'
#' This function identifies, flags, and optionally removes potential contamination
#' between samples processed and/or sequenced together in a IgScan-annotated bulk
#' NGS dataset. Contamination is assessed by analyzing the frequency of identical
#' IG sequences identified in different samples A contamination cutoff can be specified,
#' indicating that if a clonotype is observed X times more in one sample than another,
#' it will be flagged as contamination in the latter sample.
#'
#' @param igscan_data_frame An IgScan output dataframe coming from the annotation
#' of bulk NGS data.
#' @param batch_col A string with the name of the column identifying the `batch` of processing
#' and/or sequencing for each sample. Each `batch` will be evaluated independently. Default is NULL,
#' meaning that all sequences in the data frame will be considered from the same batch.
#' @param case_col A string with the name of the column identifying the `case` or `individual` for
#' each sample. Default is NULL, meaning that each sample will be considered from different individuals.
#' @param contamination_cutoff A numeric value specifying the contamination ratio threshold.
#' Default is 10, higher values result in stricter contamination detection.
#' @param contamination_clone_cutoff A numeric value defining what percentage of contaminated sequences within
#' a clonotype is allowed (i.e., sequences sharing the same clonotype ID). If the proportion of contaminated
#' sequences exceeds this threshold, all sequences beloning to the clonptype will be flagged as `OUT_CLONOTYPE`.
#' Default is 80 percent.
#' @param remove_contamination Logical. If `TRUE`, sequences flagged as contaminated are
#' removed from the returned data frame (default is `FALSE`).
#' @param recalc_column A string/vector of strings with the name of the column/s to use for ID
#' recalculation. Needed when remove_contamination=TRUE. Default is `SampleID`.
#' @param threads The number of threads to perform BCR ID recalculation. Default is 1.
#'
#' @return The input `igscan_data_frame` including three new columns:
#' \itemize{
#'   \item Contamination_FLAG: Indicates whether a cell passed contamination checks (`PASS`), was flagged as contaminated (`OUT`),
#'   or was involved in crossed contamination (`OUT_X`)
#'   \item Contamination_Freq: The frequency of the contaminated sequence in the original sample (dominant).
#'   \item Contamination_Sample: The dominant sample responsible for contamination.
#' }
#'
#'   If `remove_contamination = TRUE`, sequences flagged as `OUT` or `OUT_X` are removed from
#'   the data frame and immunogenetic IDs are recalculated with `recalculate_IDs_bulk()`.
#'
#' @export
#'
#' @importFrom parallel mclapply
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' cont_filtered_igscan <- filter_BCR_contam_bulk(igscan_data_frame = igscan_df, bath_col = "batchID", case_col = "case", contamination_cutoff = 30, threads = 4)
#' }
#'
filter_BCR_contam_bulk <- function(igscan_data_frame, batch_col = NULL, sample_col = "SampleID", case_col = NULL, contamination_cutoff = 10, contamination_clone_cutoff = 80, remove_contamination = F, recalc_column = "SampleID", threads = 1){

  igscan_data_frame <- as.data.frame(igscan_data_frame)

  if(!is.null(batch_col) && !batch_col %in% colnames(igscan_data_frame)){stop(paste0("\nUnknown batch column (", batch_col, ") selected for flagging BCR contamination! Please, set a valid column name."))}
  if(!is.null(case_col) && !case_col %in% colnames(igscan_data_frame)){stop(paste0("\nUnknown case column (", case_col, ") selected for flagging BCR contamination! Please, set a valid column name."))}
  if(!all(recalc_column %in% colnames(igscan_data_frame)) & remove_contamination){stop(paste0("\nUnknown recalc_column column (", recalc_column[!recalc_column %in% colnames(igscan_data_frame)], ") selected for flagging BCR contamination! Please, set a valid column name."))}

  if(is.null(batch_col)){
    igscan_data_frame$BatchID <- "BatchX"
  } else{
    colnames(igscan_data_frame)[colnames(igscan_data_frame) == batch_col] <- "BatchID"
  }

  if(is.null(case_col)){
    igscan_data_frame$CaseID <- paste0("Case_", igscan_data_frame$SampleID)
  } else{
    colnames(igscan_data_frame)[colnames(igscan_data_frame) == case_col] <- "CaseID"
  }

  contam_df_list <- mclapply(unique(igscan_data_frame$BatchID), function(col_v){

    samples_in_batch <- unique(igscan_data_frame$SampleID[igscan_data_frame$BatchID == col_v])
    all_subclones_in_batch <- igscan_data_frame[igscan_data_frame$SampleID %in% samples_in_batch,]

    ## Calculate total ClonotypeVariant frequency
    all_subclones_in_batch <- all_subclones_in_batch %>%
      group_by(SampleID) %>%
      mutate(ClonotypeVariant_freqTotal = ClonotypeVariant_nReads/sum(ClonotypeVariant_nReads[!duplicated(ClonotypeVariantID)])*100)

    all_subclones_in_batch$SubcloneBarcode <- paste0(all_subclones_in_batch$VDJ_sequence, "-", all_subclones_in_batch$InDels)

    all_subclones_in_batch$Contamination_FLAG <- "PASS"
    all_subclones_in_batch$Contamination_Freq <- NA
    all_subclones_in_batch$Contamination_Sample <- NA

    mm <- all_subclones_in_batch[, c("CaseID", "SubcloneBarcode")]
    mm <- mm[!duplicated(mm),]
    contamination_sequences <- unique(mm$SubcloneBarcode[duplicated(mm$SubcloneBarcode)])

    message(paste0("There is a total of ", length(contamination_sequences), " sequences that could be contamination in batch ", col_v, ". Starting to process...\n"))

    for(cont_seq in contamination_sequences){

      cont_seq_df <- all_subclones_in_batch[all_subclones_in_batch$SubcloneBarcode == cont_seq,]

      maxFreq <- max(cont_seq_df$ClonotypeVariant_freqTotal)
      maxSample <- cont_seq_df$SampleID[cont_seq_df$ClonotypeVariant_freqTotal == maxFreq]
      maxCase <- cont_seq_df$CaseID[cont_seq_df$ClonotypeVariant_freqTotal == maxFreq]

      cont_seq_df <- cont_seq_df[cont_seq_df$CaseID != maxCase,]

      v <- maxFreq/cont_seq_df$ClonotypeVariant_freqTotal >= contamination_cutoff
      if(any(v)){
        all_subclones_in_batch$Contamination_FLAG[all_subclones_in_batch$SampleID %in% cont_seq_df$SampleID[v] & all_subclones_in_batch$SubcloneBarcode == cont_seq] <- "OUT"
        all_subclones_in_batch$Contamination_Freq[all_subclones_in_batch$SampleID %in% cont_seq_df$SampleID[v] & all_subclones_in_batch$SubcloneBarcode == cont_seq] <- maxFreq
        all_subclones_in_batch$Contamination_Sample[all_subclones_in_batch$SampleID %in% cont_seq_df$SampleID[v] & all_subclones_in_batch$SubcloneBarcode == cont_seq] <- maxSample
      }

      v <- maxFreq/cont_seq_df$ClonotypeVariant_freqTotal < contamination_cutoff
      if(any(v)){
        all_subclones_in_batch$Contamination_FLAG[all_subclones_in_batch$SampleID %in% cont_seq_df$SampleID[v] & all_subclones_in_batch$SubcloneBarcode == cont_seq] <- "OUT_X"
        all_subclones_in_batch$Contamination_Freq[all_subclones_in_batch$SampleID %in% cont_seq_df$SampleID[v] & all_subclones_in_batch$SubcloneBarcode == cont_seq] <- maxFreq
        all_subclones_in_batch$Contamination_Sample[all_subclones_in_batch$SampleID %in% cont_seq_df$SampleID[v] & all_subclones_in_batch$SubcloneBarcode == cont_seq] <- maxSample

        all_subclones_in_batch$Contamination_FLAG[all_subclones_in_batch$SampleID == maxSample & all_subclones_in_batch$SubcloneBarcode == cont_seq] <- "OUT_X"
        all_subclones_in_batch$Contamination_Freq[all_subclones_in_batch$SampleID == maxSample & all_subclones_in_batch$SubcloneBarcode == cont_seq] <- max(cont_seq_df$ClonotypeVariant_freqTotal[v])
        all_subclones_in_batch$Contamination_Sample[all_subclones_in_batch$SampleID == maxSample & all_subclones_in_batch$SubcloneBarcode == cont_seq] <- cont_seq_df$SampleID[cont_seq_df$ClonotypeVariant_freqTotal == max(cont_seq_df$ClonotypeVariant_freqTotal[v])]
      }
    }

    batch_annot_cont_df <- data.frame()
    ## Once contamination has been determined, we apply it
    for(sample in unique(all_subclones_in_batch$SampleID)){

      sample_df <- all_subclones_in_batch[all_subclones_in_batch$SampleID == sample, !colnames(all_subclones_in_batch) %in% c("SubcloneBarcode", "ClonotypeVariant_freqTotal")]
      cont_clonotypeVariants <- unique(sample_df$ClonotypeVariantID[sample_df$Contamination_FLAG != "PASS"])
      potential_cont_clonotypes <- unique(sapply(cont_clonotypeVariants, function(x) paste(strsplit(x, "\\.")[[1]][1:2], collapse = ".")))

      for(clon in potential_cont_clonotypes){
        clonVars <- cont_clonotypeVariants[startsWith(cont_clonotypeVariants, paste0(clon, "."))]
        clonVarsFreqInClon <- sum(sample_df$ClonotypeVariant_freq[!duplicated(sample_df$ClonotypeVariantID) & sample_df$ClonotypeVariantID %in% clonVars])
        if(clonVarsFreqInClon >= contamination_clone_cutoff){
          sample_df$Contamination_FLAG[sample_df$Contamination_FLAG != "PASS" & sample_df$ClonotypeID == clon] <- "OUT_CLONOTYPE"
        }else{
          sample_df$Contamination_FLAG[sample_df$Contamination_FLAG != "PASS" & sample_df$ClonotypeVariantID %in% clonVars] <- "OUT_CLONOTYPEVARIANT"
        }
      }
      batch_annot_cont_df <- rbind(batch_annot_cont_df, sample_df)
    }
    return(as.data.frame(batch_annot_cont_df))
  }, mc.cores = threads)

  contam_df <- do.call(rbind, contam_df_list)

  if(remove_contamination){
    contam_df <- contam_df[contam_df$Contamination_FLAG == "PASS",]
    contam_df <- recalculate_IDs_bulk(contam_df, group_col = recalc_column)
  }
  return(contam_df)
}

