#' Filter BCR-based contamination in single-cell data
#'
#' This function identifies, flags, and optionally removes potential contamination
#' between samples processed and/or sequenced together fors cRNA-seq+scBCR-seq.
#' Contamination is assessed by analyzing the distribution of identicalBCR clonotypes
#' between samples. A contamination cutoff can be specified, indicating that if a
#' clonotype is observed X times more in one sample tag than another, it will be
#' flagged as contamination in the latter sample. The function works with both
#' `SingleCellExperiment` and `Seurat` objects.
#'
#' @param single_cell_list A list of `SingleCellExperiment` or `Seurat` objects containing IgScan
#' annotation for single cell data.
#' @param sample_col A string with the name of the column identifying `sample`.
#' Default is 'orig.ident'.
#' @param batch_col A string with the name of the column identifying the `batch` of processing
#' and/or sequencing for each sample. Each `batch` will be evaluated independently. Default is NULL,
#' meaning that all samples in the single_cell_list list will be considered from the same batch.
#' @param case_col A string with the name of the column identifying the `case` or `individual` for
#' each sample. Default is NULL, meaning that each sample will be considered from different individuals.
#' @param contamination_cutoff A numeric value specifying the contamination ratio threshold.
#' Default is 10, higher values result in stricter contamination detection.
#' @param contamination_clone_cutoff A numeric value defining what percentage of contaminated cells within
#' a clonotype is allowed (i.e., cells sharing the same clonotype ID). If the proportion of contaminated
#' cells exceeds this threshold, all cells beloning to the clonptype will be flagged as `OUT_CLONOTYPE`.
#' Default is 80 percent.
#' @param remove_contamination Logical. If `TRUE`, cells flagged as contaminated are
#' removed from the returned object (default is `FALSE`).
#' @param recalc_column A string/vector of strings with the name of the column/s to use for ID
#' recalculation. Needed when remove_contamination=TRUE. Default is `orig.ident`.
#' @param threads The number of threads to perform BCR ID recalculation. Default is 1.
#'
#' @return The input `single_cell_list` with each object including three new columns:
#' \itemize{
#'   \item Contamination_FLAG: Indicates whether a cell passed contamination checks (`PASS`), was flagged as contaminated (`OUT`),
#'   or was involved in crossed contamination (`OUT_X`)
#'   \item Contamination_Ratio: The ratio of clonal overlap between the contaminated and dominant samples.
#'   \item Contamination_Sample: The dominant sample responsible for contamination.
#' }
#'
#'   If `remove_contamination = TRUE`, cells flagged as `OUT` or `OUT_X` are removed from
#'   the objects and immunogenetic IDs are recalculated with `recalculate_IDs_single_cell()`.
#'
#' @export
#'
#' @importFrom parallel mclapply
#' @import dplyr
#' @import Seurat
#' @import SeuratObject
#' @import SingleCellExperiment
#'
#' @examples
#' \dontrun{
#' # Example with a SingleCellExperiment object:
#' cont_filtered_sce <- filter_BCR_contam_single_cell(single_cell_list = sc_list, bath_col = "batchID", case_col = "case", contamination_cutoff = 30, threads = 4)
#' }
#'
filter_BCR_contam_single_cell <- function(single_cell_list, sample_col = "orig.ident", batch_col = NULL, case_col = NULL, contamination_cutoff = 10, contamination_clone_cutoff = 80, remove_contamination = F, recalc_column = "orig.ident", threads = 1){

  meta_data <- data.frame()
  for(i in 1:length(single_cell_list)){
    if(class(single_cell_list[[i]])[1] == "SingleCellExperiment"){
      tmp <- colData(single_cell_list[[i]])
    } else if(class(single_cell_list[[i]])[1] == "Seurat"){
      tmp <- single_cell_list[[i]]@meta.data
    }
    meta_data <- rbind(meta_data, tmp)
  }

  if(!sample_col %in% colnames(meta_data)){stop(paste0("\nUnknown sample_col column (", sample_col, ") selected for flagging contamination! Please, set a valid column name."))}
  if(remove_contamination & !all(recalc_column %in% colnames(meta_data))){stop(paste0("\nUnknown recalc_column column (", recalc_column[!recalc_column %in% colnames(meta_data)], ") selected for flagging SampleTag contamination! Please, set a valid column name."))}

  if(is.null(batch_col)){
    meta_data$BatchID <- "BatchX"
  } else{
    if(batch_col %in% colnames(meta_data)){stop(paste0("\nUnknown batch_col column (", batch_col, ") selected for flagging contamination! Please, set a valid column name."))}
    colnames(meta_data)[colnames(meta_data) == batch_col] <- "BatchID"
  }

  if(is.null(case_col)){
    meta_data$CaseID <- paste0("Case_", meta_data[[sample_col]])
  } else{
    if(!case_col %in% colnames(meta_data)){stop(paste0("\nUnknown case_col column (", case_col, ") selected for flagging contamination! Please, set a valid column name."))}
    colnames(meta_data)[colnames(meta_data) == case_col] <- "CaseID"
  }

  contam_df_list <- mclapply(unique(meta_data$BatchID), function(col_v){

    all_subclones_in_batch <- data.frame()
    tmp_df <- meta_data[meta_data$BatchID == col_v,]
    samples_in_batch <- as.character(unique(tmp_df[[sample_col]]))

    samples_in_batch <- as.character(unique(meta_data[meta_data$BatchID == col_v, sample_col]))

    count_df <- tmp_df %>%
      group_by(!!sym(sample_col), igSubcloneID_all) %>%
      summarize(count = n(), .groups = "drop") %>%
      as.data.frame()

    tmp_df <- tmp_df %>%
      left_join(count_df, by = c(sample_col, "igSubcloneID_all")) %>%
      filter(!duplicated(paste(.[[sample_col]], igSubcloneID_all, sep = "_"))) %>%
      filter(!is.na(completeBCR)) %>%
      mutate(SubcloneBarcode = paste0(igVDJ_sequence, "-", igInDels),
             Contamination_FLAG = "PASS", Contamination_Ratio = NA, Contamination_Sample = NA)

    contamination_sequences <- unique(tmp_df$SubcloneBarcode[duplicated(tmp_df$SubcloneBarcode)])

    message(paste0("There is a total of ", length(contamination_sequences), " sequences that could be contamination in batch ", col_v, ". Starting to process..."))

    for(cont_seq in contamination_sequences){

      sub_df <- tmp_df[tmp_df$SubcloneBarcode == cont_seq,]
      max_nCell <- max(as.numeric(sub_df$count))
      if(max_nCell == 1){next}

      maxSample <- as.character(sub_df[sub_df$count == max_nCell, sample_col])[1]
      maxCase <- as.character(sub_df$CaseID[sub_df$count == max_nCell])[1]

      if(all(sub_df$CaseID == maxCase)) {next}
      sub_df <- sub_df[sub_df$CaseID != maxCase,]

      v <- (max_nCell/as.numeric(sub_df$count)) >= contamination_cutoff
      if(any(v)){
        tmp_df$Contamination_FLAG[tmp_df[[sample_col]] %in% as.character(sub_df[[sample_col]][v]) & tmp_df$SubcloneBarcode == cont_seq] <- "OUT"
        tmp_df$Contamination_Ratio[tmp_df[[sample_col]] %in% as.character(sub_df[[sample_col]][v]) & tmp_df$SubcloneBarcode == cont_seq] <- paste0(tmp_df$count[tmp_df[[sample_col]] %in% as.character(sub_df[[sample_col]][v]) & tmp_df$SubcloneBarcode == cont_seq],"|",max_nCell)
        tmp_df$Contamination_Sample[tmp_df[[sample_col]] %in% as.character(sub_df[[sample_col]][v]) & tmp_df$SubcloneBarcode == cont_seq] <- maxSample
      }

      v <- (max_nCell/as.numeric(sub_df$count)) < contamination_cutoff
      if(any(v)){
        tmp_df$Contamination_FLAG[tmp_df[[sample_col]] %in% as.character(sub_df[[sample_col]][v]) & tmp_df$SubcloneBarcode == cont_seq] <- "OUT_X"
        tmp_df$Contamination_Ratio[tmp_df[[sample_col]] %in% as.character(sub_df[[sample_col]][v]) & tmp_df$SubcloneBarcode == cont_seq] <- paste0(tmp_df$count[tmp_df[[sample_col]] %in% as.character(sub_df[[sample_col]][v]) & tmp_df$SubcloneBarcode == cont_seq],"|",max_nCell)
        tmp_df$Contamination_Sample[tmp_df[[sample_col]] %in% as.character(sub_df[[sample_col]][v]) & tmp_df$SubcloneBarcode == cont_seq] <- maxSample

        tmp_df$Contamination_FLAG[tmp_df$CaseID == maxCase & tmp_df$SubcloneBarcode == cont_seq] <- "OUT_X"
        tmp_df$Contamination_Ratio[tmp_df$CaseID == maxCase & tmp_df$SubcloneBarcode == cont_seq] <- paste0(max_nCell, "|", paste0(tmp_df$count[tmp_df[[sample_col]] %in% as.character(sub_df[[sample_col]][v]) & all_subclones_in_batch$SubcloneBarcode == cont_seq], collapse = "-"))
        tmp_df$Contamination_Sample[tmp_df$CaseID == maxCase & tmp_df$SubcloneBarcode == cont_seq] <- paste0(tmp_df[[sample_col]][tmp_df[[sample_col]] %in% as.character(sub_df[[sample_col]][v]) & all_subclones_in_batch$SubcloneBarcode == cont_seq], collapse = "-")
      }
    }

    if(sum(tmp_df$Contamination_FLAG != "PASS") == 0){return(NULL)}

    return(tmp_df[tmp_df$Contamination_FLAG != "PASS",])
  }, mc.cores = threads)

  contam_df_list <- Filter(Negate(is.null), contam_df_list)
  contam_df <- do.call(rbind, contam_df_list)

  ## Assign contamination back to Single-cell objects
  for(i in 1:length(single_cell_list)){
    if(class(single_cell_list[[i]])[1] == "SingleCellExperiment"){
      tmp_meta_data <- colData(single_cell_list[[i]])
      tmp_meta_data <- .flagCont_singlecell(tmp_meta_data, contam_df, sample_col, contamination_clone_cutoff)
      colData(single_cell_list[[i]]) <- DataFrame(tmp_meta_data)
      if(remove_contamination){
        single_cell_list[[i]] <- single_cell_list[[i]][,colData(single_cell_list[[i]])$Contamination_FLAG == "PASS"]
        single_cell_list[[i]] <- recalculate_IDs_single_cell(single_cell_list[[i]], recalc_column, threads)
      }

    } else if(class(single_cell_list[[i]])[1] == "Seurat"){
      tmp_meta_data <- single_cell_list[[i]]@meta.data
      tmp_meta_data <- .flagCont_singlecell(tmp_meta_data, contam_df, sample_col, contamination_clone_cutoff)
      single_cell_list[[i]]@meta.data <- tmp_meta_data
      if(remove_contamination){
        single_cell_list[[i]] <- subset(single_cell_list[[i]], subset = Contamination_FLAG == "PASS")
        single_cell_list[[i]] <- recalculate_IDs_single_cell(single_cell_list[[i]], recalc_column, threads)
      }
    }
  }
  return(single_cell_list)
}

.flagCont_singlecell <- function(object_meta_data, contam_df, sample_col, contamination_clone_cutoff){

  obj_id_comb <- paste0(object_meta_data[[sample_col]], "_", object_meta_data$igSubcloneID_all)
  contam_id_comb <- paste0(contam_df[[sample_col]], "_", contam_df$igSubcloneID_all)

  if(length(intersect(obj_id_comb, contam_id_comb)) == 0){
    object_meta_data$Contamination_FLAG <- "PASS"
    object_meta_data$Contamination_Ratio <- NA
    object_meta_data$Contamination_Sample <- NA

  } else{
    object_meta_data$Contamination_FLAG <- contam_df$Contamination_FLAG[match(obj_id_comb, contam_id_comb)]
    object_meta_data$Contamination_Ratio <- contam_df$Contamination_Ratio[match(obj_id_comb, contam_id_comb)]
    object_meta_data$Contamination_Sample <- contam_df$Contamination_Sample[match(obj_id_comb, contam_id_comb)]
    object_meta_data$Contamination_FLAG[is.na(object_meta_data$Contamination_FLAG)] <- "PASS"

    clone_cont_df <- object_meta_data[object_meta_data$Contamination_FLAG != "PASS", c(sample_col, "igClonotypeID")]
    for(row in 1:nrow(clone_cont_df)){
      table_out <- table(object_meta_data$Contamination_FLAG[object_meta_data[[sample_col]] == clone_cont_df[row, 1] & object_meta_data$igClonotypeID == clone_cont_df[row, 2]])
      if(table_out[["OUT"]]/sum(table_out)*100 >= contamination_clone_cutoff){
        object_meta_data$Contamination_FLAG[object_meta_data[[sample_col]] == clone_cont_df[row, 1] & object_meta_data$igClonotypeID == clone_cont_df[row, 2]] <- "OUT_CLONOTYPE"
      }
    }
  }
  return(object_meta_data)
}
