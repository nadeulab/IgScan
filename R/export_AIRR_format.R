#' Export IgScan outputs in AIRR format
#'
#' This function converts the output from IgScan (or compatible objects single cell
#' objects with IgScan annotation) into the standardized AIRR format and writes it to
#' a `.tsv` file. It supports input objects from `IgScan output data.frame`, `Seurat`,
#' or `SingleCellExperiment` classes.
#'
#' @param object A data frame, `Seurat` object, or `SingleCellExperiment` object
#' containing IgScan annotation fields. For data frames, it is expected to match
#' the IgScan format (specific column names will be validated).
#' @param dir A character string specifying the output directory. If it does not exist, it will be created.
#' @param fileName A character string for the output file name. Defaults to `"IgScan_AIRR_formatted.tsv"` if not provided.
#' @param germline_aln Type of germline alignment to be included in the AIRR file.
#' Options are `masked` (uses IgBlast masked germline alignment) and `consensus` (uses the IgScan custom-built
#' consensus germline sequence).
#'
#' @return A `.tsv` file written to the specified directory in AIRR format as well
#' as the AIRR dataframe as an R object.
#'
#' @export
#'
#' @importFrom data.table fwrite
#' @import Seurat
#' @import SingleCellExperiment
#' @importFrom stringr str_split_fixed
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' # Export from a data.frame
#' export_AIRR_format(igscan_df, dir = "results/")
#'
#' # Export from Seurat object
#' export_AIRR_format(seurat_obj, dir = "results/", germline_aln = "consensus")
#' }
#'
export_AIRR_format <- function(object, dir, fileName = NULL, germline_aln = "masked"){

  germline_aln <- tolower(germline_aln)
  if(!germline_aln %in% c("masked", "consensus")){stop("Invalid value for `germline_aln`. Please set an option between `masked` or `consensus`.")}
  if(!endsWith(dir, "/")){dir <- paste0(dir, "/")}
  if(!dir.exists(dir)){
    message("The indicated output directory does not exist. Creating it...")
    system(paste0("mkdir ", dir))
  }

  converted_object <- .convert_IgScan_to_AIRR(object, germline_aln)

  if(is.null(fileName)){
    message("No fileName specified. Output file will be named as `IgScan_AIRR_formatted.tsv`.")
    fileName <- "IgScan_AIRR_formatted.tsv"
  }

  fwrite(converted_object, file = paste0(dir, fileName), quote = F, sep = "\t", col.names = T, row.names = F)
  return(converted_object)
}

.convert_IgScan_to_AIRR <- function(object, germline_aln){

  if(class(object)[1] == "SingleCellExperiment"){
    data_frame <- colData(object)
    data_frame <- .extract_IgScanDf_from_SingleCell(data_frame)

  } else if(class(object)[1] == "Seurat"){
    data_frame <- object@meta.data
    data_frame <- .extract_IgScanDf_from_SingleCell(data_frame)

  } else{
    data_frame <- object
    colnames_dictionary_path <- system.file("colnames_dictionary.RData", package = "IgScan", mustWork = T)
    load(colnames_dictionary_path)
    if(!all(colnames_dictionary$FinalNames[c(1:18, 20:21, 25, 28, 51)] %in% colnames(data_frame))){stop("The provided object does not contain the expected IgScan fields. Please provide a valid IgScan output file as input.")}
  }

  converted_object <- data.frame(sequence_id = data_frame[, "contig_id"])
  converted_object$sequence <- data_frame$Raw_sequence
  converted_object$locus <- substr(data_frame$VDJ_genes, start = 1, stop = 3)
  converted_object$productive <- ifelse(data_frame$Functionality == "productive", "T", "F")
  converted_object$rev_comp <- "F"

  converted_object$sequence_alignment <- data_frame$VDJ_sequence_correctedCDR3
  seq_aa_translate <- sapply(unique(data_frame$VDJ_sequence_correctedCDR3), .translate_sequence)
  seq_aa_df <- data.frame(nt = names(seq_aa_translate), aa = unname(seq_aa_translate))
  converted_object$sequence_alignment_aa <- seq_aa_df$aa[match(converted_object$sequence_alignment, seq_aa_df$nt)]

  if(germline_aln == "masked"){
    converted_object$germline_alignment <- data_frame$IgBlast_Germline_alignment
    germ_aa_translate <- sapply(unique(data_frame$IgBlast_Germline_alignment), .translate_sequence)
    germ_aa_df <- data.frame(nt = names(germ_aa_translate), aa = unname(germ_aa_translate))
    converted_object$germline_alignment_aa <- germ_aa_df$aa[match(converted_object$germline_alignment, germ_aa_df$nt)]

  } else if(germline_aln == "consensus"){
    converted_object$germline_alignment <- data_frame$Consensus_Germline
    converted_object$germline_alignment_aa <- data_frame$Consensus_Germline_aa
  }

  converted_object$junction_aa <- data_frame$Junction_aa

  split_vdj_df <- as.data.frame(str_split_fixed(string = data_frame$VDJ_genes, pattern = "/", n = 3))
  split_vdj_df$V3[grepl("IG[HKL]J", split_vdj_df$V2)] <- split_vdj_df$V2[grepl("IG[HKL]J", split_vdj_df$V2)]
  split_vdj_df$V2[grepl("IG[HKL]J", split_vdj_df$V2)] <- ""

  converted_object$v_call <- split_vdj_df$V1
  converted_object$d_call <- split_vdj_df$V2
  converted_object$j_call <- split_vdj_df$V3
  converted_object$c_call <- data_frame$C_gene
  converted_object[, c("v_cigar", "d_cigar", "j_cigar")] <- ""

  split_pos_df <- as.data.frame(str_split_fixed(string = data_frame$VDJ_positions, pattern = "-", n = 7))
  colnames(split_pos_df) <- c("fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4")

  converted_object$fwr1_start <- 1
  converted_object$cdr1_start <- converted_object$fwr1_start + as.numeric(split_pos_df$fr1)
  converted_object$fwr2_start <- converted_object$cdr1_start + as.numeric(split_pos_df$cdr1)
  converted_object$cdr2_start <- converted_object$fwr2_start + as.numeric(split_pos_df$fr2)
  converted_object$fwr3_start <- converted_object$cdr2_start + as.numeric(split_pos_df$cdr2)
  converted_object$cdr3_start <- converted_object$fwr3_start + as.numeric(split_pos_df$fr3)
  converted_object$fwr4_start <- converted_object$cdr3_start + as.numeric(split_pos_df$cdr3)

  converted_object$fwr1_end <- converted_object$cdr1_start - 1
  converted_object$cdr1_end <- converted_object$fwr2_start - 1
  converted_object$fwr2_end <- converted_object$cdr2_start - 1
  converted_object$cdr2_end <- converted_object$fwr3_start - 1
  converted_object$fwr3_end <- converted_object$cdr3_start - 1
  converted_object$cdr3_end <- converted_object$fwr4_start - 1
  converted_object$fwr4_end <- nchar(converted_object$sequence_alignment)

  subset_df <- converted_object[!duplicated(converted_object$sequence_alignment),]
  subset_df <- subset_df %>%
    rowwise() %>%
    mutate(
      fwr1 = substr(sequence_alignment, fwr1_start, fwr1_end),
      fwr1_aa = substr(sequence_alignment_aa, fwr1_start, fwr1_end / 3),
      cdr1 = substr(sequence_alignment, cdr1_start, cdr1_end),
      cdr1_aa = substr(sequence_alignment_aa, fwr1_end / 3 + 1, cdr1_end / 3),
      fwr2 = substr(sequence_alignment, fwr2_start, fwr2_end),
      fwr2_aa = substr(sequence_alignment_aa, cdr1_end / 3 + 1, fwr2_end / 3),
      cdr2 = substr(sequence_alignment, cdr2_start, cdr2_end),
      cdr2_aa = substr(sequence_alignment_aa, fwr2_end / 3 + 1, cdr2_end / 3),
      fwr3 = substr(sequence_alignment, fwr3_start, fwr3_end),
      fwr3_aa = substr(sequence_alignment_aa, cdr2_end / 3 + 1, fwr3_end / 3),
      cdr3 = substr(sequence_alignment, cdr3_start, cdr3_end),
      cdr3_aa = substr(sequence_alignment_aa, fwr3_end / 3 + 1, cdr3_end / 3),
      fwr4 = substr(sequence_alignment, fwr4_start, fwr4_end),
      fwr4_aa = substr(sequence_alignment_aa, cdr3_end / 3 + 1, nchar(sequence_alignment_aa)),
      junction = substr(sequence_alignment, cdr3_start - 3, cdr3_end + 3)
    ) %>%
    ungroup()

  merge_data <- subset_df[match(converted_object$sequence_alignment, subset_df$sequence_alignment),
                          c("fwr1", "fwr1_aa", "cdr1", "cdr1_aa", "fwr2", "fwr2_aa", "cdr2", "cdr2_aa", "fwr3", "fwr3_aa", "cdr3", "cdr3_aa", "fwr4", "fwr4_aa", "junction")]

  converted_object <- cbind(converted_object, merge_data)
  converted_object$junction_length <- nchar(converted_object$junction)

  if("Subclone_nReads" %in% colnames(data_frame)){
    converted_object$duplicate_count <- data_frame$Subclone_nReads[match(converted_object$sequence_alignment, data_frame$VDJ_sequence_correctedCDR3)]
    converted_object$clone_id <- data_frame$ClonotypeID[match(converted_object$sequence_alignment, data_frame$VDJ_sequence_correctedCDR3)]
    converted_object$sample_id <- data_frame$SampleID[match(converted_object$sequence_alignment, data_frame$VDJ_sequence_correctedCDR3)]

  } else{ ## WORK ON THIS!!!
    converted_object$duplicate_count <- 1
    converted_object$clone_id <- data_frame$ClonotypeID[match(converted_object$sequence_alignment, data_frame$VDJ_sequence_correctedCDR3)]
    converted_object$sample_id <- data_frame$SampleID[match(converted_object$sequence_alignment, data_frame$VDJ_sequence_correctedCDR3)]
  }

  colnames_AIRR <- c(
    "sequence_id", "sequence", "locus", "productive", "rev_comp",
    "junction", "junction_aa", "junction_length",
    "v_call", "d_call", "j_call", "c_call",
    "sequence_alignment", "sequence_alignment_aa", "germline_alignment", "germline_alignment_aa",
    "v_cigar", "d_cigar", "j_cigar",
    "fwr1", "fwr1_aa", "cdr1", "cdr1_aa", "fwr2", "fwr2_aa", "cdr2", "cdr2_aa", "fwr3", "fwr3_aa", "cdr3", "cdr3_aa", "fwr4", "fwr4_aa",
    "fwr1_start", "fwr1_end", "cdr1_start", "cdr1_end", "fwr2_start", "fwr2_end", "cdr2_start", "cdr2_end", "fwr3_start", "fwr3_end", "cdr3_start", "cdr3_end", "fwr4_start", "fwr4_end",
    "duplicate_count", "clone_id", "sample_id")

  converted_object <- converted_object[, colnames_AIRR]

  return(converted_object)
}
