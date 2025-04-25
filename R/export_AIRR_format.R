#' Export IgScan outputs in AIRR format
#'
#' @param object
#' @param dir
#' @param fileName
#'
#' @return
#' @export
#'
#' @importFrom data.table fwrite
#'
#' @examples
#'
export_AIRR_format <- function(object, dir, fileName = NULL){

  if(!endsWith(dir, "/")){dir <- paste0(dir, "/")}
  if(!dir.exists(dir)){stop("The indicated output directory does not exist. Please, set a valid directory for output.")}

  converted_object <- .convert_IgScan_to_AIRR(object)

  if(is.null(fileName)){
    fileName <- "IgScan_AIRR_formatted.tsv"
  }

  fwrite(converted_object, file = paste0(dir, fileName), quote = F, sep = "\t", col.names = T, row.names = F)
}

.convert_IgScan_to_AIRR <- function(object, data_mode = "complete"){

  if(class(object)[1] == "SingleCellExperiment"){
    data_frame <- object@colData
    data_frame <- .extract_IgScanDf_from_SingleCell(data_frame)
    id_col <- "contig_id"

  } else if(class(object)[1] == "Seurat"){
    data_frame <- object@meta.data
    data_frame <- .extract_IgScanDf_from_SingleCell(data_frame)
    id_col <- "contig_id"

  } else{ ## PENSAR SI AÑADIR LA CLASE IgSCAN a estas cosas!
    data_frame <- object
    id_col <- "sequence_id"
  }

  converted_object <- data.frame(sequence_id = data_frame[,id_col])
  converted_object$sequence <- data_frame$Raw_sequence
  converted_object$locus <- substr(data_frame$VDJ_genes, start = 1, stop = 3)
  converted_object$productive <- ifelse(data_frame$Functionality == "productive", "T", "F")
  converted_object$rev_comp <- "F"

  converted_object$junction <- ###########
  converted_object$junction_aa <- data_frame$Junction_aa

    c("sequence_id", "sequence", "locus", "productive", "rev_comp", "junction", "junction_aa",
    "v_call", "d_call", "j_call", "c_call", "sequence_alignment", "sequence_alignment_aa", "germline_alignment", "germline_alignment_aa",
    "v_cigar", "d_cigar", "j_cigar",
    "fwr1", "fwr1_aa", "cdr1", "cdr1_aa", "fwr2", "fwr2_aa", "cdr2", "cdr2_aa", "fwr3", "fwr3_aa", "cdr3", "cdr3_aa", "fwr4", "fwr4_aa",
    "fwr1_start", "fwr1_end", "cdr1_start", "cdr1_end", "fwr2_start", "fwr2_end", "cdr2_start", "cdr2_end", "fwr3_start", "fwr3_end", "cdr3_start", "cdr3_end", "fwr4_start", "fwr4_end")


  if(data_mode == "complete"){

  } else if(data_mode == "for_phylo"){

  }

  return(converted_object)
}
