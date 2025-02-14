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
  } else if(class(object)[1] == "Seurat"){
    data_frame <- object@meta.data
    data_frame <- .extract_IgScanDf_from_SingleCell(data_frame)
  } else{ ## PENSAR SI AÑADIR LA CLASE IgSCAN a estas cosas!
    data_frame <- object
  }

  if(data_mode == "complete"){

  } else if(data_mode == "for_phylo"){

  }

  return(converted_object)
}
