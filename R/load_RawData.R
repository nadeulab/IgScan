#' Load IG contigs from NGS data
#'
#' This function loads raw data from various common input formats and
#' returns a list of elements that can be used as the 'samples_list'
#' parameter in other IgScan functions.
#'
#' It supports data from multiple sequencing platforms and file types:
#' \itemize{
#'   \item 10x BCR =  "filtered_contig_annotations.csv / filtered_contig.fa"
#'   \item Parse BCR =  "bcr_annotation_airr.tsv"
#'   \item BD Rhapsody BCR =  "Contigs_AIRR.tsv"
#'   \item MiXCR =  "clonotypes.IGX.txt"
#'   \item TRUST4 = "barcode_report.tsv"
#'   \item AIRR = "airr_rearrangement.tsv"
#'   \item IMGT AIRR = "vquest_airr.tsv"
#'   \item fasta
#' }
#'
#' @param sample_paths A vector of paths to input files or a vector of
#' directories containing the input file. Supported formats depend on the
#' `input_format` parameter.
#' @param input_format A string specifying the format of the input data, currently supporting:
#' '10xBCR_fasta', '10xBCR_csv', 'ParseBCR', 'BDRhapsodyBCR', 'MiXCR',
#' 'TRUST4', 'AIRR', 'IMGT_AIRR' and 'fasta'.
#'
#' @return A list of input elementsfor further functions in IgScan.
#'   For 10xBCR_fasta and fasta formats, the raw file paths are returned as elements.
#'   For other sources, data frames containing the loaded and parsed data are returned.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' paths <- c("path/to/bcr_annotation_airr_file1.tsv", "path/to/bcr_annotation_airr_file2.tsv")
#' samples_list <- load_RawData(paths, "ParseBCR")
#' }
#'
load_RawData <- function(sample_paths, input_format){

  input_format <- tolower(input_format)
  if(!input_format %in% c("10xbcr_fasta", "10xbcr_csv", "parsebcr", "bdrhapsodybcr", "mixcr", "trust4", "airr", "imgt_airr", "fasta")){stop("Invalid value for 'input_format'. It should be either '10xbcr_fasta', '10xbcr_csv', 'parsebcr', 'bdrhapsodybcr', 'mixcr', 'trust4', 'airr', 'imgt_airr' or 'fasta'.")}

  info <- file.info(sample_paths)

  sample_paths_read <- c()
  for(row in 1:nrow(info)){
    if(info$isdir[row]){
      sample_paths_read <- c(sample_paths_read, list.files(rownames(info)[row], full.names = T, recursive = F))
    } else{
      if(!file.exists(rownames(info)[row])){stop(paste0("File ", rownames(info)[row], " does not exist!"))}
      sample_paths_read <- c(sample_paths_read, rownames(info)[row])
    }
  }

  samples_list <- list()
  for(i in 1:length(sample_paths_read)){
    if(input_format %in% c("10xbcr_fasta", "fasta")){
      tmp <- sample_paths_read[i]
    }
    else if(input_format %in% c("mixcr", "bdrhapsodybcr", "parsebcr", "trust4", "airr", "imgt_airr")){
      tmp <- fread(sample_paths_read[i], header = T, sep = "\t", stringsAsFactors = F, data.table = F)
    }
    else if(input_format == "10x_csv"){
      tmp <- fread(sample_paths_read[i], header = T, sep = ",", stringsAsFactors = F, data.table = F)
    }
    samples_list[[i]] <- tmp
  }
  return(samples_list)
}
