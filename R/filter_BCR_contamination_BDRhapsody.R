#' Filter BCR Contamination in BD Rhapsody Data
#'
#' This function identifies, flags, and optionally removes cell contamination between sample tags
#' in an IgScan-annotated single-cell object derived from BD Rhapsody data. Contamination is assessed
#' by analyzing the distribution of BCR clonotypes between sample tags. A contamination cutoff can be
#' specified, indicating that if a clonotype is observed X times more in one sample tag than another,
#' it will be flagged as contamination in the latter sample tag. The function works with both
#' `SingleCellExperiment` and `Seurat` objects.
#'
#' @param single_cell_object A `SingleCellExperiment` or `Seurat` object containing IgScan
#' annotation for single cell data. The object must also include the `Sample_Tag` identifier.
#' @param contamination_cutoff A numeric value specifying the contamination ratio threshold.
#'   Default is 10, higher values result in stricter contamination detection.
#' @param remove_contamination Logical. If `TRUE`, cells flagged as contaminated are
#'   removed from the returned object (default is `FALSE`).
#'
#' @return The input `single_cell_object` with updated metadata, which includes three new columns:
#' \itemize{
#'   \item Contamination_FLAG: Indicates whether a cell passed contamination checks (`PASS`), was flagged as contaminated (`OUT`),
#'   or was involved in crossed contamination (`OUT_X`)
#'   \item Contamination_Ratio: The ratio of clonal overlap between the contaminated and dominant samples.
#'   \item Contamination_Sample: The dominant sample responsible for contamination.
#' }
#'
#'   If `remove_contamination = TRUE`, cells flagged as `OUT` or `OUT_X` are removed from the object.
#'
#' @export
#'
#' @importFrom tidyr pivot_wider
#'
#' @examples
#' \dontrun{
#' # Example with a SingleCellExperiment object:
#' cont_filtered_sce <- filter_BCR_contamination_BDRhapsody(single_cell_object = sce, contamination_cutoff = 30)
#' }
#'
filter_BCR_contamination_BDRhapsody <- function(single_cell_object, contamination_cutoff = 10, remove_contamination = F){

  if(class(single_cell_object)[1] == "SingleCellExperiment"){
    meta_data <- single_cell_object@colData
  } else if(class(single_cell_object)[1] == "Seurat"){
    meta_data <- single_cell_object@meta.data
  }

  meta_data <- single_cell_object@meta.data[!single_cell_object@meta.data$Sample_Tag %in% c("Multiplet", "Undetermined"),]
  meta_data$Sample_Tag <- as.character(meta_data$Sample_Tag)

  meta_data$Contamination_FLAG <- "PASS"
  meta_data$Contamination_Ratio <- NA
  meta_data$Contamination_Sample <- NA

  tab <- as.data.frame(table(meta_data$igClonotypeID, meta_data$Sample_Tag))
  tab <- as.data.frame(tab %>% pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0))
  colnames(tab)[1] <- "igClonotypeID"
  tab$igClonotypeID <- as.character(tab$igClonotypeID)
  tab$total <- rowSums(tab[,2:ncol(tab)])
  tab <- tab[tab$total > 1,]

  sampleTags <- colnames(tab)[2:(ncol(tab)-1)]

  for(clone in tab$igClonotypeID){

    max_sampleTag <- sampleTags[which.max(tab[tab$igClonotypeID == clone, sampleTags])]
    max_Value <- tab[tab$igClonotypeID == clone, max_sampleTag]

    for(other_ST in sampleTags[sampleTags != max_sampleTag]){

      other_Value <- tab[tab$igClonotypeID == clone, other_ST]
      if(other_Value == 0){next}

      cont_ratio <- max_Value / other_Value

      if(cont_ratio > contamination_cutoff){
        meta_data$Contamination_FLAG[meta_data$Sample_Tag == other_ST & meta_data$igClonotypeID == clone] <- "OUT"
        meta_data$Contamination_Ratio[meta_data$Sample_Tag == other_ST & meta_data$igClonotypeID == clone] <- paste0(other_Value,":",max_Value)
        meta_data$Contamination_Sample[meta_data$Sample_Tag == other_ST & meta_data$igClonotypeID == clone] <- max_sampleTag

      } else{
        meta_data$Contamination_FLAG[meta_data$Sample_Tag %in% c(other_ST, max_sampleTag) & meta_data$igClonotypeID == clone] <- "OUT_X"
        meta_data$Contamination_Ratio[meta_data$Sample_Tag %in% c(other_ST, max_sampleTag) & meta_data$igClonotypeID == clone] <- paste0(other_Value,":",max_Value)
      }
    }
  }

  if(class(single_cell_object)[1] == "SingleCellExperiment"){
    single_cell_object@colData$Contamination_FLAG <- meta_data$Contamination_FLAG[match(rownames(single_cell_object@colData), rownames(meta_data))]
    single_cell_object@colData$Contamination_Ratio <- meta_data$Contamination_FLAG[match(rownames(single_cell_object@colData), rownames(meta_data))]
    single_cell_object@colData$Contamination_Sample <- meta_data$Contamination_FLAG[match(rownames(single_cell_object@colData), rownames(meta_data))]

  } else if(class(single_cell_object)[1] == "Seurat"){
    single_cell_object@meta.data$Contamination_FLAG <- meta_data$Contamination_FLAG[match(rownames(single_cell_object@meta.data), rownames(meta_data))]
    single_cell_object@meta.data$Contamination_Ratio <- meta_data$Contamination_FLAG[match(rownames(single_cell_object@meta.data), rownames(meta_data))]
    single_cell_object@meta.data$Contamination_Sample <- meta_data$Contamination_FLAG[match(rownames(single_cell_object@meta.data), rownames(meta_data))]
  }

  if(remove_contamination){single_cell_object <- subset(single_cell_object, subset = Contamination_FLAG == "PASS")}

  return(single_cell_object)
}
