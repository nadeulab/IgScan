#' Filter SampleTag Contamination in BD Rhapsody Data
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
#' Default is 10, higher values result in stricter contamination detection.
#' @param remove_contamination Logical. If `TRUE`, cells flagged as contaminated are
#' removed from the returned object (default is `FALSE`).
#' @param recalc_column A string/vector of strings with the name of the column/s to use for ID
#' recalculation. Needed when remove_contamination=TRUE. Default is `orig.ident`.
#'
#' @return The input `single_cell_object` with updated metadata, which includes three new columns:
#' \itemize{
#'   \item Contamination_FLAG: Indicates whether a cell passed contamination checks (`PASS`), was flagged as contaminated (`OUT`),
#'   or was involved in crossed contamination (`OUT_X`)
#'   \item Contamination_Ratio: The ratio of clonal overlap between the contaminated and dominant samples.
#'   \item Contamination_Sample: The dominant sample responsible for contamination.
#' }
#'
#'   If `remove_contamination = TRUE`, cells flagged as `OUT` or `OUT_X` are removed from
#'   the data frame and immunogenetic IDs are recalculated with `recalculate_IDs_single_cell()`.
#'
#' @export
#'
#' @importFrom tidyr pivot_wider
#' @import Seurat
#' @import SeuratObject
#' @import SingleCellExperiment
#'
#' @examples
#' \dontrun{
#' # Example with a SingleCellExperiment object:
#' cont_filtered_sce <- filter_SampleTag_contam_Rhapsody(single_cell_object = sce, contamination_cutoff = 30)
#' }
#'
filter_SampleTag_contam_Rhapsody <- function(single_cell_object, contamination_cutoff = 10, remove_contamination = F, recalc_column = "orig.ident"){

  if(class(single_cell_object)[1] == "SingleCellExperiment"){
    meta_data <- colData(single_cell_object)
  } else if(class(single_cell_object)[1] == "Seurat"){
    meta_data <- single_cell_object@meta.data
  }

  if(!"Sample_Tag" %in% colnames(meta_data)){stop(paste0("\n`Sample_Tag` column not found in the meta_data! Please, ensure to include the Sample_Tag column before executing this function."))}
  if(!all(recalc_column %in% colnames(meta_data)) & remove_contamination){stop(paste0("\nUnknown recalc_column column (", recalc_column[!recalc_column %in% colnames(meta_data)], ") selected for flagging SampleTag contamination! Please, set a valid column name."))}

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
      }
    }
  }

  if(class(single_cell_object)[1] == "SingleCellExperiment"){
    colData(single_cell_object)$Contamination_FLAG <- meta_data$Contamination_FLAG[match(rownames(colData(single_cell_object)), rownames(meta_data))]
    colData(single_cell_object)$Contamination_Ratio <- meta_data$Contamination_Ratio[match(rownames(colData(single_cell_object)), rownames(meta_data))]
    colData(single_cell_object)$Contamination_Sample <- meta_data$Contamination_Sample[match(rownames(colData(single_cell_object)), rownames(meta_data))]

  } else if(class(single_cell_object)[1] == "Seurat"){
    single_cell_object@meta.data$Contamination_FLAG <- meta_data$Contamination_FLAG[match(rownames(single_cell_object@meta.data), rownames(meta_data))]
    single_cell_object@meta.data$Contamination_Ratio <- meta_data$Contamination_Ratio[match(rownames(single_cell_object@meta.data), rownames(meta_data))]
    single_cell_object@meta.data$Contamination_Sample <- meta_data$Contamination_Sample[match(rownames(single_cell_object@meta.data), rownames(meta_data))]
  }

  if(remove_contamination){
    if(class(single_cell_object)[1] == "SingleCellExperiment"){
      colData(single_cell_object)$Contamination_FLAG[is.na(colData(single_cell_object)$Contamination_FLAG)] <- "."
      single_cell_object <- subset(single_cell_object, subset = Contamination_FLAG != "OUT")
      colData(single_cell_object)$Contamination_FLAG[colData(single_cell_object)$Contamination_FLAG == "."] <- NA

    } else if(class(single_cell_object)[1] == "Seurat"){
      single_cell_object@meta.data$Contamination_FLAG[is.na(single_cell_object@meta.data$Contamination_FLAG)] <- "."
      single_cell_object <- subset(single_cell_object, subset = Contamination_FLAG != "OUT")
      single_cell_object@meta.data$Contamination_FLAG[single_cell_object@meta.data$Contamination_FLAG == "."] <- NA
    }
    single_cell_object <- recalculate_IDs_single_cell(single_cell_object, group_col = recalc_column)
  }
  return(single_cell_object)
}
