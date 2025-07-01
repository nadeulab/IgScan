#' Prepare gctree inputs from IgScan annotated datasets
#'
#' This function takes as input an IgScan annotated dataframe coming from bulk
#' NGS or single cell workflows, or a single cell object with IgScan annotation,
#' and prepares the input files for running the `gctree` program for BCR phylogenetic
#' inference. Once the input files have been generated, the user can run `gctree`
#' following the instructions in XXXX jupyter notebook provided with the IgScan package.
#'
#' The function generates the following files in the specified `outDir`:
#' \itemize{
#'   \item `abundances.csv`: Abundance of each subclone including the germline (GL).
#'   \item `idmap.txt`: Mapping of subclone IDs to individual barcodes (or rownames).
#'   \item `subclones.fasta`: FASTA file of the aligned somatic and germline VDJ sequences.
#'   \item `deduplicated.phylip`: Alignment file in PHYLIP format for input to `gctree`.
#' }
#'
#' @param object A data.frame, `SingleCellExperiment`, or `Seurat` object with IgScan annotation.
#' @param outDir A string. Path to the directory where output files should be written.
#' @param mode A string. Either `"single_cell"` or `"bulk"` depending on the dataset type. Default is `"single_cell"`.
#' @param cloneID A string. The clonotype identifier to extract and process. Note that the `cloneID` parameter must
#' match a value in the clonotype identifier field in the IgScan objects.
#' @param chain A string. One of `"H"` (heavy chain), `"L"` (light chain), or `"HL"` (paired chains). Default is `"HL"`.
#' Only needed when `mode` is set to `"single_cell"`.
#'
#' @return The function returns no value in R but writes the necessary input files for `gctree` into the `outDir`.
#'
#' @export
#'
#' @importFrom ape read.dna write.dna
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' # Example with a SingleCellExperiment object
#' prepare_gctree_input(sce_object, outDir = "results/clone1/", mode = "single_cell", cloneID = "C1", chain = "HL")
#'
#' # Example with a Seurat object
#' prepare_gctree_input(seurat_object, outDir = "results/clone2/", mode = "single_cell", cloneID = "C5", chain = "H")
#'
#' # Example with a bulk IgScan-annotated dataframe
#' prepare_gctree_input(bulk_df, outDir = "results/clone3/", mode = "bulk", cloneID = "CloneA")
#' }
#'
prepare_gctree_input <- function(object, outDir, mode = "single_cell", cloneID = "C1", chain = "HL"){

  if(!mode %in% c("bulk", "single_cell")){stop("Invalid value for `mode`. Please set an option between `bulk` or `single_cell`.")}
  if(!chain %in% c("H", "L", "HL")){stop("Invalid value for `chain`. Please set an option between `H`, `L` or `HL`.")}

  if(class(object)[1] == "SingleCellExperiment"){
    data_frame <- colData(object)
    if(!mode == "single_cell"){stop("Unmatched object class with the `mode` option specified. Provided a `single_cell` object but specified `bulk` mode.")}

  } else if(class(object)[1] == "Seurat"){
    if(!mode == "single_cell"){stop("Unmatched object class with the `mode` option specified. Provided a `single_cell` object but specified `bulk` mode.")}
    data_frame <- object@meta.data

  } else{
    data_frame <- object
  }

  if(!endsWith(outDir, "/")){outDir <- paste0(outDir, "/")}
  if(!dir.exists(outDir)){
    message("The indicated output directory does not exist. Creating it...")
    system(paste0("mkdir ", outDir))
  }

  if(mode == "single_cell"){

    if(!"completeBCR" %in% colnames(data_frame)){stop("The provided object does not contain the expected IgScan fields. Please provide a valid IgScan output file as input.")}

    if(!cloneID %in% unique(data_frame$igClonotypeID_num)){stop("The specified `cloneID` is not present in the input data.")}
    data_frame_subset <- data_frame[data_frame$igClonotypeID_num %in% cloneID,]
    if(nrow(data_frame_subset) < 2){stop("There are not enough sequences for phylogenetic inference in the chosen clonotype.")}

    if(chain == "HL"){
      data_frame_subset <- data_frame_subset[data_frame_subset$completeBCR == "Yes",]

      abundance_df <- rbind(c("GL", 0), data_frame_subset %>%
                                          count(igSubcloneID_in_Clonotype_num, name = "Abundance"))

      write.table(abundance_df, file = paste0(outDir, "abundances.csv"), quote = F, sep = ",", row.names = F, col.names = F)

      sink(paste0(outDir, "idmap.txt"))
      for(sbc in abundance_df$igSubcloneID_in_Clonotype_num){
        bcd <- paste(rownames(data_frame_subset[data_frame_subset$igSubcloneID_in_Clonotype_num == sbc,]), collapse = ":")
        cat(paste0(sbc, ",", bcd, "\n"))
      }
      sink()

      data_frame_subset <- data_frame_subset[!duplicated(data_frame_subset$igSubcloneID_in_Clonotype_num),]

      somatic_seqs <- data_frame_subset[, paste0(c("IGH1", "IGH2", "IGL1", "IGL2"), "_VDJ_sequence_correctedCDR3")]
      somatic_seqs[is.na(somatic_seqs)] <- ""
      data_frame_subset$seq_write <- as.character(apply(somatic_seqs, 1, paste, collapse = ""))

      sink(paste0(outDir, "subclones_", chain, "_", cloneID, ".fasta"))
      germline <- gsub("_", "", unique(data_frame_subset$igClonotype_Consensus_Germline))
      cat(paste0(">GL      \n", germline, "\n"))
      for(row in 1:nrow(data_frame_subset)){
        sbc_id <- data_frame_subset$igSubcloneID_in_Clonotype_num[row]
        sbc_id_fasta <- paste0(sbc_id, strrep(" ", 9-nchar(sbc_id)))
        cat(paste0(">", sbc_id_fasta, "\n", data_frame_subset$seq_write[row], "\n"))
      }
      sink()

      aln <- read.dna(paste0(outDir, "subclones_", chain, "_", cloneID, ".fasta"), format = "fasta")
      write.dna(aln, file = paste0(outDir, "deduplicated.phylip"), format = "interleaved", nbcol = 5, indent = 10)

    } else if(chain == "H"){
      data_frame_subset <- data_frame_subset[data_frame_subset$completeBCR %in% c("Yes", "Yes_rescue"),]

      data_frame_subset$Chain_SbcID <- gsub("_NA", "", paste(data_frame_subset$IGH1_SubcloneID, data_frame_subset$IGH2_SubcloneID, sep = "_"))

      ## Write Abundance file
      abundance_df <- rbind(c("GL", 0), data_frame_subset %>%
                              count(Chain_SbcID, name = "Abundance"))

      abundance_df$Chain_SbcID_OLD <- abundance_df$Chain_SbcID
      abundance_df$Chain_SbcID[2:nrow(abundance_df)] <- paste0(gsub("\\.CV.\\.S.", "", abundance_df$Chain_SbcID[2:nrow(abundance_df)]), ".", 1:(nrow(abundance_df)-1))

      write.table(abundance_df[,1:2], file = paste0(outDir, "abundances.csv"), quote = F, sep = ",", row.names = F, col.names = F)

      data_frame_subset$Chain_SbcID <- abundance_df$Chain_SbcID[match(data_frame_subset$Chain_SbcID, abundance_df$Chain_SbcID_OLD)]

      ## Write IDMAP file
      sink(paste0(outDir, "idmap.txt"))
      for(sbc in abundance_df$Chain_SbcID){
        bcd <- paste(rownames(data_frame_subset[data_frame_subset$Chain_SbcID == sbc,]), collapse = ":")
        cat(paste0(sbc, ",", bcd, "\n"))
      }
      sink()

      data_frame_subset <- data_frame_subset[!duplicated(data_frame_subset$Chain_SbcID),]

      somatic_seqs <- data_frame_subset[, paste0(c("IGH1", "IGH2"), "_VDJ_sequence_correctedCDR3")]
      somatic_seqs[is.na(somatic_seqs)] <- ""
      data_frame_subset$seq_write <- as.character(apply(somatic_seqs, 1, paste, collapse = ""))

      sink(paste0(outDir, "subclones_", chain, "_", cloneID, ".fasta"))
      germline <- gsub("_NA", "", paste0(unique(data_frame_subset$IGH1_Consensus_Germline), "_", unique(data_frame_subset$IGH2_Consensus_Germline)))
      cat(paste0(">GL      \n", germline, "\n"))
      for(row in 1:nrow(data_frame_subset)){
        sbc_id <- data_frame_subset$Chain_SbcID[row]
        sbc_id_fasta <- paste0(sbc_id, strrep(" ", 9-nchar(sbc_id)))
        cat(paste0(">", sbc_id_fasta, "\n", data_frame_subset$seq_write[row], "\n"))
      }
      sink()

      aln <- read.dna(paste0(outDir, "subclones_", chain, "_", cloneID, ".fasta"), format = "fasta")
      write.dna(aln, file = paste0(outDir, "deduplicated.phylip"), format = "interleaved", nbcol = 5, indent = 10)

    } else if(chain == "L"){
      data_frame_subset <- data_frame_subset[data_frame_subset$completeBCR %in% c("Yes", "Yes_rescue"),]

      data_frame_subset$Chain_SbcID <- gsub("_NA", "", paste(data_frame_subset$IGL1_SubcloneID, data_frame_subset$IGL2_SubcloneID, sep = "_"))

      ## Write Abundance file
      abundance_df <- rbind(c("GL", 0), data_frame_subset %>%
                              count(Chain_SbcID, name = "Abundance"))

      abundance_df$Chain_SbcID_OLD <- abundance_df$Chain_SbcID
      abundance_df$Chain_SbcID[2:nrow(abundance_df)] <- paste0(gsub("\\.CV.\\.S.", "", abundance_df$Chain_SbcID[2:nrow(abundance_df)]), ".", 1:(nrow(abundance_df)-1))

      write.table(abundance_df[,1:2], file = paste0(outDir, "abundances.csv"), quote = F, sep = ",", row.names = F, col.names = F)

      data_frame_subset$Chain_SbcID <- abundance_df$Chain_SbcID[match(data_frame_subset$Chain_SbcID, abundance_df$Chain_SbcID_OLD)]

      ## Write IDMAP file
      sink(paste0(outDir, "idmap.txt"))
      for(sbc in abundance_df$Chain_SbcID){
        bcd <- paste(rownames(data_frame_subset[data_frame_subset$Chain_SbcID == sbc,]), collapse = ":")
        cat(paste0(sbc, ",", bcd, "\n"))
      }
      sink()

      data_frame_subset <- data_frame_subset[!duplicated(data_frame_subset$Chain_SbcID),]

      somatic_seqs <- data_frame_subset[, paste0(c("IGL1", "IGL2"), "_VDJ_sequence_correctedCDR3")]
      somatic_seqs[is.na(somatic_seqs)] <- ""
      data_frame_subset$seq_write <- as.character(apply(somatic_seqs, 1, paste, collapse = ""))

      sink(paste0(outDir, "subclones_", chain, "_", cloneID, ".fasta"))
      germline <- gsub("_NA", "", paste0(unique(data_frame_subset$IGL1_Consensus_Germline), "_", unique(data_frame_subset$IGL2_Consensus_Germline)))
      cat(paste0(">GL      \n", germline, "\n"))
      for(row in 1:nrow(data_frame_subset)){
        sbc_id <- data_frame_subset$Chain_SbcID[row]
        sbc_id_fasta <- paste0(sbc_id, strrep(" ", 9-nchar(sbc_id)))
        cat(paste0(">", sbc_id_fasta, "\n", data_frame_subset$seq_write[row], "\n"))
      }
      sink()

      aln <- read.dna(paste0(outDir, "subclones_", chain, "_", cloneID, ".fasta"), format = "fasta")
      write.dna(aln, file = paste0(outDir, "deduplicated.phylip"), format = "interleaved", nbcol = 5, indent = 10)
    }


  } else if(mode == "bulk"){
    if(!"Clonotype_nReads" %in% colnames(data_frame)){stop("The provided object does not contain the expected IgScan fields. Please provide a valid IgScan output file as input.")}

    if(!cloneID %in% unique(data_frame$ClonotypeID)){stop("The specified `cloneID` is not present in the input data.")}
    data_frame_subset <- data_frame[data_frame$ClonotypeID %in% cloneID,]
    if(nrow(data_frame_subset) < 2){stop("There are not enough sequences for phylogenetic inference in the chosen clonotype.")}

    abundance_df <- rbind(c("GL", 0), data_frame_subset[, c("SubcloneID_in_Clonotype", "Subclone_nReads")])
    write.table(abundance_df[,1:2], file = paste0(outDir, "abundances.csv"), quote = F, sep = ",", row.names = F, col.names = F)

    ## Write IDMAP file
    sink(paste0(outDir, "idmap.txt"))
    for(sbc in abundance_df$SubcloneID_in_Clonotype){
      bcd <- strsplit(data_frame_subset$contig_id[data_frame_subset$SubcloneID_in_Clonotype == sbc], split = "_")[[1]][1]
      cat(paste0(sbc, ",", bcd, "\n"))
    }
    sink()

    sink(paste0(outDir, "subclones_", cloneID, ".fasta"))
    germline <- unique(data_frame_subset$Consensus_Germline)
    cat(paste0(">GL      \n", germline, "\n"))
    for(row in 1:nrow(data_frame_subset)){
      sbc_id <- data_frame_subset$SubcloneID_in_Clonotype[row]
      sbc_id_fasta <- paste0(sbc_id, strrep(" ", 9-nchar(sbc_id)))
      cat(paste0(">", sbc_id_fasta, "\n", data_frame_subset$VDJ_sequence_correctedCDR3[row], "\n"))
    }
    sink()

    aln <- read.dna(paste0(outDir, "subclones_", cloneID, ".fasta"), format = "fasta")
    write.dna(aln, file = paste0(outDir, "deduplicated.phylip"), format = "interleaved", nbcol = 5, indent = 10)
  }
  message(paste0("gctree input files have been succesfully generated in ", outDir))
}
