#' Run IgScan Annotation
#'
#' @description This function performs the IgScan annotation on Run_IgBlast_from_RawData output data.
#' It groups sequences into clonotypic clusters based on a highly customizable hierarchical
#' clustering approach. The function processes clonotypes to provide key immunogenetic
#' information, such as germline sequences, somatic hypermutation levels, CLL stereotype
#' subsets, R110 mutation, Acquired N-Glycosylation sites (AGS) and more. Then, it organizes
#' this information according to the input data type: for single-cell RNA-seq datasets, it groups heavy-light
#' chains by cell barcodes, while for bulk-NGS datasets, it reports the cumulative number of reads and
#' clonotype frequencies. In addition, the workflow is also compatible with single-cell V(D)J data
#' generated using the Mission Bio platform, provided that the IgScan input has been produced
#' using our dedicated raw-data preprocessing pipeline (see IgScan GitHub).
#'
#' @param analysis_mode Defines the mode of analysis to be performed. Options are
#'   "single" and "joint". In "single" mode, each sample is annotated independently.
#'   In "joint" mode, samples are grouped by case based on the "case_list" vector
#'   for a case-level joint annotation. Default is 'single'.
#' @param sample_labels A vector of sample labels to be annotated. Must match the
#'   sample labels of the Run_IgBlast_from_RawData output files. This parameter can be also
#'   set to 'all_samples', which will analyze all the samples in the igblast_out directory
#'   in 'single' analysis_mode. Default is 'all_samples'.
#' @param case_labels A vector of case labels matching the indexes of `sample_labels`.
#'   Required when analysis mode is set to 'joint'. Default is NULL.
#' @param input_format A string specifying the format of the input data, currently supporting:
#'   '10xBCR_fasta', '10xBCR_csv', 'ParseBCR', 'BDRhapsodyBCR', 'MiXCR',
#'   'TRUST4', 'AIRR', 'IMGT_AIRR' and 'fasta'.
#' @param data_type The type of data. Options: 'single_cell', 'bulk' or 'missionbio'. Default is 'single_cell'.
#' @param material_type The biological source of material. Options: 'dna' and 'rna'. Note that for
#'   matrial_type='rna', unproductive sequences are not expected, and will be directly removed. Default is 'rna'.
#' @param v_primer The primer sequence used for the V-region amplification.
#'   Options: 'full_length', 'fr1', 'fr2' and 'fr3'. Only required for `data_type='bulk'`. Default is
#'   'full_length'. Note that sequences with unexpected length pattern based on the chosen primer
#'   will be directly excluded from the analysis.
#' @param min_reads The minimum number of reads required to keep a sequence in the downstream analyses.
#'   Only needed for bulk NGS data. Default is 2.
#' @param hc_similarity_cutoff Hierarchical clustering distance in CDR3 identity to consider two sequences
#'   clonally related. Default is 0.2, meaning that up to 20 percent of mismatch is accepted.
#' @param hc_mode The aglomeration method to be used in the hierarchical clustering. Options are: 'single',
#'   'complete' and 'average' (UPGMA). Default is 'average'.
#' @param cdr3_mode The type of CDR3 sequence to be used for the hierarchical clustering.
#'   Options are 'nt' (nucleotides) and 'aa' (amino acids). Default is 'nt'.
#' @param cdr3_InDel_correction_mode The degree of stringency in correcting InDels in the CDR3 region. Options
#'   are 'no' (no CDR3 InDel correction, recommended for highly polyclonal samples), 'hard_filter' (high stringency
#'   but computationally slower) and 'soft_filter' (less stringent but much faster). Default is 'soft_filter'.
#' @param annotate_CLL_immGen Logical value indicating whether to annotate CLL-related immunogenetic
#'   information such as CLL stereotyped subsets and the R110 mutation. Default is FALSE.
#' @param annotate_satellite_subsets Logical value indicating whether to annotate Satellite CLL stereotyped subsets.
#' Only needed if `annotate_CLL_immGen` is set to TRUE. Default is TRUE.
#' @param annotate_ags Logical value indicating whether to annotate IGH Acquired N-Glycosylation Sites (AGS).
#' Default is FALSE.
#' @param rescue_single_chain Logical value indicating whether to hard assign cells with one chain detected to the most
#' plausible complete clonotype (see `rescue_single_chain_cells` function documentation for more details). Default is FALSE.
#' @param relaxed_rescue Logical value indicating whether to apply a relaxed clonotype rescue mode.
#' Only needed if `rescue_single_chain` is set to `TRUE`.The default (recommended) is `FALSE`, requiring exact
#' V(D)J nucleotide sequence matches to rescue single-chain cells. If set to `TRUE`, which we recommend only
#' for the analysis of Mission Bio single-cell V(D)J data, rescue is performed at the clonotype level rather
#' than by exact nucleotide sequence identity.
#' @param outputDir Path to the directory containing the previous IgScan outputs coming from the
#'   Run_IgBlast_from_RawData function. There is NO default value for this parameter.
#' @param remove_tmp Logical value indicating whether to remove the temporary files (all files
#'   except the annotation_results directory). Default is TRUE.
#' @param threads An integer specifying the number of threads for processing samples.
#'   Default is 1.
#'
#' @return A list containing an IgScan annotation dataframe for every sample analyzed, which are also
#' saved in the output directory. These dataframes can be directly passed to further IgScan functions.
#'
#' @export
#'
#' @import dplyr
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom purrr list_flatten
#' @importFrom parallel mclapply
#' @import stringr
#' @importFrom data.table fread fwrite
#' @importFrom Seurat CreateSeuratObject
#'
#' @examples
#' \dontrun{
#' # Example 1: XXXXX
#' result <- Run_IgScan_Annotation(analysis_mode = "joint",
#'    sample_labels = c("Sample1", "Sample2"),
#'    case_labels = c("Case1", "Case1"),
#'    input_format = "fasta",
#'    material_type = "dna",
#'    v_primer = "full_length",
#'    data_type = "bulk",
#'    min_reads = 2,
#'    remove_tmp = TRUE,
#'    outputDir = "path/output/dir/",
#'    hc_similarity_cutoff = 0.2,
#'    hc_mode = "average",
#'    cdr3_mode = "nt",
#'    cdr3_InDel_correction_mode = "soft_filter",
#'    annotate_CLL_immGen = TRUE,
#'    annotate_satellite_subsets = FALSE,
#'    annotate_ags = TRUE,
#'    threads = 4)
#'
#' # Example 2: XXXXX
#' result <- Run_IgScan_Annotation(
#'    analysis_mode = "joint",
#'    sample_labels = c("Sample1", "Sample2"),
#'    case_labels = c("Case1", "Case1"),
#'    input_format = "fasta",
#'    material_type = "dna",
#'    v_primer = "full_length",
#'    data_type = "bulk",
#'    min_reads = 2,
#'    remove_tmp = TRUE,
#'    outputDir = "path/output/dir/",
#'    hc_similarity_cutoff = 0.2,
#'    hc_mode = "average",
#'    cdr3_mode = "nt",
#'    cdr3_InDel_correction_mode = "soft_filter",
#'    annotate_CLL_immGen = TRUE,
#'    annotate_satellite_subsets = FALSE,
#'    annotate_ags = TRUE,
#'    threads = 4)
#' }
#'
Run_IgScan_Annotation <- function(sample_labels = "all_samples", case_labels = NULL, input_format, outputDir, analysis_mode = "single", material_type = "rna", v_primer = "full_length", data_type = "single_cell", min_reads = 2, remove_tmp = TRUE, hc_similarity_cutoff = 0.2, hc_mode = "average", cdr3_mode = "nt", cdr3_InDel_correction_mode = "soft_filter", annotate_CLL_immGen = FALSE, annotate_satellite_subsets = TRUE, annotate_ags = FALSE, rescue_single_chain = FALSE, relaxed_rescue = FALSE, summary_file = NULL, threads = 1){


  ## First checks
  if(is.null(sample_labels)){
    sample_labels <- "all_samples"
  } else{
    analysis_mode <- tolower(analysis_mode)
    if(!analysis_mode %in% c("single","joint")){stop("Invalid value for 'analysis_mode'. It should be either 'single' or 'joint'.")}
    if(is.null(case_labels) & analysis_mode == "joint"){stop("'case_labels' is required when 'analyisis_mode' is set to 'joint'. Please, set a valid combination of parameters!")}
    if(analysis_mode == "joint" & length(case_labels) != length(sample_labels)){stop("Invalid list of cases. There is different number of samples and case indexes!")}
  }

  if(!is.logical(remove_tmp)){stop("Invalid value for 'remove_tmp'. It should be either TRUE or FALSE.")}
  if((!is.numeric(hc_similarity_cutoff) | hc_similarity_cutoff > 1) & hc_similarity_cutoff != "automatic"){stop("Invalid value for 'hc_similarity_cutoff'. It should be either a float between 0 and 1 or 'automatic'.")}
  hc_mode <- tolower(hc_mode)
  if(!hc_mode %in% c("average", "complete", "single")){stop("Invalid value for 'hc_mode'. It should be either 'average', 'complete' or 'single'.")}
  cdr3_mode <- tolower(cdr3_mode)
  if(!cdr3_mode %in% c("nt", "aa")){stop("Invalid value for 'cdr3_mode'. It should be either 'nt' or 'aa'.")}
  cdr3_InDel_correction_mode <- tolower(cdr3_InDel_correction_mode)
  if(!cdr3_InDel_correction_mode %in% c("soft_filter", "hard_filter", "no")){stop("Invalid value for 'cdr3_mode'. It should be either 'soft_filter', 'hard_filter' or 'no'.")}
  if(!dir.exists(outputDir)){stop("outputDir does not exist. Please, provide a valid IgScan output directory!")}
  if(!endsWith(outputDir, "/")){outputDir <- paste0(outputDir, "/")}
  input_format <- tolower(input_format)
  material_type <- tolower(material_type)
  if(!material_type %in% c("dna","rna")){stop("Invalid value for 'material_type'. It should be either 'dna' or 'rna'.")}
  data_type <- tolower(data_type)
  if(!data_type %in% c("single_cell","bulk","missionbio")){stop("Invalid value for 'data_type'. It should be either 'single_cell' 'missionbio', or 'bulk'.")}
  if(data_type == "single_cell"){
    v_primer <- "full_length"
  } else if(data_type == "missionbio"){
    v_primer <- "missionbio"
  } else{
    v_primer <- tolower(v_primer)
    if(!v_primer %in% c("full_length","fr1", "fr2", "fr3")){stop("Invalid value for 'v_primer'. It should be either 'full_length', 'fr1', 'fr2' or 'fr3'.")}
  }
  if(!is.logical(annotate_CLL_immGen)){stop("Invalid value for 'annotate_CLL_immGen'. It should be either TRUE or FALSE.")}
  if(!is.logical(annotate_satellite_subsets)){stop("Invalid value for 'annotate_satellite_subsets'. It should be either TRUE or FALSE.")}
  if(!is.logical(annotate_ags)){stop("Invalid value for 'annotate_ags'. It should be either TRUE or FALSE.")}
  if(!is.logical(rescue_single_chain)){stop("Invalid value for 'rescue_single_chain'. It should be either TRUE or FALSE.")}
  if(!is.logical(relaxed_rescue)){stop("Invalid value for 'relaxed_rescue'. It should be either TRUE or FALSE.")}

  ## Create the annotation result directory
  system(paste0("mkdir ",outputDir,"annotation_results/"))

  if(analysis_mode == "single"){
    if(length(sample_labels) == 1 && sample_labels == "all_samples"){
      files_to_annotate <- list.files(paste0(outputDir,"igblast_outs"), full.names = T, recursive = F, pattern = "_igblast_out\\.tsv$")
    } else{
      files_to_annotate <- list.files(paste0(outputDir,"igblast_outs"), full.names = T, recursive = F, pattern = paste0("^(",paste(sample_labels, collapse = "|"), ")_igblast_out\\.tsv"))
    }

  } else if (analysis_mode == "joint"){
    system(paste0("mkdir ",outputDir,"per_case_merge/"))
    for(case in unique(case_labels)){
      merge_df <- data.frame()
      for(sample in sample_labels[which(case_labels == case)]){
        sample_df <- fread(paste0(outputDir,"igblast_outs/",sample,"_igblast_out.tsv"), header = T, sep = "\t", stringsAsFactors = F, data.table = F)
        sample_df$sample <- sample
        merge_df <- rbind(merge_df, sample_df)
      }
      fwrite(merge_df, paste0(outputDir,"per_case_merge/",case,"_merge_df.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)
    }
    files_to_annotate <- paste0(outputDir, "per_case_merge/", unique(case_labels), "_merge_df.tsv")
  }

  all_out_list <- .run_Core_IgScan_annotation(files_to_annotate, outputDir, analysis_mode, material_type, v_primer, data_type, min_reads, sample_labels, hc_similarity_cutoff, hc_mode, cdr3_mode, cdr3_InDel_correction_mode, annotate_CLL_immGen, annotate_satellite_subsets, annotate_ags, rescue_single_chain, relaxed_rescue, summary_file, threads, input_format)
  flat_out_list <- list_flatten(all_out_list)
  final_list <- .name_final_list(flat_out_list)

  if(remove_tmp == TRUE){
    system(paste0("rm -r ",outputDir, "fasta_inputs/ ", outputDir, "igblast_outs/"))
    if(analysis_mode == "joint"){system(paste0("rm -r ",outputDir, "per_case_merge/"))}
  }

  return(final_list)
}

## Function to parallelize the Core_IgScan_annotation function
.run_Core_IgScan_annotation <- function(files_to_annotate, outputDir, analysis_mode, material_type, v_primer, data_type, min_reads, sample_labels, hc_similarity_cutoff, hc_mode, cdr3_mode, cdr3_InDel_correction_mode, annotate_CLL_immGen, annotate_satellite_subsets, annotate_ags, rescue_single_chain, relaxed_rescue, summary_file, threads, input_format){
  mclapply(files_to_annotate, function(i) {
    name <- gsub("_igblast_out.tsv|_merge_df.tsv", "", basename(i))
    input_df <- fread(i, header = T, sep = "\t", stringsAsFactors = F, data.table = F)
    if(endsWith(i, "_igblast_out.tsv")){input_df$sample <- name}
    output_df <- .Core_IgScan_annotation(input_df, outputDir, name, analysis_mode, material_type, v_primer, data_type, min_reads, sample_labels, hc_similarity_cutoff, hc_mode, cdr3_mode, cdr3_InDel_correction_mode, annotate_CLL_immGen, annotate_satellite_subsets, annotate_ags, rescue_single_chain, relaxed_rescue, summary_file, input_format, threads)
    return(output_df)
  }, mc.cores = threads)
}

## General IgScan annotation function
.Core_IgScan_annotation <- function(input_df, outputDir, name, analysis_mode, material_type, v_primer, data_type, min_reads, sample_labels, hc_similarity_cutoff, hc_mode, cdr3_mode, cdr3_InDel_correction_mode, annotate_CLL_immGen, annotate_satellite_subsets, annotate_ags, rescue_single_chain, relaxed_rescue, summary_file, input_format, threads){

  total_tasks <- ifelse(data_type == "bulk", 5, 9)
  completed_tasks <- 0

  ## Start the summary file for each run
  if(is.null(summary_file)){
    summary_file <- paste0(outputDir,"annotation_results/", ifelse(analysis_mode == "joint", "case_", "sample_"), name, "_reporting_summary.log")
    write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Preparing IgBlast outputs for IG Clonotype annotation."), file = summary_file)
  } else{
    write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Preparing IgBlast outputs for IG Clonotype annotation."), file = summary_file, append = T)
  }
  message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Preparing IgBlast outputs for IG Clonotype annotation."))

  tidy_dataset <- .tidy_input_df(input_df, material_type, v_primer)
  tidy_dataset <- .correct_InDels(tidy_dataset)

  tidy_dataset$junction <- apply(tidy_dataset, MARGIN = 1, FUN = .update_junction)
  for(jun in unique(tidy_dataset$junction)){
    tidy_dataset$junction_aa[tidy_dataset$junction == jun] <- .translate_sequence(jun)
  }
  tidy_dataset$junction_aa_length <- nchar(tidy_dataset$junction_aa)

  if(data_type == "bulk"){
    read_names <- tidy_dataset$sequence_id
    read_count <- sapply(read_names, function(x) strsplit(x, split = "=")[[1]][2])
    if(all(grepl("=", read_names)) && all(grepl("^[0-9]+$", read_count))){
      tidy_dataset$n_reads <- as.numeric(read_count)
    } else{
      warning("Read counts were not found in fasta read tags (expected tag format ==> readName_n=reacCount).\nAssuming that each sequence has 1 read count.")
      tidy_dataset$n_reads <- 1
    }

    tidy_dataset$VDJseq_indels <- paste0(tidy_dataset$VDJseq, "-", tidy_dataset$indel)
    dict_unique_seq <- aggregate(x = tidy_dataset$n_reads, by = list(tidy_dataset$VDJseq_indels), FUN = sum)
    dict_unique_seq <- dict_unique_seq[order(dict_unique_seq$x, decreasing = T), ]
    dict_unique_seq$Unique_SequenceID <- paste0("UniqueSeq_", 1:nrow(dict_unique_seq))
    tidy_dataset$Unique_SequenceID <- dict_unique_seq$Unique_SequenceID[match(tidy_dataset$VDJseq_indels, dict_unique_seq$Group.1)]

    ## Calculate the total number of reads from every unique sequence
    tidy_dataset <- tidy_dataset %>%
      group_by(Unique_SequenceID) %>%
      mutate(total_reads_unique_seq = sum(n_reads))

    rows_to_remove <- which(tidy_dataset$total_reads_unique_seq < min_reads)
    if(length(rows_to_remove) > 0){
      tidy_dataset <- tidy_dataset[-rows_to_remove,]
    }
    tidy_dataset <- tidy_dataset[!duplicated(tidy_dataset$Unique_SequenceID),]
  }

  ## Label clonotype labels (V gene + CDR3 amino acid sequence)
  list_unique_Vs <- unique(sapply(tidy_dataset$VDJ, function(x) strsplit(x, split = "/")[[1]][1]))
  list_unique_Vs <- unname(sapply(list_unique_Vs, function(x) str_replace_all(x, "\\*\\d+", "")))
  list_unique_Vs <- unique(unname(sapply(list_unique_Vs, function(x) paste0(sort(unique(strsplit(x, ",")[[1]])), collapse = ","))))
  list_unique_Vs <- list_unique_Vs[order(nchar(list_unique_Vs), decreasing = TRUE)]

  final_IGHV <- sapply(.combine_IGHV_genes(strsplit(list_unique_Vs, ",")), function(group) {
    paste(sort(group), collapse = ",")
  })

  v_gene <- unname(sapply(tidy_dataset$VDJ, function(x) strsplit(x, split = "\\*")[[1]][1]))
  u_v_gene <- unique(v_gene)
  u_v_gene_match <- unname(unlist(sapply(u_v_gene, function(v) { final_IGHV[sapply(final_IGHV, function(x) v %in% strsplit(x, ",")[[1]])] })))
  v_gene <- u_v_gene_match[match(v_gene, u_v_gene)]

  if(cdr3_mode == "aa"){
    tidy_dataset$clonotypeLabel <- paste(v_gene, sapply(tidy_dataset$junction_aa, function(x) paste0(strsplit(x, "")[[1]][2:(nchar(x)-1)], collapse = "")), sep = "_")
  }else if(cdr3_mode == "nt"){
    tidy_dataset$clonotypeLabel <- paste(v_gene, sapply(tidy_dataset$junction, function(x) paste0(strsplit(x, "")[[1]][4:(nchar(x)-3)], collapse = "")), sep = "_")
  }

  ## Removing incomplete sequences (not if missionbio)
  if(data_type != "missionbio"){
    removed_rows <- 0
    write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Iterate ", length(unique(tidy_dataset$clonotypeLabel)), " clonotypes to remove incomplete sequences."), file = summary_file, append = T)
    message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Iterate ", length(unique(tidy_dataset$clonotypeLabel)), " clonotypes to remove incomplete sequences."))
    rows_to_remove <- as.vector(unlist(sapply(unique(tidy_dataset$clonotypeLabel), function(clonotype){
      return(which(tidy_dataset$clonotypeLabel == clonotype &
                     nchar(tidy_dataset$VDJseq) != as.numeric(names(which.max(table(nchar(tidy_dataset$VDJseq[tidy_dataset$clonotypeLabel == clonotype])))))))
    })))
    if(length(rows_to_remove) > 0){
      write(x = paste0("Sequence with ID ", tidy_dataset[rows_to_remove,"sequence_id"], " was removed from sample ", tidy_dataset[rows_to_remove,"sample"], " since VDJseq lenght was in disagreement with most recurrent sequence."), file = summary_file, append = T)
      message(paste0("Sequence with ID ", tidy_dataset[rows_to_remove,"sequence_id"], " was removed from sample ", tidy_dataset[rows_to_remove,"sample"], " since VDJseq lenght was in disagreement with most recurrent sequence.\n"))
      tidy_dataset <- tidy_dataset[-rows_to_remove,]
      removed_rows <- length(rows_to_remove)
    }
    write(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": ", removed_rows, " sequences have been removed from this sample."), file = summary_file, append = T)
    message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": ", removed_rows, " sequences have been removed from this sample."))
  }

  write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Starting IG Clonotype annotation... ",completed_tasks,"/",total_tasks," tasks completed."), file = summary_file, append = T)
  message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Starting IG Clonotype annotation... ",completed_tasks,"/",total_tasks," tasks completed."))
  completed_tasks <- completed_tasks + 1

  assign_Clonotypes_out <- .assign_Clonotypes(tidy_dataset, hc_similarity_cutoff, hc_mode, cdr3_mode, cdr3_InDel_correction_mode, summary_file, analysis_mode, name, data_type, total_tasks, completed_tasks)
  tidy_dataset <- assign_Clonotypes_out[[1]]
  completed_tasks <- assign_Clonotypes_out[[2]]

  if(data_type %in% c("single_cell", "missionbio")){
    ## Combine Heavy-Chain + Light-Chain annotations
    write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Combining clonotypeIDs by cell barcodeID... ",completed_tasks,"/",total_tasks," tasks completed."), file = summary_file, append = T)
    message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Combining clonotypeIDs by cell barcodeID... ",completed_tasks,"/",total_tasks," tasks completed."))
    completed_tasks <- completed_tasks + 1

    tidy_dataset$barcode <- matrix(unlist(strsplit(tidy_dataset$sequence_id, split = "_")), ncol=3, byrow = T)[,1]
    tidy_dataset$barcode <- paste0(tidy_dataset$barcode, "_", tidy_dataset$sample)

    ## Correct cells with repeated clonotypes of the same chain
    rm_rep_chain <- unlist(sapply(unique(tidy_dataset$barcode), function(bc) tidy_dataset$sequence_id[tidy_dataset$barcode == bc][duplicated(tidy_dataset$clonotypeID[tidy_dataset$barcode == bc])]))
    tidy_dataset <- tidy_dataset[!tidy_dataset$sequence_id %in% rm_rep_chain,]

    tidy_dataset <- .combine_clonotypeID_by_chains(tidy_dataset)

    ## Adjust IGH - IGK/IGL rearrangements:
    write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Assessing BCR completeness at single cell level... ",completed_tasks,"/",total_tasks," tasks completed."), file = summary_file, append = T)
    message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Assessing BCR completeness at single cell level... ",completed_tasks,"/",total_tasks," tasks completed."))
    completed_tasks <- completed_tasks + 1

    tidy_dataset$completeBCR <- "Yes"
    for(x in 1:nrow(tidy_dataset)){
      chains <- as.character(sapply(unique(strsplit(tidy_dataset$merge_clonotypeID[x], split = "-")[[1]]), function(x){strsplit(x, split = "\\.")[[1]]})[1,])
      if(length(strsplit(tidy_dataset$merge_clonotypeID[x], split = "-")[[1]]) == 1){
        tidy_dataset$completeBCR[x] <- "Single_chain_1"
      }else if((length(strsplit(tidy_dataset$merge_clonotypeID[x], split = "-")[[1]]) == 2) & !(("IGH" %in% chains) & ("IGK" %in% chains | "IGL" %in% chains))){
        tidy_dataset$completeBCR[x] <- "Single_chain_2"
      }else if((length(strsplit(tidy_dataset$merge_clonotypeID[x], split = "-")[[1]]) == 3 | length(strsplit(tidy_dataset$merge_clonotypeID[x], split = "-")[[1]]) == 4) & !(("IGH" %in% chains) & ("IGK" %in% chains | "IGL" %in% chains))){
        tidy_dataset$completeBCR[x] <- "Not_supported"
      }else if(length(strsplit(tidy_dataset$merge_clonotypeID[x], "-")[[1]]) == 4 && (sum(chains == "IGH") > 2 || sum(chains %in% c("IGK","IGL")) > 2)){
        tidy_dataset$completeBCR[x] <- "Not_supported"
      } else if((length(strsplit(tidy_dataset$merge_clonotypeID[x], split = "-")[[1]]) > 4)){
        tidy_dataset$completeBCR[x] <- "Not_supported"
      }
    }

    if(data_type == "missionbio"){
      tidy_dataset$completeBCR[sapply(tidy_dataset$merge_clonotypeID, function(x) sum(grepl("IGH\\.", strsplit(x, "-")[[1]])) > 1) & tidy_dataset$completeBCR != "Not_supported"] <- "Potential_MB_doublet"
    }

    tidy_dataset <- .correct_completeBCR(tidy_dataset, data_type)

    ## Correct clonotypeID (chain-specific clonotypeID) considering IG clonotype (ie combination of chains)
    write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Correcting clonotype labels... ",completed_tasks,"/",total_tasks," tasks completed."), file = summary_file, append = T)
    message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Correcting clonotype labels... ",completed_tasks,"/",total_tasks," tasks completed."))
    completed_tasks <- completed_tasks + 1

    tidy_dataset$clonotypeID <- .correct_Clonotype_labels(tidy_dataset)

    ## Re-calculate CV and subclones based on new clonotypeIDs
    write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Recalculating clonotype and subclone abundances... ",completed_tasks,"/",total_tasks," tasks completed."), file = summary_file, append = T)
    message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Recalculating clonotype and subclone abundances... ",completed_tasks,"/",total_tasks," tasks completed."))
    completed_tasks <- completed_tasks + 1

    tidy_dataset <- .make_CV_and_subclones_single_cell(tidy_dataset)
    tidy_dataset <- .combine_clonotypeID_by_chains(tidy_dataset)
  }

  ## Annotate immunogenetic data such as Ig V gene identity, stereotype subsets, R110, etc
  write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Annotating immunogenetic data... ",completed_tasks,"/",total_tasks," tasks completed."), file = summary_file, append = T)
  message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Annotating immunogenetic data... ",completed_tasks,"/",total_tasks," tasks completed."))
  completed_tasks <- completed_tasks + 1

  tidy_dataset <- .annotate_Immunogenetics(tidy_dataset, annotate_CLL_immGen, annotate_satellite_subsets, annotate_ags, data_type)

  if(data_type %in% c("single_cell", "missionbio")){
    tidy_dataset <- .combine_IG_metadata_by_chain(tidy_dataset, annotate_CLL_immGen)
  }

  ## Making a IGH-IGK/IGL clonotypeVariant ID for further analysis
  write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Preparing table to write... ",completed_tasks,"/",total_tasks," tasks completed."), file = summary_file, append = T)
  message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Preparing table to write... ",completed_tasks,"/",total_tasks," tasks completed."))

  if(data_type %in% c("single_cell", "missionbio")){

    dict_CV <- aggregate(x = tidy_dataset$merge_clonotypeVariantID[!is.na(tidy_dataset$igclonotypeID)], by = list(tidy_dataset$merge_clonotypeVariantID[!is.na(tidy_dataset$igclonotypeID)]), FUN = length)
    dict_CV <- merge(dict_CV, tidy_dataset[!duplicated(tidy_dataset$merge_clonotypeVariantID), c("merge_clonotypeVariantID","igclonotypeID")], by.x = "Group.1", by.y = "merge_clonotypeVariantID", all.x = T, all.y = F)
    dict_CV <- dict_CV[order(dict_CV$x, decreasing = T), ]
    for(clonotype in unique(dict_CV$igclonotypeID)){
      dict_CV$igClonotypeVariantID[dict_CV$igclonotypeID == clonotype] <- paste0(clonotype,".CV",1:nrow(dict_CV[dict_CV$igclonotypeID == clonotype,]))
    }
    tidy_dataset <- merge(tidy_dataset, dict_CV[,c(1,4)], by.x = "merge_clonotypeVariantID", by.y = "Group.1", all.x = TRUE)

    ## Making a IGH-IGK/IGL subclone ID for further analysis (per clonotype variant and per clonotype)
    dict_sc2 <- aggregate(x = tidy_dataset$merge_subcloneID[!is.na(tidy_dataset$igclonotypeID)], by = list(tidy_dataset$merge_subcloneID[!is.na(tidy_dataset$igclonotypeID)]), FUN = length)
    dict_sc2 <- merge(dict_sc2, tidy_dataset[!duplicated(tidy_dataset$merge_subcloneID), c("merge_subcloneID","igClonotypeVariantID")], by.x = "Group.1", by.y = "merge_subcloneID", all.x = T, all.y = F)
    dict_sc2 <- dict_sc2[order(dict_sc2$x, decreasing = T), ]
    for(CV in unique(dict_sc2$igClonotypeVariantID)){
      dict_sc2$igsubcloneID_in_clonotypeVariant[dict_sc2$igClonotypeVariantID == CV] <- paste0(CV,".S",1:nrow(dict_sc2[dict_sc2$igClonotypeVariantID == CV,]))
    }
    tidy_dataset <- merge(tidy_dataset, dict_sc2[,c(1,4)], by.x = "merge_subcloneID", by.y = "Group.1", all.x = TRUE)

    dict_sc3 <- aggregate(x = tidy_dataset$merge_subcloneID[!is.na(tidy_dataset$igclonotypeID)], by = list(tidy_dataset$merge_subcloneID[!is.na(tidy_dataset$igclonotypeID)]), FUN = length)
    dict_sc3 <- merge(dict_sc3, tidy_dataset[!duplicated(tidy_dataset$merge_subcloneID), c("merge_subcloneID","igclonotypeID")], by.x = "Group.1", by.y = "merge_subcloneID", all.x = T, all.y = F)
    dict_sc3 <- dict_sc3[order(dict_sc3$x, decreasing = T), ]
    for(clonotype in unique(dict_sc3$igclonotypeID)){
      dict_sc3$igsubcloneID_in_clonotype[dict_sc3$igclonotypeID == clonotype] <- paste0(clonotype,".",1:nrow(dict_sc3[dict_sc3$igclonotypeID == clonotype,]))
    }
    tidy_dataset <- merge(tidy_dataset, dict_sc3[,c(1,4)], by.x = "merge_subcloneID", by.y = "Group.1", all.x = TRUE)

  } else if(data_type == "bulk"){

    tidy_dataset <- tidy_dataset %>%
      group_by(clonotypeID) %>%
      mutate(clonotypeID_nreads = sum(total_reads_unique_seq))

    tidy_dataset <- tidy_dataset %>%
      group_by(clonotypeVariantID_in_Cltp) %>%
      mutate(clonotypeVariant_nreads = sum(total_reads_unique_seq))

    tidy_dataset$clonotypeID_freq <- (tidy_dataset$clonotypeID_nreads / sum(tidy_dataset$total_reads_unique_seq))*100

    tidy_dataset <- tidy_dataset %>%
      group_by(clonotypeID) %>%
      mutate(clonotypeVariant_freq = clonotypeVariant_nreads/sum(total_reads_unique_seq)*100,
             subclone_freq = total_reads_unique_seq/sum(total_reads_unique_seq)*100)

    tidy_dataset <- tidy_dataset %>%
      group_by(clonotypeVariantID_in_Cltp) %>%
      mutate(subclone_freq_in_CV = total_reads_unique_seq/sum(total_reads_unique_seq)*100)

    tidy_dataset <- as.data.frame(tidy_dataset[order(tidy_dataset$total_reads_unique_seq, decreasing = T),])
  }

  contig_data <- c("sequence_id", "sequence", "raw_VDJseq", "VDJseq", "germline_alignment", "VDJseq_correctedCDR3",
                   "VDJseq_correctedCDR3_aa", "CorrectClonotypes_Consensus_Germline",
                   "CorrectClonotypes_Consensus_Germline_aa", "VDJ", "c_call", "productive", "junction_aa",
                   "junction_aa_length", "v_identity", "arch", "len_no_CDR3", "indel",
                   ifelse(annotate_CLL_immGen, "IGHsubset", NA))

  if(data_type == "bulk"){
    contig_ID_data <- c("clonotypeID", "Clonotype_ConsensusCDR3", ifelse(annotate_CLL_immGen, "Clonotype_Subset", NA),
                        "clonotypeID_nreads", "clonotypeID_freq", "clonotypeVariantID_in_Cltp", "clonotypeVariant_nreads", "clonotypeVariant_freq",
                        "subcloneID_in_CV", "total_reads_unique_seq", "subclone_freq_in_CV", "subcloneID_in_Clt", "subclone_freq")
    perCell_data <- NA

  } else if(data_type %in% c("single_cell", "missionbio")){
    contig_ID_data <- c("clonotypeID", "Clonotype_ConsensusCDR3", ifelse(annotate_CLL_immGen, "Clonotype_Subset", NA),
                        "clonotypeVariantID_in_Cltp", "Unique_SequenceID")
    perCell_data <- c("completeBCR", "igclonotypeID", "merge_clonotypeID", "igClonotypeVariantID", "merge_clonotypeVariantID",
                      "igsubcloneID_in_clonotypeVariant", "igsubcloneID_in_clonotype", "merge_subcloneID",
                      "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", ifelse(annotate_CLL_immGen, "V9", NA))
  }

  write_fields <- c(contig_data, contig_ID_data, perCell_data, "sample")
  table_to_write <- tidy_dataset[,write_fields[!is.na(write_fields)]]

  ## Load colnames dictionary from IgScan inst/
  colnames_dictionary_path <- system.file("colnames_dictionary.RData", package = "IgScan", mustWork = T)
  load(colnames_dictionary_path)

  colnames(table_to_write) <- colnames_dictionary$FinalNames[match(colnames(table_to_write), colnames_dictionary$RawNames)]
  table_to_write[table_to_write == ""] <- NA

  if(data_type %in% c("single_cell", "missionbio")){

    if(rescue_single_chain){
      table_to_write <- rescue_single_chain_cells(single_cell_object = table_to_write, group_col = "SampleID", threads = threads, relaxed_rescue)
    }

    if(data_type == "missionbio"){
      barcodes <- sapply(table_to_write$contig_id[table_to_write$completeBCR != "Not_supported"], function(x) strsplit(x, "_")[[1]][1])
      n_cells <- length(unique(barcodes))
      fake_counts <- matrix(data = 0, nrow = 2, ncol = n_cells, dimnames = list(c("fake-gene1", "fake-gene2"), unique(barcodes)))
      fake_seurat <- CreateSeuratObject(counts = as(fake_counts, "dgCMatrix"), assay = "RNA", project = "fake_seurat")
      fake_seurat@meta.data$orig.ident <- table_to_write$SampleID[match(Cells(fake_seurat), barcodes)]
      fake_seurat_igscan <- combine_IgScan_Seurat(igscan_out = table_to_write[table_to_write$completeBCR != "Not_supported",], seurat_object = fake_seurat, seurat_sample_col = "orig.ident", igscan_sample_col = "SampleID")
      table_to_write_cell_mb <- fake_seurat_igscan@meta.data[,!colnames(fake_seurat_igscan@meta.data) %in% c("nCount_RNA", "nFeature_RNA")]

      if(any(table_to_write$completeBCR == "Not_supported")){
        barcodes_notsup <- unique(sapply(table_to_write$contig_id[table_to_write$completeBCR == "Not_supported"], function(x) strsplit(x, "_")[[1]][1]))
        tmp_df <- as.data.frame(matrix(NA, nrow = length(barcodes_notsup), ncol = ncol(table_to_write_cell_mb), dimnames = list(barcodes_notsup, colnames(table_to_write_cell_mb))))
        tmp_df$orig.ident <- unique(table_to_write_cell_mb$orig.ident)
        tmp_df$completeBCR <- "Not_supported"
        table_to_write_cell_mb <- rbind(table_to_write_cell_mb, tmp_df)
      }
    }

    if(analysis_mode == "single"){
      output_df <- table_to_write
      fwrite(output_df, file = paste0(outputDir,"annotation_results/",name,"_IGannotation.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
      if(data_type == "missionbio"){fwrite(table_to_write_cell_mb, file = paste0(outputDir,"annotation_results/",name,"_IGannotation_per_MissionBio_barcode.tsv"), sep = "\t", quote = F, col.names = T, row.names = T)}

    } else if(analysis_mode == "joint"){
      output_df <- list()
      for(sample in unique(table_to_write$SampleID)){
        wt <- table_to_write[table_to_write$SampleID == sample,]
        fwrite(wt, file = paste0(outputDir,"annotation_results/case_",name,"_sample_",sample,"_IGannotation.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
        output_df <- c(output_df, list(wt))
        if(data_type == "missionbio"){
          wt_cell_mb <- table_to_write_cell_mb[table_to_write_cell_mb$orig.ident == sample,]
          fwrite(wt_cell_mb, file = paste0(outputDir,"annotation_results/case_",name,"_sample_",sample,"_IGannotation_per_MissionBio_barcode.tsv"), sep = "\t", quote = F, col.names = T, row.names = T)
        }
      }
    }

  } else if(data_type == "bulk"){

    if(analysis_mode == "single"){
      output_df <- table_to_write
      fwrite(output_df, file = paste0(outputDir,"annotation_results/",name,"_IGannotation.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)

    } else if(analysis_mode == "joint"){
      output_df <- list()
      for(sample in unique(table_to_write$SampleID)){
        wt <- table_to_write[table_to_write$SampleID == sample,]
        fwrite(wt, file = paste0(outputDir,"annotation_results/case_",name,"_sample_",sample,"_IGannotation.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
        output_df <- c(output_df, list(wt))
      }
    }
  }

  write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": IgScan annotation succesfully completed!"), file = summary_file, append = T)
  message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": IgScan annotation succesfully completed!"))
  return(output_df)
}

## Core functions:
.tidy_input_df <- function(input_df, material_type, v_primer){

  if(!"c_call" %in% colnames(input_df)){input_df$c_call <- NA}

  fields_igblast <- c("sequence_id", "productive", "complete_vdj", "v_call", "d_call", "j_call", "c_call", "fwr1", "cdr1",
                      "fwr2", "cdr2", "fwr3", "cdr3", "fwr4", "v_identity", "junction", "junction_aa",
                      "junction_aa_length", "sequence", "sequence_alignment", "germline_alignment","sample")

  sub_data <- input_df[,fields_igblast]

  if(material_type == "rna" | v_primer == "missionbio"){
    sub_data <- sub_data[sub_data$productive %in% c("T", T),]
  }

  if(v_primer == "full_length"){
    sub_data <- sub_data[sub_data$fwr1 != "" & sub_data$cdr1 != "" & sub_data$fwr2 != "" & sub_data$cdr2 != "" & sub_data$fwr3 != "" & sub_data$cdr3 != "" & sub_data$fwr4 != "" & sub_data$complete_vdj %in% c("T", T),]
  } else if(v_primer == "fr1"){
    sub_data <- sub_data[sub_data$fwr1 != "" & sub_data$cdr1 != "" & sub_data$fwr2 != "" & sub_data$cdr2 != "" & sub_data$fwr3 != "" & sub_data$cdr3 != "" & sub_data$fwr4 != "" & sub_data$complete_vdj %in% c("F", F),]
  } else if(v_primer == "fr2"){
    sub_data <- sub_data[sub_data$fwr1 == "" & sub_data$cdr1 == "" & sub_data$fwr2 != "" & sub_data$cdr2 != "" & sub_data$fwr3 != "" & sub_data$cdr3 != "" & sub_data$fwr4 != "" & sub_data$complete_vdj %in% c("F", F),]
  } else if(v_primer == "fr3"){
    sub_data <- sub_data[sub_data$fwr1 == "" & sub_data$cdr1 == "" & sub_data$fwr2 == "" & sub_data$cdr2 == "" & sub_data$fwr3 != "" & sub_data$cdr3 != "" & sub_data$fwr4 != "" & sub_data$complete_vdj %in% c("F", F),]
  } else if(v_primer == "missionbio"){
    sub_data <- sub_data[sub_data$v_call != "" | sub_data$cdr3 != "",]
  }

  ## Extracting the V, D and J annotated genes and creating the VDJ column
  v <- gsub("/", "-", as.character(sub_data$v_call))
  d <- gsub("/", "-", as.character(sub_data$d_call))
  j <- gsub("/", "-", as.character(sub_data$j_call))

  sub_data$VDJ <- paste(v,"/",d,"/",j, sep = "")
  sub_data$VDJ <- gsub("//", "/", sub_data$VDJ)

  sub_data$raw_VDJseq <- gsub("-", "", sub_data$sequence_alignment)
  sub_data$VDJseq <- sub_data$sequence_alignment
  sub_data$len_no_CDR3 <- nchar(paste0(sub_data$fwr1, sub_data$cdr1, sub_data$fwr2, sub_data$cdr2, sub_data$fwr3))
  sub_data$len_yes_CDR3 <- nchar(paste0(sub_data$fwr1, sub_data$cdr1, sub_data$fwr2, sub_data$cdr2, sub_data$fwr3, sub_data$cdr3))
  sub_data$arch <- paste0(nchar(sub_data$fwr1),"-", nchar(sub_data$cdr1),"-", nchar(sub_data$fwr2),"-", nchar(sub_data$cdr2),"-",
                          nchar(sub_data$fwr3),"-", nchar(sub_data$cdr3),"-",
                          nchar(sub_data$VDJseq)-sub_data$len_yes_CDR3)

  tidy_fields <- c("sequence_id", "productive", "VDJ", "c_call", "raw_VDJseq", "VDJseq","v_identity", "junction", "junction_aa",
                   "junction_aa_length","sequence","germline_alignment", "len_no_CDR3","arch","sample")

  tidy_dataset <- subset(sub_data, select = tidy_fields)
  tidy_dataset$productive <- ifelse(tidy_dataset$productive %in% c(T,"T"), "productive", "unproductive")
  tidy_dataset$clonotypeLabel <- NA
  tidy_dataset <- tidy_dataset[order(tidy_dataset$sequence_id),]

  return(tidy_dataset)
}
.correct_InDels <- function(tidy_dataset){

  tidy_dataset$indel <- NA

  for(row in (1:nrow(tidy_dataset))[grepl("[ACTG]-+[ACTG]", tidy_dataset$VDJseq) | grepl("[ACTG]-+[ACTG]", tidy_dataset$germline_alignment)]){

    ## Find patterns of gaps in IgBlast output indicating the presence of InDels
    deletion_list <- unlist(apply(str_locate_all(tidy_dataset$VDJseq[row], "[ACTG]-+[ACTG]")[[1]], 1, list))
    insertion_list <- unlist(apply(str_locate_all(tidy_dataset$germline_alignment[row], "[ACTG]-+[ACTG]")[[1]], 1, list))

    ## Correct Deletions
    if(length(deletion_list) > 0){
      removed_nt <- 0
      for(start in seq(1, length(deletion_list), 2)){
        somatic <- str_split_fixed(tidy_dataset$VDJseq[row], pattern = "", nchar(tidy_dataset$VDJseq[row]))
        germline <- str_split_fixed(tidy_dataset$germline_alignment[row], pattern = "", nchar(tidy_dataset$VDJseq[row]))

        indel_start <- deletion_list[start]-removed_nt
        indel_end <- deletion_list[start+1]-removed_nt

        ref <- paste0(germline[indel_start:(indel_end-1)], collapse = "")
        alt <- somatic[indel_start]
        tidy_dataset$indel[row] <- paste(tidy_dataset$indel[row], paste(indel_start, ref, alt, sep = "_"), sep = ",")

        if(indel_start <= tidy_dataset$len_no_CDR3[row] | indel_start > sum(as.numeric(strsplit(tidy_dataset$arch[row], "-")[[1]][1:6]))){
          somatic[(indel_start+1):(indel_end-1)] <- germline[(indel_start+1):(indel_end-1)]
          tidy_dataset$VDJseq[row] <- paste0(somatic, collapse = "")

          ## Update IG gene architecture removing InDels
          arch_list <- rep(c("FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4"), as.numeric(strsplit(tidy_dataset$arch[row], "-")[[1]]))
          arch_deleted <- arch_list[(indel_start+1)]
          n_del <- as.numeric(indel_end-indel_start-1)

          arch_df <- as.data.frame(t(as.data.frame(table(c(arch_list, rep(arch_deleted, n_del))))))
          colnames(arch_df) <- arch_df[1,]
          for(reg in c("FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4")){
            if(!reg %in% colnames(arch_df)){arch_df[[reg]][2] <- "0"}
          }
          arch_df$FR4[2] <- nchar(tidy_dataset$VDJseq[row]) - sum(as.numeric(arch_df[2, c("FR1","CDR1","FR2","CDR2","FR3","CDR3")])) ## Recalculate FR4 after InDel correction

          tidy_dataset$arch[row] <- paste0(as.numeric(arch_df[2, c("FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4")]), collapse = "-")
          tidy_dataset$len_no_CDR3[row] <- sum(as.numeric(arch_df[2, c("FR1","CDR1","FR2","CDR2","FR3")]))

        } else {
          somatic <- somatic[-((indel_start+1):(indel_end-1))]
          tidy_dataset$VDJseq[row] <- paste0(somatic, collapse = "")

          germline <- germline[-((indel_start+1):(indel_end-1))]
          tidy_dataset$germline_alignment[row] <- paste0(germline, collapse = "")

          removed_nt <- removed_nt + length((indel_start+1):(indel_end-1))
        }
      }
    }

    ## Correct Insertions
    if(length(insertion_list) > 0){
      removed_nt <- 0
      for(start in seq(1, length(insertion_list), 2)){
        somatic <- str_split_fixed(tidy_dataset$VDJseq[row], pattern = "", nchar(tidy_dataset$VDJseq[row]))
        germline <- str_split_fixed(tidy_dataset$germline_alignment[row], pattern = "", nchar(tidy_dataset$germline_alignment[row]))

        indel_start <- insertion_list[start]-removed_nt
        indel_end <- insertion_list[start+1]-removed_nt

        ref <- germline[indel_start]
        alt <- paste0(somatic[indel_start:(indel_end-1)], collapse = "")
        tidy_dataset$indel[row] <- paste(tidy_dataset$indel[row], paste(indel_start, ref, alt, sep = "_"), sep = ",")

        if(indel_start <= tidy_dataset$len_no_CDR3[row] | indel_start > sum(as.numeric(strsplit(tidy_dataset$arch[row], "-")[[1]][1:6]))){

          tidy_dataset$VDJseq[row] <- paste0(somatic[,-((indel_start+1):(indel_end-1))], collapse = "")
          tidy_dataset$germline_alignment[row] <- paste0(germline[,-((indel_start+1):(indel_end-1))], collapse = "")

          ## Update IG gene architecture removing InDels
          arch_list <- rep(c("FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4"), as.numeric(strsplit(tidy_dataset$arch[row], "-")[[1]]))
          arch_df <- as.data.frame(t(as.data.frame(table(arch_list[-((indel_start+1):(indel_end-1))]))))
          colnames(arch_df) <- arch_df[1,]
          for(reg in c("FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4")){
            if(!reg %in% colnames(arch_df)){arch_df[[reg]][2] <- "0"}
          }
          arch_df$FR4[2] <- nchar(tidy_dataset$VDJseq[row]) - sum(as.numeric(arch_df[2, c("FR1","CDR1","FR2","CDR2","FR3","CDR3")])) ## Recalculate FR4 after InDel correction

          tidy_dataset$arch[row] <- paste0(as.numeric(arch_df[2, c("FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4")]), collapse = "-")
          tidy_dataset$len_no_CDR3[row] <- sum(as.numeric(arch_df[2, c("FR1","CDR1","FR2","CDR2","FR3")]))

          removed_nt <- removed_nt + length((indel_start+1):(indel_end-1))

        } else{
          germline[(indel_start+1):(indel_end-1)] <- somatic[(indel_start+1):(indel_end-1)]
          tidy_dataset$germline_alignment[row] <- paste0(germline, collapse = "")
        }
      }
    }
  }
  tidy_dataset$indel <- gsub("NA,","",tidy_dataset$indel)

  return(tidy_dataset)
}
.translate_sequence <- function(nucleotide_sequence) {
  # Define the genetic code mapping codons to amino acids
  genetic_code <- list(
    "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
    "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "TAT" = "Y", "TAC" = "Y", "TAA" = "*", "TAG" = "*",
    "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
    "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
    "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
    "TGT" = "C", "TGC" = "C", "TGA" = "*", "TGG" = "W",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
  )

  # Initialize an empty amino acid sequence
  amino_acid_sequence <- character(0)
  # Convert the nucleotide sequence to uppercase
  nucleotide_sequence <- toupper(nucleotide_sequence)

  # Split the sequence into codons and translate each codon
  for (i in seq(1, nchar(nucleotide_sequence), by = 3)) {
    codon <- substr(nucleotide_sequence, i, i + 2)
    if(nchar(codon) < 3){break}
    # Check if the codon is in the genetic code
    if (codon %in% names(genetic_code)) {
      amino_acid <- genetic_code[[codon]]
      amino_acid_sequence <- c(amino_acid_sequence, amino_acid)
    } else {
      amino_acid <- "X"
      amino_acid_sequence <- c(amino_acid_sequence, amino_acid)
    }
  }
  # Join the amino acids to form the final protein sequence
  protein_sequence <- paste(amino_acid_sequence, collapse = "")
  return(protein_sequence)
}
.update_junction <- function(row){
  start_CDR3 <- sum(as.numeric(strsplit(row[14], split = "-")[[1]][1:5]))+1
  end_CDR3 <- sum(as.numeric(strsplit(row[14], split = "-")[[1]][1:6]))
  new_cdr3_nt <- paste0(strsplit(row[6], split = "")[[1]][(start_CDR3-3):(end_CDR3+3)], collapse = "")
  return(new_cdr3_nt)
}
.assign_Clonotypes <- function(tidy_dataset, hc_similarity_cutoff, hc_mode, cdr3_mode, cdr3_InDel_correction_mode, summary_file, analysis_mode, name, data_type, total_tasks, completed_tasks){
  tidy_dataset$Vgene <- sapply(tidy_dataset$clonotypeLabel, function(x) strsplit(x, split = "_")[[1]][1])
  tidy_dataset$CDR3nt <- sapply(tidy_dataset$junction, function(x) paste0(strsplit(x, split = "")[[1]][4:(nchar(x)-3)], collapse = ""))
  tidy_dataset$CDR3aa <- sapply(tidy_dataset$junction_aa, function(x) paste0(strsplit(x, split = "")[[1]][2:(nchar(x)-1)], collapse = ""))

  if(cdr3_mode == "nt"){
    tidy_dataset$CDR3len <- nchar(tidy_dataset$CDR3nt)
  } else if(cdr3_mode == "aa"){
    tidy_dataset$CDR3len <- nchar(tidy_dataset$CDR3aa)
  }

  if(data_type == "missionbio"){
    tidy_dataset$Pre_Clonotype <- paste0(tidy_dataset$Vgene,"_",tidy_dataset$CDR3len)
  } else{
    tidy_dataset$Pre_Clonotype <- paste0(tidy_dataset$Vgene,"_",tidy_dataset$CDR3len, "_", tidy_dataset$len_no_CDR3)
  }

  write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Starting IG Clonotype hierarchical clustering... ",completed_tasks,"/",total_tasks," tasks completed."), file = summary_file, append = T)
  message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Starting IG Clonotype hierarchical clustering.... ",completed_tasks,"/",total_tasks," tasks completed."))
  completed_tasks <- completed_tasks + 1

  tmpSubset <- tidy_dataset[, c("Pre_Clonotype", ifelse(cdr3_mode == "nt", "CDR3nt", "CDR3aa"))]
  tmpSubset <- tmpSubset[!duplicated(tmpSubset),]
  colnames(tmpSubset)[2] <- "CDR3"

  if(hc_similarity_cutoff == "automatic"){
    recommended_cutoff <- ifelse(cdr3_mode == "nt", 0.2, 0.3) ## 0.2 for nucleotides and 0.3 for amino acids, to be confirmed!!

    calc_AutoThreshold <- .calculateClonotypeThreshold(tmpSubset)
    automatic_similarity_cutoff <- calc_AutoThreshold[[1]]

    ## Avoid small automatic cutoffs in extremely clonal samples!
    if(automatic_similarity_cutoff <= recommended_cutoff){
      hc_similarity_cutoff <- recommended_cutoff
    } else{
      hc_similarity_cutoff <- automatic_similarity_cutoff
    }
    .plot_automatic_threshold_density(calc_AutoThreshold, outputDir, sample, cdr3_mode, recommended_cutoff, hc_similarity_cutoff) ## IAN ARREGLAR igblast_dir!

  } else{
    hc_similarity_cutoff <- as.numeric(hc_similarity_cutoff)
  }

  tmpSubset <- tmpSubset %>%
    group_by(Pre_Clonotype) %>%
    mutate(Clonotype = as.character(paste0(Pre_Clonotype, "_C", .compute_clusters(CDR3, hc_mode, hc_similarity_cutoff))))

  if(cdr3_mode=="nt"){
    tidy_dataset$Clonotype <- tmpSubset$Clonotype[match(paste0(tidy_dataset$Pre_Clonotype, "_", tidy_dataset$CDR3nt), paste0(tmpSubset$Pre_Clonotype, "_", tmpSubset$CDR3))]
  }else{
    tidy_dataset$Clonotype <- tmpSubset$Clonotype[match(paste0(tidy_dataset$Pre_Clonotype, "_", tidy_dataset$CDR3aa), paste0(tmpSubset$Pre_Clonotype, "_", tmpSubset$CDR3))]
  }

  ## Merge clonotypes based on CDR3 indels
  write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Starting CDR3 InDel correction... ",completed_tasks,"/",total_tasks," tasks completed."), file = summary_file, append = T)
  message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Starting CDR3 InDel correction... ",completed_tasks,"/",total_tasks," tasks completed."))
  completed_tasks <- completed_tasks + 1

  if(cdr3_InDel_correction_mode == "no"){
    tidy_dataset$CorrectClt <- tidy_dataset$Clonotype
  } else{
    tidy_dataset$CorrectClt <- NA
    tidy_dataset$CorrectClt <- as.character(.correct_Clonotype_InDels(tidy_dataset, cdr3_mode, hc_similarity_cutoff, cdr3_InDel_correction_mode, data_type))
  }

  ## Undo the re-naming of V genes once the clonotype clustering has been performed
  write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Undoing IGV gene re-naming performed before clonotype clustering... ",completed_tasks,"/",total_tasks," tasks completed."), file = summary_file, append = T)
  message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - ", ifelse(analysis_mode == "joint", "Case_", "Sample_"), name, ": Undoing IGV gene re-naming performed before clonotype clustering... ",completed_tasks,"/",total_tasks," tasks completed."))
  completed_tasks <- completed_tasks + 1

  comb_name <- unique(tidy_dataset$CorrectClt)[which(grepl("\\.", unique(tidy_dataset$CorrectClt)))]
  comb_list <- unname(unlist(sapply(comb_name, function(x) strsplit(x, "\\."))))
  rep_elem <- comb_list[duplicated(comb_list)]
  if(length(rep_elem) > 0){
    for(i in 1:length(rep_elem)){
      tidy_dataset$CorrectClt[tidy_dataset$CorrectClt %in% comb_name[grep(rep_elem[i], comb_name)]] <- paste(sort(unique(unname(unlist(sapply(comb_name[grep(rep_elem[i], comb_name)], function(x) strsplit(x, "\\."))))), decreasing = F), collapse = ".")
    }
  }

  ## Make clonotypeIDs
  if(data_type %in% c("single_cell", "missionbio")){
    dict_cltp <- aggregate(x = tidy_dataset$CorrectClt, by = list(tidy_dataset$CorrectClt), FUN = length)
    dict_cltp <- dict_cltp[order(dict_cltp$x, decreasing = T), ]
    dict_cltp$IGchain <- sapply(dict_cltp$Group.1, FUN = function(x) paste0(unlist(strsplit(x, split = ""))[1:3], collapse = ""))
    for(chain in sort(unique(dict_cltp$IGchain))){
      dict_cltp$clonotypeID[dict_cltp$IGchain == chain] <- paste0(chain,".C",1:nrow(dict_cltp[dict_cltp$IGchain == chain,]),sep = "")
    }
    tidy_dataset <- merge(tidy_dataset, dict_cltp[,c(1,4)], by.x = "CorrectClt",by.y = "Group.1", all.x = TRUE)

    tidy_dataset <- .make_CV_and_subclones_single_cell(tidy_dataset)

  } else if(data_type == "bulk"){

    tidy_dataset <- tidy_dataset %>%
      group_by(Clonotype) %>%
      mutate(Clonotype_nreads = as.numeric(sum(total_reads_unique_seq)))

    tidy_dataset$Clonotype_freq <- as.numeric((tidy_dataset$Clonotype_nreads / sum(tidy_dataset$total_reads_unique_seq))*100)

    tidy_dataset$CorrectedClt_nreads=tidy_dataset$CorrectedClt_freq=NA
    for(corr_cltp in unique(tidy_dataset$CorrectClt)){
      tmp_df <- as.data.frame(tidy_dataset[tidy_dataset$CorrectClt == corr_cltp,])
      tmp_df <- tmp_df[!duplicated(tmp_df$Clonotype),]

      tidy_dataset$CorrectedClt_nreads[tidy_dataset$CorrectClt == corr_cltp] <- sum(tmp_df$Clonotype_nreads)
      tidy_dataset$CorrectedClt_freq[tidy_dataset$CorrectClt == corr_cltp] <- sum(tmp_df$Clonotype_freq)
    }
    tidy_dataset <- .make_sequence_IDs_bulk(tidy_dataset)
  }

  return_fields <- c("sample", "sequence_id", "sequence", "raw_VDJseq", "VDJseq", "c_call", "germline_alignment",
                     "productive", "VDJ", "len_no_CDR3", "CDR3aa", "CDR3len", "junction_aa",
                     "junction_aa_length", "arch", "indel", "VDJseq_indels", "Unique_SequenceID",
                     "Clonotype", "CorrectClt", "clonotypeID", "clonotypeVariantID_in_Cltp")

  if(data_type == "bulk"){
    return_fields <- c(return_fields, c("CorrectedClt_nreads", "CorrectedClt_freq", "subcloneID_in_CV", "subcloneID_in_Clt", "total_reads_unique_seq"))
  }
  return(list(as.data.frame(tidy_dataset[,return_fields]), completed_tasks))
}
.annotate_Immunogenetics <- function(tidy_dataset, annotate_CLL_immGen, annotate_satellite_subsets, annotate_ags, data_type){

  if(data_type == "missionbio"){
    tidy_dataset$CorrectClonotypes_Consensus_Germline <- tidy_dataset$germline_alignment
    tidy_dataset$VDJseq_correctedCDR3 <- tidy_dataset$VDJseq

  } else{
    tidy_dataset$CorrectClonotypes_Consensus_Germline <- NA
    correct_germline_subMat <- matrix(-1, nrow = 6, ncol = 6, dimnames = list(c("A", "T", "C", "G", "X", "N"), c("A", "T", "C", "G", "X", "N")))
    diag(correct_germline_subMat) <- 1
    correct_germline_subMat[nrow(correct_germline_subMat),] <- 1
    correct_germline_subMat[,ncol(correct_germline_subMat)] <- 1

    tidy_dataset$VDJseq_correctedCDR3 <- tidy_dataset$VDJseq
    for(cltp in unique(tidy_dataset$clonotypeID)){
      cons_germ <- .consensus_germline(tidy_dataset[tidy_dataset$clonotypeID == cltp,])
      tidy_dataset$CorrectClonotypes_Consensus_Germline[tidy_dataset$clonotypeID == cltp] <- cons_germ[[2]]

      seqs_to_correct <- unique(tidy_dataset$Unique_SequenceID[tidy_dataset$clonotypeID == cltp & tidy_dataset$CDR3len != cons_germ[[1]]])
      if(length(seqs_to_correct) > 0){
        for(seq in seqs_to_correct){
          corr_germ <- .correct_consensus_germline(tidy_dataset$VDJseq[tidy_dataset$Unique_SequenceID == seq][1], tidy_dataset$germline_alignment[tidy_dataset$Unique_SequenceID == seq][1], cons_germ[[2]], correct_germline_subMat)
          tidy_dataset$indel[tidy_dataset$Unique_SequenceID == seq] <- paste0(unique(c(strsplit(tidy_dataset$indel[tidy_dataset$Unique_SequenceID == seq][1], ",")[[1]], strsplit(corr_germ[[1]], ",")[[1]])), collapse = ",")
          tidy_dataset$VDJseq_correctedCDR3[tidy_dataset$Unique_SequenceID == seq] <- corr_germ[[2]]
        }
      }
    }
    tidy_dataset$indel <- gsub("NA,", "", tidy_dataset$indel)
    tidy_dataset$indel <- gsub("^,", "", tidy_dataset$indel)
  }

  ## Compute Videntity by Corrected Clonotypes
  for(row in 1:nrow(tidy_dataset)){
    tidy_dataset$v_identity[row] <- .calculate_v_ident(tidy_dataset$VDJseq[row], tidy_dataset$CorrectClonotypes_Consensus_Germline[row], as.numeric(tidy_dataset$len_no_CDR3[row]))
  }

  if(annotate_CLL_immGen){
    ## Assign CLL subsets
    tidy_dataset$IGHsubset <- apply(tidy_dataset, 1, function(x) .assign_CLL_subsets(x[which(colnames(tidy_dataset) == "productive")],
                                                                                    x[which(colnames(tidy_dataset) == "VDJ")],
                                                                                    x[which(colnames(tidy_dataset) == "junction_aa")],
                                                                                    annotate_satellite_subsets))
    ## Annotating CLL R110 mutation
    tidy_dataset$VDJ <- paste0(tidy_dataset$VDJ, apply(tidy_dataset, 1, function(x) .annotateR110mut(x[which(colnames(tidy_dataset) == "VDJ")],
                                                                                                    x[which(colnames(tidy_dataset) == "VDJseq")])))
  }

  ## Translate nucleotide sequence into amino-acid
  for(seq in unique(tidy_dataset$VDJseq)){
    tidy_dataset$VDJseq_aa[tidy_dataset$VDJseq == seq] <- .translate_sequence(seq)
  }

  for(seq in unique(tidy_dataset$VDJseq_correctedCDR3)){
    tidy_dataset$VDJseq_correctedCDR3_aa[tidy_dataset$VDJseq_correctedCDR3 == seq] <- .translate_sequence(seq)
  }

  for(seq in unique(tidy_dataset$CorrectClonotypes_Consensus_Germline)){
    tidy_dataset$CorrectClonotypes_Consensus_Germline_aa[tidy_dataset$CorrectClonotypes_Consensus_Germline == seq] <- .translate_sequence(seq)
  }

  if(annotate_ags){
    tidy_dataset$VDJ <- paste0(tidy_dataset$VDJ, apply(tidy_dataset, 1, function(x) .annotateAGS(x[which(colnames(tidy_dataset) == "productive")],
                                                                                                 x[which(colnames(tidy_dataset) == "VDJ")],
                                                                                                 x[which(colnames(tidy_dataset) == "arch")],
                                                                                                 x[which(colnames(tidy_dataset) == "VDJseq_aa")],
                                                                                                 x[which(colnames(tidy_dataset) == "germline_alignment")])))
  }

  ## Calculate the consensus CDR3 aa sequence without unproductive sublcones
  cons_CDR3_df <- tidy_dataset %>%
    group_by(clonotypeID) %>%
    filter(productive == "productive") %>%
    mutate(Clonotype_ConsensusCDR3 = as.character(.consensus_CDR3(CDR3aa)))
  tidy_dataset$Clonotype_ConsensusCDR3 <- cons_CDR3_df$Clonotype_ConsensusCDR3[match(tidy_dataset$VDJseq_indels, cons_CDR3_df$VDJseq_indels)]

  if(annotate_CLL_immGen){
    tidy_dataset$Clonotype_Subset <- apply(tidy_dataset, MARGIN = 1, function(x) .assign_CLL_subsets_to_Clonotypes(x[which(colnames(tidy_dataset) == "productive")],
                                                                                                                   x[which(colnames(tidy_dataset) == "VDJ")],
                                                                                                                   x[which(colnames(tidy_dataset) == "Clonotype_ConsensusCDR3")],
                                                                                                                   annotate_satellite_subsets))
  }

  return(tidy_dataset)
}

## Function for computing the hierarchical clustering of clonotypes
.compute_identity <- function(seq1, seq2){
  ident <- (sum(strsplit(seq1,"")[[1]] == strsplit(seq2,"")[[1]])/nchar(seq1))
  return(ident)
}
.compute_clusters <- function(CDR3seqs, hc_mode, hc_similarity_cutoff){
  if(length(CDR3seqs) == 1){return(1)
  } else{

    distance_matrix <- matrix(0, nrow = length(CDR3seqs), ncol = length(CDR3seqs), dimnames = list(CDR3seqs,CDR3seqs))
    indices <- which(upper.tri(distance_matrix, diag = FALSE), arr.ind = TRUE) ## Select only upper part of the matrix

    for (i in 1:nrow(indices)) {
      x <- indices[i, 1]
      y <- indices[i, 2]
      distance_matrix[x, y] <- .compute_identity(CDR3seqs[x], CDR3seqs[y])
      distance_matrix[y, x] <- distance_matrix[x, y]
    }
    diag(distance_matrix) <- 1

    hc <- hclust(as.dist(1 - distance_matrix), method = hc_mode)
    hc$height <- round(hc$height, 6)
    clusters <- cutree(hc, h = hc_similarity_cutoff)
    clusters_dict <- data.frame(CDR3 = names(clusters), Cluster = clusters, row.names = NULL)
    clusters_v <- clusters_dict$Cluster[match(CDR3seqs, clusters_dict$CDR3)]
    return(clusters_v)
  }
}

## Functions to correct CDR3 artifacts/InDels
.compute_clonotype_similarity <- function(distance_matrix, x, y, cdr3_mode, substitution_matrix){

  seq1 <- rownames(distance_matrix)[x]
  seq2 <- colnames(distance_matrix)[y]

  if(cdr3_mode == "nt"){
    aln <- pairwiseAlignment(pattern = c(seq1,seq2)[which.max(nchar(c(seq1,seq2)))],
                             subject = c(seq1,seq2)[which.min(nchar(c(seq1,seq2)))],
                             type = "global")

  } else if(cdr3_mode == "aa"){ ## Added substitution matrix for aa alignment to avoid penalty for physico-chemical property mismatch
    aln <- pairwiseAlignment(pattern = c(seq1,seq2)[which.max(nchar(c(seq1,seq2)))],
                             subject = c(seq1,seq2)[which.min(nchar(c(seq1,seq2)))],
                             type = "global", substitutionMatrix = substitution_matrix)
  }

  pattern <- strsplit(gsub("pattern: ", "", capture.output(aln)[2]), "")[[1]]
  subject = correct_subject = strsplit(gsub("subject: ", "", capture.output(aln)[3]), "")[[1]]

  if(any(grepl("-", pattern)) | aln@score <= 0){
    return(1)
  } else{
    correct_subject[which(correct_subject == "-")] <- pattern[which(correct_subject == "-")]
    similarity <- 1-((sum(pattern == correct_subject) - (sum(!diff(which(subject == "-")) == 1)+1)) / nchar(paste(pattern, collapse = "")))
    return(similarity)
  }
}
.compare_similar_clonotypes <- function(data_frame_by_sample, c1, c2, cdr3_mode, cdr3_InDel_correction_mode, hc_similarity_cutoff, substitution_matrix){
  if(cdr3_InDel_correction_mode == "hard_filter"){
    if(cdr3_mode == "nt"){
      distance_matrix <- matrix(0, nrow = length(unique(data_frame_by_sample$CDR3nt[data_frame_by_sample$Clonotype == c1])),
                                ncol = length(unique(data_frame_by_sample$CDR3nt[data_frame_by_sample$Clonotype == c2])),
                                dimnames = list(unique(data_frame_by_sample$CDR3nt[data_frame_by_sample$Clonotype == c1]),
                                                unique(data_frame_by_sample$CDR3nt[data_frame_by_sample$Clonotype == c2])))

    }else{
      distance_matrix <- matrix(0, nrow = length(unique(data_frame_by_sample$CDR3aa[data_frame_by_sample$Clonotype == c1])),
                                ncol = length(unique(data_frame_by_sample$CDR3aa[data_frame_by_sample$Clonotype == c2])),
                                dimnames = list(unique(data_frame_by_sample$CDR3aa[data_frame_by_sample$Clonotype == c1]),
                                                unique(data_frame_by_sample$CDR3aa[data_frame_by_sample$Clonotype == c2])))
    }
  } else if(cdr3_InDel_correction_mode == "soft_filter"){
    if(cdr3_mode == "nt"){
      distance_matrix <- matrix(0, nrow = length(.possibleCDR3(.consensus_CDR3(unique(data_frame_by_sample$CDR3nt[data_frame_by_sample$Clonotype == c1])))),
                                ncol = length(.possibleCDR3(.consensus_CDR3(unique(data_frame_by_sample$CDR3nt[data_frame_by_sample$Clonotype == c2])))),
                                dimnames = list(.possibleCDR3(.consensus_CDR3(unique(data_frame_by_sample$CDR3nt[data_frame_by_sample$Clonotype == c1]))),
                                                .possibleCDR3(.consensus_CDR3(unique(data_frame_by_sample$CDR3nt[data_frame_by_sample$Clonotype == c2])))))

    }else{
      distance_matrix <- matrix(0, nrow = length(.possibleCDR3(.consensus_CDR3(unique(data_frame_by_sample$CDR3aa[data_frame_by_sample$Clonotype == c1])))),
                                ncol = length(.possibleCDR3(.consensus_CDR3(unique(data_frame_by_sample$CDR3aa[data_frame_by_sample$Clonotype == c2])))),
                                dimnames = list(.possibleCDR3(.consensus_CDR3(unique(data_frame_by_sample$CDR3aa[data_frame_by_sample$Clonotype == c1]))),
                                                .possibleCDR3(.consensus_CDR3(unique(data_frame_by_sample$CDR3aa[data_frame_by_sample$Clonotype == c2])))))
    }
  }

  for(x in 1:nrow(distance_matrix)){
    for(y in 1:ncol(distance_matrix)){
      fill_value <- .compute_clonotype_similarity(distance_matrix, x, y, cdr3_mode, substitution_matrix)
      distance_matrix[x,y] <- fill_value
      if(fill_value > hc_similarity_cutoff){break}
    }
    if(fill_value > hc_similarity_cutoff){break}
  }

  return(max(distance_matrix))
}
.correct_Clonotype_InDels <- function(data_frame_by_sample, cdr3_mode, hc_similarity_cutoff, cdr3_InDel_correction_mode, data_type){

  if(cdr3_mode == "aa"){
    amino_acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X", "*")
    substitution_matrix <- matrix(-1, nrow = 22, ncol = 22, dimnames = list(amino_acids, amino_acids))
  } else if(cdr3_mode == "nt"){
    nucleotides <- c("A", "T", "C", "G", "X", "N")
    substitution_matrix <- matrix(-1, nrow = 6, ncol = 6, dimnames = list(nucleotides, nucleotides))
  }
  diag(substitution_matrix) <- 1

  for(clone in unique(data_frame_by_sample$Clonotype)){
    list_Clt <- c(clone)
    if(all(is.na(data_frame_by_sample$CorrectClt[data_frame_by_sample$Clonotype == clone]))){
      vgene <- unique(data_frame_by_sample$Vgene[data_frame_by_sample$Clonotype == clone])
      len_cdr3 <- unique(data_frame_by_sample$CDR3len[data_frame_by_sample$Clonotype == clone])
      if(data_type == "missionbio"){
        for(candidate_indel in unique(data_frame_by_sample$Clonotype[data_frame_by_sample$Vgene == vgene & data_frame_by_sample$CDR3len != len_cdr3])){
          if(.compare_similar_clonotypes(data_frame_by_sample, clone, candidate_indel, cdr3_mode, cdr3_InDel_correction_mode, hc_similarity_cutoff, substitution_matrix) <= hc_similarity_cutoff){
            list_Clt <- c(list_Clt, candidate_indel)
          }
        }
      } else{
        len_no_cdr3 <- unique(data_frame_by_sample$len_no_CDR3[data_frame_by_sample$Clonotype == clone])
        for(candidate_indel in unique(data_frame_by_sample$Clonotype[data_frame_by_sample$Vgene == vgene & data_frame_by_sample$CDR3len != len_cdr3 & data_frame_by_sample$len_no_CDR3 == len_no_cdr3])){
          if(.compare_similar_clonotypes(data_frame_by_sample, clone, candidate_indel, cdr3_mode, cdr3_InDel_correction_mode, hc_similarity_cutoff, substitution_matrix) <= hc_similarity_cutoff){
            list_Clt <- c(list_Clt, candidate_indel)
          }
        }
      }
      data_frame_by_sample$CorrectClt[data_frame_by_sample$Clonotype %in% list_Clt] <- paste(sort(list_Clt), collapse = ".")
    }
  }
  add_to_df <- data_frame_by_sample$CorrectClt
  return(add_to_df)
}

## Functions for the automatic calculation of clonotype similarity cutoff
.neighborDist <- function(seq1, seq2){
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]
  dist <- sum(s1 != s2)/length(s1)
  return(dist)
}
.calculateClonotypeThreshold <- function(tmpSubset){
  neighDist_v <- c()
  for(pC in unique(tmpSubset$Pre_Clonotype)){
    chain_CDR3 <- tmpSubset$CDR3[tmpSubset$Pre_Clonotype == pC]
    if(length(chain_CDR3) == 1){next}
    chain_len_dist <- sapply(chain_CDR3, function(seq1) min(sapply(chain_CDR3[chain_CDR3 != seq1], function(seq2) .neighborDist(seq1, seq2))))
    neighDist_v <- c(neighDist_v, as.numeric(chain_len_dist))
  }
  ## Extract the value with the valley between the two distributions
  density <- density(neighDist_v)
  peaks <- findpeaks(density$y)
  valleys <- findpeaks(-density$y)

  ## Percentage of identity of the two modes
  modes <- density$x[peaks[, 2]]

  bimodal_cutoff <- density$x[valleys[, 2][which.min(valleys[, 1])]]
  return(list(bimodal_cutoff, neighDist_v))
}
.plot_automatic_threshold_density <- function(calc_AutoThreshold, outputDir, name, cdr3_mode, recommended_cutoff, hc_similarity_cutoff){

  OutPlotPath <- paste0(outputDir,"annotation_results/", ifelse(analysis_mode == "joint", "case_", "sample_"), name, "_Automatic_ClonotypeSimilarity_cutoff.pdf")

  auto_threshold <- calc_AutoThreshold[[1]]
  identity_distr <- calc_AutoThreshold[[2]]

  data <- data.frame(identity1 = identity_distr, identity2 = identity_distr)

  min_xlim <- min(identity_distr)-0.01
  xlab <- paste0("CDR3 ", cdr3_mode, " distance")

  title_plot <- ifelse(auto_threshold <= recommended_cutoff,
                       paste0(xlab, " - Forced threshold to ", recommended_cutoff, " due to high clonality"),
                       paste0(xlab, " - Threshold found: ", hc_similarity_cutoff))

  max_xlim <- ifelse(max(identity_distr) <= recommended_cutoff, recommended_cutoff, max(identity_distr)+0.1)

  p1 <- ggplot(data, aes(x = identity1)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    labs(title = title_plot, x = paste0("\n",xlab), y = "Density\n") +
    theme_minimal() + xlim(c(min_xlim, max_xlim)) +
    geom_vline(xintercept = hc_similarity_cutoff,
               linetype = "dotted", col = "darkred", linewidth = 1.5)

  pdf(OutPlotPath, width = 10, height = 5, useDingbats = F)
  print(p1)
  dev.off()
}

## Function to create clonotypeVariant and Subclone IDs
.make_CV_and_subclones_single_cell <- function(tidy_dataset){
  tidy_dataset$clonotypeVariantID_in_Cltp <- NA
  ## Make ClonotypeVariant and Subclone IDs based on the Clonotype architecture
  for(cltp in unique(tidy_dataset$clonotypeID)){
    clonotypeVariant_dict <- aggregate(x = tidy_dataset$CDR3aa[tidy_dataset$clonotypeID == cltp], by = list(tidy_dataset$CDR3aa[tidy_dataset$clonotypeID == cltp]), FUN = length)
    clonotypeVariant_dict <- clonotypeVariant_dict[order(clonotypeVariant_dict$x, decreasing = T),]
    clonotypeVariant_dict$clonotypeVariantID_in_Cltp <- paste0(cltp, ".CV", 1:nrow(clonotypeVariant_dict))
    tidy_dataset$clonotypeVariantID_in_Cltp[tidy_dataset$clonotypeID == cltp] <- clonotypeVariant_dict$clonotypeVariantID_in_Cltp[match(tidy_dataset$CDR3aa[tidy_dataset$clonotypeID == cltp], clonotypeVariant_dict$Group.1)]
  }

  tidy_dataset$VDJseq_indels <- paste0(tidy_dataset$VDJseq, "-", tidy_dataset$indel)
  tidy_dataset$Unique_SequenceID <- NA
  for(cv in unique(tidy_dataset$clonotypeVariantID_in_Cltp)){
    dict_unique_seq <- aggregate(x = tidy_dataset$VDJseq_indels[tidy_dataset$clonotypeVariantID_in_Cltp == cv], by = list(tidy_dataset$VDJseq_indels[tidy_dataset$clonotypeVariantID_in_Cltp == cv]), FUN = length)
    dict_unique_seq <- dict_unique_seq[order(dict_unique_seq$x, decreasing = T),]
    dict_unique_seq$Unique_SequenceID <- paste0(cv, ".S", 1:nrow(dict_unique_seq))
    tidy_dataset$Unique_SequenceID[tidy_dataset$clonotypeVariantID_in_Cltp == cv] <- dict_unique_seq$Unique_SequenceID[match(tidy_dataset$VDJseq_indels[tidy_dataset$clonotypeVariantID_in_Cltp == cv], dict_unique_seq$Group.1)]
  }

  return(tidy_dataset)
}
.make_sequence_IDs_bulk <- function(tidy_dataset){
  tidy_dataset <- tidy_dataset[order(tidy_dataset$CorrectedClt_freq, decreasing = T),]

  dict_Clt <- tidy_dataset[!duplicated(tidy_dataset$CorrectClt),c("CorrectClt","CorrectedClt_freq")]
  dict_Clt$chain <- substr(dict_Clt$CorrectClt, 1, 3)
  dict_Clt$clonotypeID <- paste0("C",1:nrow(dict_Clt))
  dict_Clt$clonotypeID <- paste0(dict_Clt$chain, ".", dict_Clt$clonotypeID)
  tidy_dataset$clonotypeID <- dict_Clt$clonotypeID[match(tidy_dataset$CorrectClt, dict_Clt$CorrectClt)]

  clonotype_variant_dict <- tidy_dataset[order(tidy_dataset$total_reads_unique_seq, decreasing = T), c("Vgene", "CDR3aa","total_reads_unique_seq","clonotypeID")]
  clonotype_variant_dict$clonotypeVariantLabel <- paste0(clonotype_variant_dict$Vgene,"_",clonotype_variant_dict$CDR3aa)
  cv_nreads <- aggregate(x = clonotype_variant_dict$total_reads_unique_seq, by = list(clonotype_variant_dict$clonotypeVariantLabel), FUN = sum)
  clonotype_variant_dict <- merge(clonotype_variant_dict, cv_nreads, by.x = "clonotypeVariantLabel", by.y = "Group.1")
  clonotype_variant_dict <- clonotype_variant_dict[!duplicated(clonotype_variant_dict$clonotypeVariantLabel),]
  clonotype_variant_dict <- clonotype_variant_dict[order(clonotype_variant_dict$x, decreasing = T),]

  clonotype_variant_dict$clonotypeVariantID_in_Cltp <- NA
  for(clone in unique(clonotype_variant_dict$clonotypeID)){
    clonotype_variant_dict$clonotypeVariantID_in_Cltp[clonotype_variant_dict$clonotypeID == clone] <- paste0(clone, ".CV", 1:nrow(clonotype_variant_dict[clonotype_variant_dict$clonotypeID == clone,]))
  }
  tidy_dataset$clonotypeVariantLabel <- paste0(tidy_dataset$Vgene,"_",tidy_dataset$CDR3aa)
  tidy_dataset$clonotypeVariantID_in_Cltp <- clonotype_variant_dict$clonotypeVariantID_in_Cltp[match(tidy_dataset$clonotypeVariantLabel, clonotype_variant_dict$clonotypeVariantLabel)]

  subclone_dict <- tidy_dataset[order(tidy_dataset$total_reads_unique_seq, decreasing = T), c("VDJseq_indels", "total_reads_unique_seq","clonotypeID", "clonotypeVariantID_in_Cltp")]
  subclone_dict$subcloneID_in_CV <- NA
  subclone_dict$subcloneID_in_Clt <- NA
  for(cv in unique(subclone_dict$clonotypeVariantID_in_Cltp)){
    subclone_dict$subcloneID_in_CV[subclone_dict$clonotypeVariantID_in_Cltp == cv] <- paste0(cv, ".S", 1:nrow(subclone_dict[subclone_dict$clonotypeVariantID_in_Cltp == cv,]))
  }
  for(clt in unique(subclone_dict$clonotypeID)){
    subclone_dict$subcloneID_in_Clt[subclone_dict$clonotypeID == clt] <- paste0(clt, ".", 1:nrow(subclone_dict[subclone_dict$clonotypeID == clt,]))
  }
  tidy_dataset$subcloneID_in_CV <- subclone_dict$subcloneID_in_CV[match(tidy_dataset$VDJseq_indels, subclone_dict$VDJseq_indels)]
  tidy_dataset$subcloneID_in_Clt <- subclone_dict$subcloneID_in_Clt[match(tidy_dataset$VDJseq_indels, subclone_dict$VDJseq_indels)]

  return(tidy_dataset)
}

## Function to combine IGHV genes
.combine_IGHV_genes <- function(ighv) {
  combined_groups <- list()
  group_index <- integer(length(ighv))

  for (i in seq_along(ighv)) {
    added <- FALSE
    for (j in seq_along(combined_groups)) {
      if (length(intersect(ighv[[i]], combined_groups[[j]])) > 0) {
        combined_groups[[j]] <- union(combined_groups[[j]], ighv[[i]])
        group_index[i] <- j
        added <- TRUE
        break
      }
    }
    if (!added) {
      combined_groups[[length(combined_groups) + 1]] <- ighv[[i]]
      group_index[i] <- length(combined_groups)
    }
  }

  return(combined_groups)
}

## Function to combine IGH and IGK/IGL clonotypeIDs
.combine_clonotypeID_by_chains  <- function(tidy_dataset){
  tmp <- aggregate(x = tidy_dataset$clonotypeID, by = list(tidy_dataset$barcode), FUN = function(x) paste(sort(x), collapse = "-"))
  colnames(tmp)[2] <- "merge_clonotypeID"
  tidy_dataset$merge_clonotypeID <- tmp$merge_clonotypeID[match(tidy_dataset$barcode, tmp$Group.1)]

  tmp <- aggregate(x = tidy_dataset$clonotypeVariantID_in_Cltp, by = list(tidy_dataset$barcode), FUN = function(x) paste(sort(x), collapse = "-"))
  colnames(tmp)[2] <- "merge_clonotypeVariantID"
  tidy_dataset$merge_clonotypeVariantID <- tmp$merge_clonotypeVariantID[match(tidy_dataset$barcode, tmp$Group.1)]

  tmp <- aggregate(x = tidy_dataset$Unique_SequenceID, by = list(tidy_dataset$barcode), FUN = function(x) paste(sort(x), collapse = "-"))
  colnames(tmp)[2] <- "merge_subcloneID"
  tidy_dataset$merge_subcloneID <- tmp$merge_subcloneID[match(tidy_dataset$barcode, tmp$Group.1)]

  return(tidy_dataset)
}

## Function to combine IGH and IGK/IGL IG metadata
.combine_IG_metadata_by_chain <- function(tidy_dataset, annotate_CLL_immGen){

  n_cols <- ifelse(annotate_CLL_immGen, 9, 8)
  tidy_dataset[paste0("V", 1:n_cols)] <- NA

  for(subclone in unique(tidy_dataset$merge_subcloneID)){
    raw_somatic <- paste(unlist(lapply(strsplit(subclone, split = "-")[[1]], function(x){ tidy_dataset$raw_VDJseq[tidy_dataset$Unique_SequenceID == x][1] })), collapse = "_")
    somatic <- paste(unlist(lapply(strsplit(subclone, split = "-")[[1]], function(x){ tidy_dataset$VDJseq_correctedCDR3[tidy_dataset$Unique_SequenceID == x][1] })), collapse = "_")
    aa_somatic <- paste(unlist(lapply(strsplit(subclone, split = "-")[[1]], function(x){ tidy_dataset$VDJseq_correctedCDR3_aa[tidy_dataset$Unique_SequenceID == x][1] })), collapse = "_")
    germline <- paste(unlist(lapply(strsplit(subclone, split = "-")[[1]], function(x){ tidy_dataset$CorrectClonotypes_Consensus_Germline[tidy_dataset$Unique_SequenceID == x][1] })), collapse = "_")
    aa_germ <- paste(unlist(lapply(strsplit(subclone, split = "-")[[1]], function(x){ tidy_dataset$CorrectClonotypes_Consensus_Germline_aa[tidy_dataset$Unique_SequenceID == x][1] })), collapse = "_")
    architecture <- paste(unlist(lapply(strsplit(subclone, split = "-")[[1]], function(x){ tidy_dataset$arch[tidy_dataset$Unique_SequenceID == x][1] })), collapse = "_")
    indels <- paste(unlist(lapply(strsplit(subclone, split = "-")[[1]], function(x){ tidy_dataset$indel[tidy_dataset$Unique_SequenceID == x][1] })), collapse = "-")
    Clonotype_ConsensusCDR3 <- paste(unlist(lapply(strsplit(subclone, split = "-")[[1]], function(x){ tidy_dataset$Clonotype_ConsensusCDR3[tidy_dataset$Unique_SequenceID == x][1] })), collapse = "_")
    return_list <- c(raw_somatic, somatic, aa_somatic, germline, aa_germ, architecture, indels, Clonotype_ConsensusCDR3)

    if(annotate_CLL_immGen){
      subsets <- paste(unlist(lapply(strsplit(subclone, split = "-")[[1]], function(x){ tidy_dataset$IGHsubset[tidy_dataset$Unique_SequenceID == x][1] })), collapse = "_")
      return_list <- c(return_list, subsets)
    }

    number_of_rows <- sum(tidy_dataset$merge_subcloneID == subclone)
    tidy_dataset[tidy_dataset$merge_subcloneID == subclone, paste0("V", 1:n_cols)] <- matrix(return_list, nrow = number_of_rows, ncol = n_cols, byrow = TRUE)
  }
  tidy_dataset[tidy_dataset$completeBCR == "Not_supported", c("merge_clonotypeID", "merge_clonotypeVariantID", "merge_subcloneID", paste0("V", 1:n_cols))] <- NA

  return(tidy_dataset)
}

## Function to correct the "completeBCR" field
.correct_completeBCR <- function(tidy_dataset, data_type, flag.alternative.clonotypes = T, flag.clonotype.doublets = T){
  dict_ct2 <- aggregate(x = tidy_dataset$merge_clonotypeID[tidy_dataset$completeBCR %in% c("Yes", "Single_chain_2", "Potential_MB_doublet")],
                        by = list(tidy_dataset$merge_clonotypeID[tidy_dataset$completeBCR %in% c("Yes","Single_chain_2", "Potential_MB_doublet")]), FUN = length)
  dict_ct2$original_completeBCR <- tidy_dataset$completeBCR[match(dict_ct2$Group.1, tidy_dataset$merge_clonotypeID)]

  if(data_type == "missionbio"){
    dict_ct2 <- dict_ct2[order(factor(dict_ct2$original_completeBCR, levels = c("Yes", "Single_chain_2", "Potential_MB_doublet")), -dict_ct2$x),]
  } else{
    dict_ct2 <- dict_ct2[order(dict_ct2$x, decreasing = T), ]
  }

  dict_ct2$New_ClonotypeID <- NA
  dict_ct2$completeBCR <- NA

  numID <- 1
  for(i in 1:nrow(dict_ct2)){
    if(!is.na(dict_ct2$New_ClonotypeID[i])){next}
    if(numID == 1){
      dict_ct2$New_ClonotypeID[i] <- numID
      numID <- numID + 1
    }else{
      curr_clpt <- strsplit(dict_ct2$Group.1[i], split = "-")[[1]]
      for(j in 1:(i-1)){
        comp_cltp <- strsplit(dict_ct2$Group.1[j], split = "-")[[1]]
        comp_completeBCR <- dict_ct2$completeBCR[j]
        if(!is.na(comp_completeBCR)){ next }

        if(flag.clonotype.doublets == T & length(curr_clpt) >= 3 & sum(curr_clpt %in% comp_cltp) >= 1){
          dict_ct2$completeBCR[i] <- "Clonotype_doublet"
          break
        }

        else if(all(curr_clpt %in% comp_cltp)){
          dict_ct2$New_ClonotypeID[i] <- dict_ct2$New_ClonotypeID[j]
          if(flag.alternative.clonotypes == T){
            dict_ct2$completeBCR[i] <- "Alternative_clonotype"
            break
          }
        }
      }
      if(is.na(dict_ct2$New_ClonotypeID[i]) & is.na(dict_ct2$completeBCR[i])){
        dict_ct2$New_ClonotypeID[i] <-  numID
        numID <- numID + 1
      }
    }
  }

  if(flag.clonotype.doublets == T){
    if(nrow(dict_ct2) == 1){
      dict_ct2$completeBCR <- dict_ct2$original_completeBCR
    } else{
      for(i in 2:nrow(dict_ct2)){
        if(!is.na(dict_ct2$completeBCR[i])){ next }
        cltps <- dict_ct2$Group.1[i]
        v <- lapply(dict_ct2$Group.1[1:(i-1)][is.na(dict_ct2$completeBCR[1:(i-1)])], function(x) strsplit(x, split = "-")[[1]])
        pos <- c()
        ch <- c()
        for(cltp in strsplit(cltps, "-")[[1]]){
          for(j in 1:length(v)){
            if(cltp %in% v[[j]]){
              pos <- c(pos, j)
              ch <- c(ch, cltp)
            }
          }
        }
        if(length(unique(ch)) > 1 & length(unique(pos)) >= 1){
          dict_ct2$completeBCR[i] <- "Clonotype_doublet"
        }
      }
    }
  }

  dict_ct2$tmp <- paste(dict_ct2$New_ClonotypeID, ifelse(dict_ct2$completeBCR != "Clonotype_doublet" | is.na(dict_ct2$completeBCR), NA, "Clonotype_doublet"), sep = "_")
  dict_ct2 <- dict_ct2[!(dict_ct2$original_completeBCR == "Single_chain_2" & is.na(dict_ct2$completeBCR)),]

  dict_ct2 <- dict_ct2 %>%
    group_by(tmp) %>%
    mutate(New_x = sum(x))

  if(data_type == "missionbio"){
    fix_df <- dict_ct2[dict_ct2$original_completeBCR == "Potential_MB_doublet" & is.na(dict_ct2$completeBCR),]
    if(nrow(fix_df) > 0){
      for(row in 1:nrow(fix_df)){
        chains <- as.character(sapply(strsplit(fix_df$Group.1[row], split = "-")[[1]], function(x){strsplit(x, split = "\\.")[[1]]})[1,])
        if((length(strsplit(fix_df$Group.1[row], split = "-")[[1]]) == 2) & !(("IGH" %in% chains) & ("IGK" %in% chains | "IGL" %in% chains))){
          dict_ct2$completeBCR[dict_ct2$Group.1 == fix_df$Group.1[row]] <- "Single_chain_2"
        } else{
          dict_ct2$completeBCR[dict_ct2$Group.1 == fix_df$Group.1[row]] <- "Yes"
        }
      }
    }
  }

  dict_ct2 <- dict_ct2[order(dict_ct2$New_x, decreasing = T), ]
  dict_ct2$igclonotypeID <- NA
  dict_ct2$igclonotypeID[!dict_ct2$completeBCR %in% c("Clonotype_doublet", "Single_chain_2") | is.na(dict_ct2$completeBCR)] <- paste0("C",as.integer(factor(paste0(dict_ct2$tmp[!dict_ct2$completeBCR %in% c("Clonotype_doublet", "Single_chain_2") | is.na(dict_ct2$completeBCR)],
                                                                                                                                           dict_ct2$New_ClonotypeID[!dict_ct2$completeBCR %in% c("Clonotype_doublet", "Single_chain_2") | is.na(dict_ct2$completeBCR)]),
                                                                                                                                    levels = unique(paste0(dict_ct2$tmp[!dict_ct2$completeBCR %in% c("Clonotype_doublet", "Single_chain_2") | is.na(dict_ct2$completeBCR)],
                                                                                                                                                           dict_ct2$New_ClonotypeID[!dict_ct2$completeBCR %in% c("Clonotype_doublet", "Single_chain_2") | is.na(dict_ct2$completeBCR)])))))

  tidy_dataset <- merge(tidy_dataset, dict_ct2[,c("Group.1","igclonotypeID")], by.x = "merge_clonotypeID", by.y = "Group.1", all.x = TRUE)
  tidy_dataset$completeBCR[tidy_dataset$merge_clonotypeID %in% dict_ct2$Group.1[dict_ct2$completeBCR %in% "Yes"]] <- "Yes"
  tidy_dataset$completeBCR[tidy_dataset$merge_clonotypeID %in% dict_ct2$Group.1[dict_ct2$completeBCR %in% "Alternative_clonotype"]] <- "Alternative_clonotype"
  tidy_dataset$completeBCR[tidy_dataset$merge_clonotypeID %in% dict_ct2$Group.1[dict_ct2$completeBCR %in% "Clonotype_doublet"]] <- "Clonotype_doublet"
  tidy_dataset$completeBCR[tidy_dataset$merge_clonotypeID %in% dict_ct2$Group.1[dict_ct2$completeBCR %in% "Single_chain_2"]] <- "Single_chain_2"

  return(tidy_dataset)
}

## Function to correct clonotype labels based on cell-level information
.correct_Clonotype_labels <- function(tidy_dataset){
  tidy_dataset$clonotypeID_new <- NA
  # split chain if present in different IG clonotypes
  for(cltp in unique(tidy_dataset$clonotypeID)){
    if(length(unique(tidy_dataset$igclonotypeID[tidy_dataset$clonotypeID == cltp & tidy_dataset$completeBCR == "Yes"])) > 1){
      clonotypes <- names(sort(table(tidy_dataset$igclonotypeID[tidy_dataset$clonotypeID == cltp & tidy_dataset$completeBCR == "Yes"]), decreasing = T))
      for(i in 1:length(clonotypes)){
        tidy_dataset$clonotypeID_new[tidy_dataset$clonotypeID == cltp & tidy_dataset$igclonotypeID == clonotypes[i]] <- paste0(cltp, "x", i)
      }
    }
  }
  # Add chain if no need of split
  tidy_dataset$clonotypeID_new[is.na(tidy_dataset$clonotypeID_new) & tidy_dataset$completeBCR == "Yes"] <- tidy_dataset$clonotypeID[is.na(tidy_dataset$clonotypeID_new) & tidy_dataset$completeBCR == "Yes"]

  tidy_dataset$clonotypeID_new[is.na(tidy_dataset$clonotypeID_new)] <- sapply(tidy_dataset$clonotypeVariantID_in_Cltp[is.na(tidy_dataset$clonotypeID_new)], function(cv){
    if(length(tidy_dataset$clonotypeID_new[tidy_dataset$completeBCR == "Yes" & tidy_dataset$clonotypeVariantID_in_Cltp == cv]) == 0){
      NA
    }else{
      dict_ct3 <- aggregate(x = tidy_dataset$clonotypeID_new[tidy_dataset$completeBCR == "Yes" & tidy_dataset$clonotypeVariantID_in_Cltp == cv],
                            by = list(tidy_dataset$clonotypeID_new[tidy_dataset$completeBCR == "Yes" & tidy_dataset$clonotypeVariantID_in_Cltp == cv]), FUN = length)
      dict_ct3 <- dict_ct3[order(dict_ct3$x, decreasing = T),]
      dict_ct3$Group.1[1]
    }
  })

  # If still NA... add chain based on CV like it was an independent clonotype (from number of current clonotypes to nClones)
  for(row in 1:nrow(tidy_dataset)){

    if(!is.na(tidy_dataset$clonotypeID_new[row])){next}

    cltp <- tidy_dataset$clonotypeID[row]

    newCltps <- unique(tidy_dataset$clonotypeID_new[!is.na(tidy_dataset$clonotypeID_new) & tidy_dataset$clonotypeID == cltp])

    if(length(newCltps) == 0){
      tidy_dataset$clonotypeID_new[row] <- cltp
    }else if(length(newCltps) == 1){
      if(tidy_dataset$completeBCR[row] == "Single_chain_1"){
        tidy_dataset$clonotypeID_new[!is.na(tidy_dataset$clonotypeID_new) & tidy_dataset$clonotypeID == cltp] <- paste0(cltp, "x1")
        tidy_dataset$clonotypeID_new[row] <- paste0(cltp, "x2")
      }else{
        tidy_dataset$clonotypeID_new[row] <- cltp
      }
    }else{
      tidy_dataset$clonotypeID_new[row] <- paste0(cltp, "x", (length(newCltps)+1))
    }
  }
  return(tidy_dataset$clonotypeID_new)
}

## Function for annotating immunogenetic data
.assign_CLL_subsets <- function(functionality, VDJ, cdr3aa, annotate_satellite_subsets){
  subset <- "NA"
  vGene <- unname(sapply(strsplit(strsplit(VDJ, split = "/")[[1]][1], ",")[[1]], function(x) strsplit(x, "\\*")[[1]][1]))

  if(functionality == "productive" & all(grepl("IGH", vGene)) & cdr3aa != "" & !is.na(cdr3aa)){
    subset <- "Unassigned"
    cdr3aa <- paste0(strsplit(cdr3aa, split = "")[[1]][2:(nchar(cdr3aa)-1)], collapse = "")
    cdr3len <- nchar(cdr3aa)

    if(any(grepl("IGHV[157]", vGene)) & cdr3len == 13 & grepl("^AR.QWL....FDY$", cdr3aa)){ subset <- "CLL#1"
    } else if("IGHV3-21" %in% vGene & cdr3len == 9 & grepl("^A.[DE]...MDV$", cdr3aa)){ subset <- "CLL#2"
    } else if("IGHV4-34" %in% vGene & cdr3len == 20 & grepl("^[AV]RG.......[RK]RYYYYGMDV$", cdr3aa)){ subset <- "CLL#4"
    } else if("IGHV1-69" %in% vGene & cdr3len == 20 & grepl("^AR....GV[IV]...YYYY[GY]MDV$", cdr3aa)){ subset <- "CLL#5"
    } else if("IGHV1-69" %in% vGene & cdr3len == 21 & grepl("^ARGG.YDY[VI]WGSYR.NDAFDI$", cdr3aa)){ subset <- "CLL#6"
    } else if("IGHV4-39" %in% vGene & cdr3len == 19 & grepl("^A[RST]...YSSSWY...NWFDP$", cdr3aa)){ subset <- "CLL#8"
    } else if("IGHV4-39" %in% vGene & cdr3len == 18 & grepl("^A[RST]...YSSSWY...WFDP$", cdr3aa)){ subset <- "CLL#8B"
    } else if("IGHV4-39" %in% vGene & cdr3len == 22 & grepl("^AR[HD]R.GYCSSTSCYYYYYGMDV$", cdr3aa)){ subset <- "CLL#10"
    } else if(any(grepl("IGHV1-(2|46)$", vGene)) & cdr3len == 19 & grepl("^ARD..YYDSSGYY[ST]..FDY$", cdr3aa)){ subset <- "CLL#12"
    } else if(any(grepl("IGHV[246]", vGene)) & cdr3len == 10 & grepl("^[AV]RGG.W.FD.$", cdr3aa)){ subset <- "CLL#14"
    } else if("IGHV4-34" %in% vGene & cdr3len == 24 & grepl("^A.RFYCSG..C....YYYYYG[LM]D[VA]$", cdr3aa)){ subset <- "CLL#16"
    } else if("IGHV3-48" %in% vGene & cdr3len == 21 & grepl("^AR[DE].DFWSGYY.YYYYY[GY]MDV$", cdr3aa)){ subset <- "CLL#31"
    } else if(any(grepl("IGHV1-(58|69)$", vGene)) & cdr3len == 12 & grepl("^A...DFWSGY..$", cdr3aa)){ subset <- "CLL#59"
    } else if(any(grepl("IGHV3", vGene)) & cdr3len == 21 & grepl("^A[KR][DE][ST][PL]LVV[PV][AT]AI[FY]YYYYGMDV$", cdr3aa)){ subset <- "CLL#64B"
    } else if(any(grepl("IGHV3", vGene)) & cdr3len == 12 & grepl("^A[KR]D....[WY]..DY$", cdr3aa)){ subset <- "CLL#73"
    } else if(any(grepl("IGHV4-(4|59)$", vGene)) & cdr3len == 14 & grepl("^[AV]RG[PA][DN].[ST]GW..[FL].Y$", cdr3aa)){ subset <- "CLL#77"
    } else if(any(grepl("IGHV[157]", vGene)) & cdr3len == 14 & grepl("^AR.QWL.....FDY$", cdr3aa)){ subset <- "CLL#99"
    } else if("IGHV1-18" %in% vGene & cdr3len == 10 & grepl("^AR.SGG..[DE].$", cdr3aa)){ subset <- "CLL#111"
    } else if("IGHV3-48" %in% vGene & cdr3len == 9 & grepl("^AR[DE]......$", cdr3aa)){ subset <- "CLL#169"
    } else if("IGHV3-72" %in% vGene & cdr3len == 17 & grepl("^[AV]R..YC[ST][SG][TG][TS]CR..[FL]D.$", cdr3aa)){ subset <- "CLL#188"
    } else if("IGHV4-34" %in% vGene & cdr3len == 17 & grepl("^ARR...W.....D[AG]FD.$", cdr3aa)){ subset <- "CLL#201"
    } else if("IGHV3-30" %in% vGene & cdr3len == 19 & grepl("^AK[VI]...G.F....YYGMD[VA]$", cdr3aa)){ subset <- "CLL#252"
    } else if("IGHV2-5"  %in% vGene & cdr3len == 17 & grepl("^AHR......W..G.FDY$", cdr3aa)){ subset <- "CLL#148B"
    } else if("IGHV1-2"  %in% vGene & cdr3len == 17 & grepl("^AR.[YL]SGSYYYYYYGMDV$", cdr3aa)){ subset <- "CLL#28A"
    } else if("IGHV1-69" %in% vGene & cdr3len == 22 & grepl("^A...DIVVVPAA..YYYYGMDV$", cdr3aa)){ subset <- "CLL#3C2"
    } else if("IGHV1-69" %in% vGene & cdr3len == 22 & grepl("^AR..PDIVVVPAAI.[YR]YYGMDV$", cdr3aa)){ subset <- "CLL#3C3"
    } else if("IGHV1-69" %in% vGene & cdr3len == 23 & grepl("^A[RST]....DFWSGYYPNYYYYGMDV$", cdr3aa)){ subset <- "CLL#7C2"
    } else if("IGHV1-69" %in% vGene & cdr3len == 24 & grepl("^A.....[YGD]DFWSGYYPNYYYY[GY]MDV$", cdr3aa)){ subset <- "CLL#7D3"
    } else if(any(grepl("IGHV3", vGene)) & cdr3len == 14 & grepl("^ARG..GDY...FD[YIV]$", cdr3aa)){ subset <- "CLL#202" }

    if(subset == "Unassigned" & annotate_satellite_subsets){

      if(any(grepl("IGHV[157]", vGene)) & cdr3len %in% 11:16 & grepl("QWL", substr(cdr3aa, 2, 8)) ){ subset <- "Satellite#1/99"
      } else if(any(grepl("IGHV3", vGene)) & cdr3len %in% 7:11 & grepl("[DE]", substr(cdr3aa, 1, 5)) ){ subset <- "Satellite#2/169"
      } else if("IGHV4-34" %in% vGene & cdr3len %in% 18:22 & grepl("KR|RR", substr(cdr3aa, 9, 18)) ){ subset <- "Satellite#4"
      } else if(any(grepl("IGHV[157]", vGene)) & cdr3len %in% 18:22 & grepl("GVV|GVI", substr(cdr3aa, 5, 11)) ){ subset <- "Satellite#5"
      } else if(any(grepl("IGHV[157]", vGene)) & cdr3len %in% 19:23 & grepl("WGSYK|WGSYR", substr(cdr3aa, 8, 16)) ){ subset <- "Satellite#6"
      } else if(any(grepl("IGHV[246]", vGene)) & cdr3len %in% 16:21 & grepl("YSSSWY", substr(cdr3aa, 4, 13)) ){ subset <- "Satellite#8/8B"
      } else if(any(grepl("IGHV[246]", vGene)) & cdr3len %in% 20:24 & grepl("GYCSSTSC", substr(cdr3aa, 4, 15)) ){ subset <- "Satellite#10"
      } else if(any(grepl("IGHV[157]", vGene)) & cdr3len %in% 17:21 & grepl("YYDSSGYY", substr(cdr3aa, 4, 15)) ){ subset <- "Satellite#12"
      } else if(any(grepl("IGHV[246]", vGene)) & cdr3len %in% 8:12 & grepl("RGG", substr(cdr3aa, 1, 6)) ){ subset <- "Satellite#14"
      } else if(any(grepl("IGHV[246]", vGene)) & cdr3len %in% 22:26 & grepl("FYCS", substr(cdr3aa, 2, 9)) ){ subset <- "Satellite#16"
      } else if(any(grepl("IGHV3", vGene)) & cdr3len %in% 19:23 & grepl("FWSGY", substr(cdr3aa, 4, 12)) ){ subset <- "Satellite#31"
      } else if(any(grepl("IGHV[157]", vGene)) & cdr3len %in% 10:14 & grepl("FWSGY", substr(cdr3aa, 4, 11)) ){ subset <- "Satellite#59"
      } else if(any(grepl("IGHV3", vGene)) & cdr3len %in% 19:23 & grepl("LVV", substr(cdr3aa, 4, 10)) ){ subset <- "Satellite#64B"
      } else if(any(grepl("IGHV3", vGene)) & cdr3len %in% 10:14 & grepl("KD|RD", substr(cdr3aa, 1, 5)) &  grepl("DY", substr(cdr3aa, 9, 13)) ){ subset <- "Satellite#73"
      } else if(any(grepl("IGHV[246]", vGene)) & cdr3len %in% 12:16 & grepl("GW", substr(cdr3aa, 6, 11)) ){ subset <- "Satellite#77"
      } else if(any(grepl("IGHV[157]", vGene)) & cdr3len %in% 8:12 & grepl("SGG", substr(cdr3aa, 2, 8)) ){ subset <- "Satellite#111"
      } else if(any(grepl("IGHV3", vGene)) & cdr3len %in% 15:19 & grepl("YC|FC", substr(cdr3aa, 3, 8)) &  grepl("CR|CY", substr(cdr3aa, 9, 14)) ){ subset <- "Satellite#188"
      } else if(any(grepl("IGHV[246]", vGene)) & cdr3len %in% 15:19 & grepl("ARR", substr(cdr3aa, 1, 5)) &  grepl("W", substr(cdr3aa, 5, 9)) ){ subset <- "Satellite#201"
      } else if(any(grepl("IGHV3", vGene)) & cdr3len %in% 17:21 & grepl("G", substr(cdr3aa, 5, 9)) &  grepl("F", substr(cdr3aa, 7, 11)) ){ subset <- "Satellite#252"
      } else if(any(grepl("IGHV[246]", vGene)) & cdr3len %in% 15:19 & grepl("HR", substr(cdr3aa, 1, 5)) &  grepl("W", substr(cdr3aa, 8, 12)) ){ subset <- "Satellite#148B"
      } else if(any(grepl("IGHV[157]", vGene)) & cdr3len %in% 15:19 & grepl("SGS", substr(cdr3aa, 3, 9)) ){ subset <- "Satellite#28A"
      } else if(any(grepl("IGHV[157]", vGene)) & cdr3len %in% 20:24 & grepl("DIVVVPAA", substr(cdr3aa, 3, 15)) ){ subset <- "Satellite#3C2/3C3"
      } else if(any(grepl("IGHV[157]", vGene)) & cdr3len %in% 21:26 & grepl("DFWSGY", substr(cdr3aa, 5, 15)) ){ subset <- "Satellite#7C2/7D3"
      } else if(any(grepl("IGHV3", vGene)) & cdr3len %in% 12:16 & grepl("GDY", substr(cdr3aa, 4, 10)) ){ subset <- "Satellite#202" }
    }
  }
  return(subset)
}
.annotateR110mut <- function(VDJ, VDJseq){
  r110 <- ""
  vGene <- strsplit(strsplit(VDJ, split = "/")[[1]][1], split = "\\*")[[1]][1]
  if(vGene == "IGLV3-21" & endsWith(VDJseq, "C")){
    r110 <- "[R110]"
  }
  return(r110)
}
.annotateAGS <- function(functionality, VDJ, gene_pos, VDJ_aa_seq, Germline_nt_seq){
  ags <- ""
  vGene <- unname(sapply(strsplit(strsplit(VDJ, split = "/")[[1]][1], ",")[[1]], function(x) strsplit(x, "\\*")[[1]][1]))

  if(functionality == "productive" & all(grepl("IGH", vGene))){

    pos_v <- as.numeric(strsplit(gene_pos, split = "-")[[1]])%/%3
    fr1 <- pos_v[1]
    cdr1 <- pos_v[2]
    fr2 <- pos_v[3]
    cdr2 <- pos_v[4]
    fr3 <- pos_v[5]
    cdr3 <- pos_v[6]
    fr4 <- pos_v[7]

    motif <- "(?=(N[^P][TS]))"
    matches <- str_match_all(VDJ_aa_seq, motif)[[1]][,2]
    locs <- str_locate_all(VDJ_aa_seq, motif)[[1]][,1]

    matches_df <- data.frame(start = locs, match = matches)

    if(nrow(matches_df) > 0){

      matches_df[, c("location", "gene", "filter")] <- NA
      Germline_aa_seq <- .translate_sequence(Germline_nt_seq)

      for(i in 1:nrow(matches_df)){

        match_start <- matches_df$start[i]
        match_end <- match_start + 2

        if(substr(VDJ_aa_seq, match_start, match_end) == substr(Germline_aa_seq, match_start, match_end)){
          matches_df$filter[i] <- "OUT"
          next

        } else{
          matches_df$filter[i] <- "PASS"

          if(match_end <= fr1){
            matches_df$location[i] <- "FR"
            matches_df$gene[i] <- "FR1"

          } else if(match_start <= (fr1+cdr1)){
            matches_df$location[i] <- "CDR"
            matches_df$gene[i] <- "CDR1"

          } else if(match_end <= (fr1+cdr1+fr2)){
            matches_df$location[i] <- "FR"
            matches_df$gene[i] <- "FR2"

          } else if(match_start <= (fr1+cdr1+fr2+cdr2)){
            matches_df$location[i] <- "CDR"
            matches_df$gene[i] <- "CDR2"

          } else if(match_end <= (fr1+cdr1+fr2+cdr2+fr3)){
            matches_df$location[i] <- "FR"
            matches_df$gene[i] <- "FR3"

          } else if(match_start <= (fr1+cdr1+fr2+cdr2+fr3+cdr3)){
            matches_df$location[i] <- "CDR"
            matches_df$gene[i] <- "CDR3"

          } else {
            matches_df$location[i] <- "FR"
            matches_df$gene[i] <- "FR4"
          }
        }
      }
      matches_df <- matches_df[matches_df$filter == "PASS",]

      if(nrow(matches_df) > 0){
        agsType <- ifelse("CDR" %in% matches_df$location, "CDR-AGS", "FR-AGS")
        agsList <- paste(paste(matches_df$match, matches_df$gene, sep = ":"), collapse = ";")
        ags <- paste0("[", agsType, " ", "(", agsList, ")", "]")
      }
    }
  }
  return(ags)
}
.assign_CLL_subsets_to_Clonotypes <- function(functionality, VDJ, CDR3_regexp, annotate_satellite_subsets){
  cdr3_list <- .possibleCDR3(CDR3_regexp)
  Clonotype_subsets <- paste0(unique(sapply(cdr3_list, function(x) .assign_CLL_subsets(functionality, VDJ, paste0("X",x,"X"), annotate_satellite_subsets))), collapse = ",")
  return(Clonotype_subsets)
}
.calculate_v_ident <- function(seq, germ, len){
  seq_v <- unlist(strsplit(seq, split = ""))[1:len]
  germ_v <- unlist(strsplit(germ, split = ""))[1:len]
  v_iden <- (sum(seq_v == germ_v)/len)*100
  return(v_iden)
}
.consensus_germline <- function(clonotype_df){

  major_length <- as.numeric(names(sort(table(clonotype_df$CDR3len), decreasing = T)[1]))

  clonotype_df <- clonotype_df[clonotype_df$CDR3len == major_length,]

  clonotype_df$perc_ident <- apply(clonotype_df[, c("germline_alignment", "VDJseq")], MARGIN = 1, function(x){
    germline_df <- str_split_fixed(x[1], pattern = "", nchar(x[1]))
    somatic_df <- str_split_fixed(x[2], pattern = "", nchar(x[2]))
    perc_ident <- sum(somatic_df[germline_df != "N"] == germline_df[germline_df != "N"])/length(somatic_df[germline_df != "N"])
    perc_ident })

  clonotype_df$N_in_seq <- apply(clonotype_df[, c("germline_alignment", "VDJseq")], MARGIN = 1, function(x){
    germline_df <- str_split_fixed(x[1], pattern = "", nchar(x[1]))
    somatic_df <- str_split_fixed(x[2], pattern = "", nchar(x[2]))
    N_in_seq <- length(somatic_df[germline_df == "N"])
    N_in_seq })

  closest_to_germline <- clonotype_df[clonotype_df$perc_ident == max(clonotype_df$perc_ident),]
  closest_to_germline <- closest_to_germline[closest_to_germline$N_in_seq == min(closest_to_germline$N_in_seq),]

  con_germline_df <- str_split_fixed(closest_to_germline$germline_alignment[1], pattern = "", nchar(closest_to_germline$germline_alignment[1]))
  con_somatic_df <- str_split_fixed(closest_to_germline$VDJseq[1], pattern = "", nchar(closest_to_germline$VDJseq[1]))

  con_germline_df[con_germline_df == "N"] <- con_somatic_df[con_germline_df == "N"]
  con_germ <- paste0(con_germline_df, collapse = "")
  return(list(major_length, con_germ))
}
.correct_consensus_germline <- function(somatic_seq, germline_seq, consensus_germline_seq, correct_germline_subMat){

  aln <- pairwiseAlignment(pattern = consensus_germline_seq, subject = germline_seq, substitutionMatrix = correct_germline_subMat)

  pattern <- strsplit(as.character(aln@pattern), split = "")[[1]]
  subject <- strsplit(as.character(aln@subject), split = "")[[1]]
  subject_somatic <- strsplit(somatic_seq, split = "")[[1]]

  deletion_list <- unlist(apply(str_locate_all(as.character(aln@subject), "[ACTG]-+[ACTG]")[[1]], 1, list))
  insertion_list <- unlist(apply(str_locate_all(as.character(aln@pattern), "[ACTG]-+[ACTG]")[[1]], 1, list))

  indel_v <- c()
  if(length(deletion_list) > 0){
    for(start in seq(1, length(deletion_list), 2)){
      ref <- paste0(pattern[deletion_list[start]:(deletion_list[start+1]-1)], collapse = "")
      alt <- subject[deletion_list[start]]

      indel_v <- c(indel_v, paste(deletion_list[start], ref, alt, sep = "_"))
      subject[(deletion_list[start]+1):(deletion_list[start+1]-1)] <- pattern[(deletion_list[start]+1):(deletion_list[start+1]-1)]
      subject_somatic <- c(subject_somatic[1:deletion_list[start]], pattern[(deletion_list[start]+1):(deletion_list[start+1]-1)],  subject_somatic[(deletion_list[start]+1):length(subject_somatic)])
    }
  }

  if(length(insertion_list) > 0){
    removed_nt <- 0
    for(start in seq(1, length(insertion_list), 2)){
      indel_start <- insertion_list[start]-removed_nt
      indel_end <- insertion_list[start+1]-removed_nt

      ref <- pattern[indel_start]
      alt <- paste0(subject[indel_start:(indel_end-1)], collapse = "")
      indel_v <- c(indel_v, paste(indel_start, ref, alt, sep = "_"))

      subject <- paste0(subject[-((indel_start+1):(indel_end-1))], collapse = "")
      subject_somatic <- paste0(subject_somatic[-((indel_start+1):(indel_end-1))], collapse = "")
      removed_nt <- removed_nt + length((indel_start+1):(indel_end-1))
    }
  }
  return(list(paste(indel_v, collapse = ","), paste0(subject_somatic, collapse = "")))
}

## Functions for building consensus CDR3 sequences
.consensus_CDR3 <- function(CDR3_list){

  major_length <- as.numeric(names(sort(table(nchar(CDR3_list)), decreasing = T)[1]))
  CDR3_list <- CDR3_list[nchar(CDR3_list) == major_length]

  con_cdr3 <- ""
  cdr3_mat <- matrix(unlist(strsplit(CDR3_list, "")), nrow = length(CDR3_list), byrow = TRUE)
  for(col in 1:ncol(cdr3_mat)){
    rank <- table(cdr3_mat[,col])
    mainLetter <- names(rank[which(rank == max(rank))])
    if(length(mainLetter) > 1){mainLetter <- paste0("[",paste(mainLetter, collapse = ""),"]")}
    con_cdr3 <- paste0(con_cdr3, mainLetter)
  }
  return(con_cdr3)
}
.possibleCDR3 <- function(CDR3_regexp){
  cdr3_list <- c()
  if(!grepl("\\[", CDR3_regexp)){return(CDR3_regexp)}
  sp <- strsplit(CDR3_regexp, split = "")[[1]]
  start <- grep("\\[", sp)[1]
  end <- grep("\\]", sp)[1]

  possibleLetter <- sp[(start+1):(end-1)]
  for(let in possibleLetter){
    nextCDR3 <- sp[-((start+1):end)]
    nextCDR3[start] <- let
    cdr3_list <- c(cdr3_list, .possibleCDR3(paste0(nextCDR3, collapse = "")))
  }
  return(cdr3_list)
}

## Function to rename the final list
.name_final_list <- function(all_out_list){
  list_names <- c()
  for(i in 1:length(all_out_list)){
    name <- unique(all_out_list[[i]]$SampleID)
    list_names <- c(list_names, paste0(name,"_annot"))
  }
  names(all_out_list) <- list_names
  return(all_out_list)
}

