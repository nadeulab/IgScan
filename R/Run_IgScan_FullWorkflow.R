#' Run IgScan Full Workflow
#'
#' @description This function provides a streamlined pipeline to process raw
#' sequencing data, perform IgBLAST reannotation, and run IgScan to generate
#' advanced immunogenetic analyses. It combines the functionalities of
#' `Run_IgBlast_from_RawData` and `Run_IgScan_Annotation` into a single step,
#' allowing users to process raw scRNA-seq/bulk NGS data, annotate immunoglobulin
#' sequences, and analyze clonotypes in a unified workflow. In addition, the workflow is
#' also compatible with single-cell V(D)J data generated using the Mission Bio platform,
#' provided that the IgScan input has been produced using our dedicated raw-data
#' preprocessing pipeline (see IgScan GitHub).
#'
#' It supports data from multiple sequencing platforms and file types:
#' \itemize{
#'   \item 10x BCR = "filtered_contig_annotations.csv / filtered_contig.fa"
#'   \item Parse BCR = "bcr_annotation_airr.tsv"
#'   \item BD Rhapsody BCR = "Contigs_AIRR.tsv"
#'   \item MiXCR = "clonotypes.IGX.txt"
#'   \item TRUST4 = "barcode_report.tsv"
#'   \item AIRR = "airr_rearrangement.tsv"
#'   \item IMGT AIRR = "vquest_airr.tsv"
#'   \item fasta
#' }
#'
#' @param sample_paths A vector of paths to input files or a vector of
#' directories containing the input file. Supported formats depend on the
#' `input_format` parameter.
#' @param sample_labels A vector of sample labels corresponding to `sample_paths`.
#'   Must have the same length as `sample_paths`.
#' @param analysis_mode Defines the mode of analysis to be performed. Options are
#'   "single" and "joint". In "single" mode, each sample is annotated independently.
#'   In "joint" mode, samples are grouped by case based on the "case_list" vector
#'   for a case-level joint annotation. Default is 'single'.
#' @param case_labels A vector of case labels matching the indexes of `sample_labels`.
#'   Required when analysis mode is set to 'joint'. Default is NULL.
#' @param input_format A string specifying the format of the input data, currently supporting:
#'   '10xBCR_fasta', '10xBCR_csv', 'ParseBCR', 'BDRhapsodyBCR', 'MiXCR',
#'   'TRUST4', 'AIRR', 'IMGT_AIRR' and 'fasta'.
#' @param data_type The type of data. Options: 'single_cell', 'bulk' or 'missionbio'. Default is 'single_cell'.
#' @param material_type The biological source of material. Options: 'dna' and 'rna'. Note that for
#'   matrial_type='rna', unproductive sequences are not expected, and will be directly removed. Default is 'rna'.
#' @param v_primer The primer sequence used for the V-region amplification.
#'   Options: 'full_length', 'fr1', 'fr2' and 'fr3'. Note that sequences with unexpected length pattern based on
#'   the chosen primer will be directly excluded from the analysis. Default is 'full_length'.
#' @param annotate_C Logical. If `TRUE`, the function performs an IgBlast annotation of
#'   constant regions using the NCBI IG C region database.
#' @param Evalue_cutoff A numeric value specifying the E-value cutoff for IgBLAST alignment.
#'   Default is `NULL` (no cutoff applied).
#' @param threads_IgBlast An integer specifying the number of threads for running IgBLAST.
#'   Default is 1.
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
#' @details
#' This function performs the following steps:
#'
#'    1. Validates the input data and format.
#'
#'    2. Converts input data to FASTA format if needed.
#'
#'    3. Runs IgBLAST for V(D)J-(C) reannotation and alignment.
#'
#'    4. Processes IgBLAST output with the IgScan workflow to perform an advanced
#'    immunogenetic annotation of the data, cluster sequences into clonotypes and
#'    formats the output based on the type of experiment (bulk or single cell).
#'
#' @return A list containing an IgScan annotation dataframe for every sample analyzed, which are also
#' saved in the output directory. These dataframes can be directly passed to further IgScan functions.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Run the complete pipeline on FASTA inputs
#' results <- Run_IgScan_FullWorkflow(
#'   sample_paths = list("path/to/sample1.fasta", "path/to/sample2.fasta"),
#'   sample_labels = c("Sample1", "Sample2"),
#'   analysis_mode = "single",
#'   Evalue_cutoff = 100,
#'   annotate_C = T,
#'   threads_IgBlast = 2,
#'   input_format = "fasta",
#'   material_type = "dna",
#'   v_primer = "full_length",
#'   data_type = "bulk",
#'   min_reads = 2,
#'   remove_tmp = TRUE,
#'   outputDir = "path/output/dir/",
#'   hc_similarity_cutoff = 0.2,
#'   hc_mode = "average",
#'   cdr3_mode = "nt",
#'   cdr3_InDel_correction_mode = "soft_filter",
#'   annotate_CLL_immGen = TRUE,
#'   annotate_satellite_subsets = FALSE,
#'   threads = 4)
#' }
#'
Run_IgScan_FullWorkflow <- function(sample_paths, sample_labels, Evalue_cutoff = NULL, annotate_C = TRUE, threads = 1, threads_IgBlast = 1, case_labels = NULL, input_format, analysis_mode = "single", material_type = "rna", v_primer = "full_length", data_type = "single_cell", min_reads = 2, remove_tmp = TRUE, outputDir = NULL, hc_similarity_cutoff = 0.2, hc_mode = "average", cdr3_mode = "nt", cdr3_InDel_correction_mode = "soft_filter", annotate_CLL_immGen = FALSE, annotate_satellite_subsets = TRUE, annotate_ags = FALSE, rescue_single_chain = FALSE, relaxed_rescue = FALSE){

  if(is.null(outputDir)){
    outputDir <- paste0("./", format(Sys.time(), "%d-%m-%Y_%H%M%S-IgScanResults/"))
    n <- 1
    while(TRUE){
      if(dir.exists(outputDir)){
        outputDir <- gsub("-IgScanResults.*/", paste0("-IgScanResults_",n,"/"), outputDir)
        n <- n+1
      } else{break}
    }
    system(paste0("mkdir ",outputDir))
    warning(paste0("Output directory has been set to ", outputDir), call. = FALSE)

  } else{
    if(!endsWith(outputDir, "/")){outputDir <- paste0(outputDir, "/")}
    if(!dir.exists(outputDir)){system(paste0("mkdir ",outputDir))}
  }

  summary_file <- paste0(outputDir, format(Sys.time(), "%d-%m-%Y_%H-%M-%S_IgScan_reporting_summary.log"))

  write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - Running IgBlast..."), file = summary_file)
  message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - Running IgBlast..."))
  Run_IgBlast_from_RawData(sample_paths = sample_paths, sample_labels = sample_labels, input_format = input_format, data_type = data_type, Evalue_cutoff = Evalue_cutoff, annotate_C = annotate_C, threads_IgBlast = threads_IgBlast, outputDir = outputDir, threads = threads, run_IgBlast_report = F)

  write(x = paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - Starting IgScan annotation..."), file = summary_file, append = T)
  message(paste0("[", format(Sys.time(), "%d-%m-%Y %H:%M:%S"), "] - Starting IgScan annotation..."))
  annotated_df_list <- Run_IgScan_Annotation(sample_labels = sample_labels, case_labels = case_labels, input_format = input_format, outputDir = outputDir, analysis_mode = analysis_mode, material_type = material_type, v_primer = v_primer, data_type = data_type, min_reads = min_reads, remove_tmp = remove_tmp, hc_similarity_cutoff = hc_similarity_cutoff, hc_mode = hc_mode, cdr3_mode = cdr3_mode, cdr3_InDel_correction_mode = cdr3_InDel_correction_mode, annotate_CLL_immGen = annotate_CLL_immGen, annotate_satellite_subsets = annotate_satellite_subsets, annotate_ags = annotate_ags, rescue_single_chain = rescue_single_chain, relaxed_rescue = relaxed_rescue, summary_file = summary_file, threads = threads)

  return(annotated_df_list)
}
