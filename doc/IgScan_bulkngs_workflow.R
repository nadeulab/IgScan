## ----echo=FALSE---------------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(devtools)
load_all("~/Desktop/IgScan/")
library(ggplot2)
library(dplyr)
library(dowser)

## ----igblast_reannot, message=FALSE, warning=FALSE, tidy=FALSE----------------

s1 <- system.file("extdata/igscan_test_bulkNGS_sample1.fasta", package = "IgScan", mustWork = T)
s2 <- system.file("extdata/igscan_test_bulkNGS_sample2.fasta", package = "IgScan", mustWork = T)

IgScan::Run_IgBlast_from_RawData(
  sample_paths = c(s1, s2), 
  sample_labels = c("Sample1", "Sample2"), 
  input_format = "fasta", 
  data_type = "bulk", 
  annotate_C = TRUE, 
  outputDir = "~/Desktop/TEST_IGSCAN_VINGETTE_NGS/", ## Change by your own outputDir
  run_IgBlast_report = TRUE)

## ----igscan_workflow, tidy = FALSE, message=FALSE, warning=FALSE--------------
igscan_out <- IgScan::Run_IgScan_Annotation(
  sample_labels = c("Sample1", "Sample2"),
  outputDir = "~/Desktop/TEST_IGSCAN_VINGETTE_NGS/", 
  input_format = "fasta", 
  analysis_mode = "single",
  data_type = "bulk", 
  min_reads = 2,
  material_type = "dna",
  v_primer = "full_length", 
  remove_tmp = FALSE, 
  hc_similarity_cutoff = 0.2, 
  hc_mode = "average", 
  cdr3_mode = "nt", 
  cdr3_InDel_correction_mode = "soft_filter")

## ----remove_cont, tidy = FALSE, message=FALSE, warning=FALSE------------------
merge_igscan_out <- do.call(rbind, igscan_out)

merge_igscan_out <- IgScan::filter_BCR_contam_bulk(igscan_data_frame = merge_igscan_out,
                                                   remove_contamination = FALSE)

merge_igscan_out[, c("SampleID", "VDJ_genes", "Junction_aa", 
                     "Contamination_FLAG", "Contamination_Freq", "Contamination_Sample")]

## ----recalc_cont, tidy = FALSE, message=FALSE, warning=FALSE------------------
merge_igscan_out <- merge_igscan_out[merge_igscan_out$Contamination_FLAG == "PASS",]    
merge_igscan_out <- recalculate_IDs_bulk(merge_igscan_out, group_col = "SampleID")

## ----export_airr, tidy = FALSE, message=FALSE, warning=FALSE------------------
airr_object <- export_AIRR_format(object = merge_igscan_out, 
                                  dir = "~/Desktop/TEST_IGSCAN_VINGETTE_NGS/", 
                                  fileName = "exportAIRR_testNGS.tsv", 
                                  germline_aln = "masked")

## ----airr_gen, tidy = FALSE, eval=FALSE---------------------------------------
# airr_object <- export_AIRR_format(object = IgScan_Df,
#                                   dir = OutputDir,
#                                   fileName = FileName,
#                                   germline_aln = "consensus",
#                                   metadata = "SubcloneID_in_Clonotype") ## Added metadata

## ----dowser_interact, tidy = FALSE, message=FALSE, warning=FALSE--------------
## Load example AIRR file
airr_path <- system.file("extdata/igscan_test_AIRR_file_sample3.rds", package = "IgScan", mustWork = T)
example_airr <- readRDS(airr_path)

## Format clones for phylogenetic inference
clones <- formatClones(example_airr, germ = "germline_alignment", 
                       traits = "SubcloneID_in_Clonotype")

# Build maxmimum parsimony trees
trees <- getTrees(clones)
plots <- plotTrees(trees, tipsize=2, tips = "SubcloneID_in_Clonotype")

## ----plot_tree, tidy = FALSE, message=FALSE, warning=FALSE, fig.height=4, fig.width=7----
plots[[1]]

