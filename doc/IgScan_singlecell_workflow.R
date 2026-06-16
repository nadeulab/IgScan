## ----echo=FALSE---------------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(devtools)
load_all("~/Desktop/IgScan/")
library(Seurat)
library(ggplot2)
library(dplyr)

## ----igblast_reannot, message=FALSE, warning=FALSE, tidy=FALSE----------------

s1 <- system.file("extdata/igscan_test_10xBCR_sample1.fasta", package = "IgScan", mustWork = T)
s2 <- system.file("extdata/igscan_test_10xBCR_sample2.fasta", package = "IgScan", mustWork = T)

IgScan::Run_IgBlast_from_RawData(
  sample_paths = c(s1, s2), 
  sample_labels = c("Sample1", "Sample2"), 
  input_format = "10xbcr_fasta", 
  data_type = "single_cell", 
  annotate_C = TRUE, 
  outputDir = "~/Desktop/TEST_IGSCAN_VINGETTE/", ## Change by your own outputDir
  run_IgBlast_report = TRUE)

## ----igscan_workflow, tidy = FALSE, message=FALSE, warning=FALSE--------------
igscan_out <- IgScan::Run_IgScan_Annotation(
  sample_labels = c("Sample1", "Sample2"),
  case_labels = c("CaseA", "CaseA"),
  outputDir = "~/Desktop/TEST_IGSCAN_VINGETTE/", 
  input_format = "10xbcr_fasta", 
  analysis_mode = "joint",
  data_type = "single_cell",
  material_type = "rna",
  v_primer = "full_length", 
  remove_tmp = FALSE, 
  hc_similarity_cutoff = 0.2, 
  hc_mode = "average", 
  cdr3_mode = "nt", 
  cdr3_InDel_correction_mode = "soft_filter")

## ----igscan_combine_sc, tidy = FALSE, message=FALSE, warning=FALSE------------
o1 <- system.file("extdata/igscan_test_10xSeurat_sample1.rds", 
                  package = "IgScan", mustWork = T)

o2 <- system.file("extdata/igscan_test_10xSeurat_sample2.rds", 
                  package = "IgScan", mustWork = T)

seurat_1 <- readRDS(o1)
seurat_2 <- readRDS(o2)

seurat_1 <- IgScan::combine_IgScan_Seurat(igscan_out = igscan_out$Sample1_annot, 
                                          seurat_object = seurat_1)

seurat_2 <- IgScan::combine_IgScan_Seurat(igscan_out = igscan_out$Sample2_annot, 
                                          seurat_object = seurat_2)

## ----color_pallete, tidy = FALSE, message=FALSE, warning=FALSE----------------
color <- c("C1" = "#1B9E77", "C2" = "#D95F09", "C3" = "#7570B9", "C4" = "#E6AB02")

## ----igscan_anal_1, tidy = FALSE, message=FALSE, fig.height=5, fig.width=6, fig.align='center'----
seurat_merge <- merge(seurat_1, seurat_2)

clonotype_df <- seurat_merge@meta.data %>%
  group_by(orig.ident, igClonotypeID_num) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(freq_rel = count / sum(count)*100)

freq_plt <- ggplot(clonotype_df, aes(x = orig.ident, y = freq_rel, fill = igClonotypeID_num)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color) + theme_bw(base_size = 15) +
  labs(x = NULL, y = "Clonotype frequency (%)", fill = "ClonotypeID")

freq_plt

## ----igscan_anal_2, message=FALSE, warning=FALSE------------------------------
seurat_merge <- JoinLayers(seurat_merge)
seurat_merge <- NormalizeData(object = seurat_merge, verbose = F)
seurat_merge <- FindVariableFeatures(object = seurat_merge, verbose = F)
seurat_merge <- ScaleData(object = seurat_merge, verbose = F)
seurat_merge <- RunPCA(object = seurat_merge, verbose = F)
seurat_merge <- FindNeighbors(object = seurat_merge, dims = 1:30, verbose = F)
seurat_merge <- FindClusters(object = seurat_merge, resolution = 1, verbose = F)
seurat_merge <- RunUMAP(object = seurat_merge, dims = 1:30, verbose = F)

sam_dim <- DimPlot(seurat_merge, group.by = "orig.ident",
                   cols = c("dodgerblue", "darkred"), pt.size = 2) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2") + ggtitle("SampleID")

clone_dim <- DimPlot(seurat_merge, group.by = "igClonotypeID_num",
                     cols = color, pt.size = 2) +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2") + ggtitle("ClonotypeID")

## ----igscan_anal_3, fig.height=3, fig.width=7, fig.align='center'-------------
print(sam_dim + clone_dim)

## ----igscan_anal_4, message=FALSE, warning=FALSE------------------------------
table(seurat_merge$igClonotypeID_num, seurat_merge$completeBCR, useNA = "ifany")

## ----igscan_anal_5, message=FALSE, warning=FALSE------------------------------
seurat_merge <- subset(seurat_merge, 
                       subset = !completeBCR %in% c("Clonotype_doublet", "Repeated_chain"))

## ----igscan_anal_6, message=FALSE, warning=FALSE------------------------------
seurat_merge$case <- "Case1"
seurat_merge <- IgScan::rescue_single_chain_cells(single_cell_object = seurat_merge, 
                                                  group_col = "case")

## ----igscan_anal_7, message=FALSE, warning=FALSE------------------------------
table(seurat_merge$igClonotypeID_num, seurat_merge$completeBCR, useNA = "ifany")

## ----igscan_anal_8, message=FALSE, warning=FALSE------------------------------
seurat_merge_filter <- subset(seurat_merge, subset = igClonotypeID_num != "C1")
table(seurat_merge_filter$igClonotypeID_num, seurat_merge_filter$orig.ident)

## ----igscan_anal_9, message=FALSE, warning=FALSE------------------------------
seurat_merge_filter <- recalculate_IDs_single_cell(single_cell_object = seurat_merge_filter,
                                                   group_col = "case")
table(seurat_merge_filter$igClonotypeID_num, seurat_merge_filter$orig.ident)

## ----igscan_anal_10, tidy = FALSE, message=FALSE, fig.height=5, fig.width=6, fig.align='center'----
clonotype_df2 <- seurat_merge_filter@meta.data %>%
  group_by(orig.ident, igClonotypeID_num) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  group_by(orig.ident) %>% 
  mutate(freq_rel = count / sum(count)*100)

freq_plt2 <- ggplot(clonotype_df2, 
                   aes(x = orig.ident, y = freq_rel, fill = igClonotypeID_num)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = color) + theme_bw(base_size = 15) + 
  labs(x = NULL, y = "Clonotype frequency (%)", fill = "ClonotypeID")

freq_plt2

