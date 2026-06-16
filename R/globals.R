# Global variables to satisfy R CMD check NOTES (NSE / dplyr / data.table)

utils::globalVariables(c(

  # =========================
  # Core clonotypes / immunology
  # =========================
  "clonotypeID",
  "Clonotype",
  "ClonotypeID",
  "ClonalID",
  "ClonotypeVariantID",
  "ClonotypeVariant_nReads",
  "clonotypeVariant_nreads",
  "Subclone_nReads",
  "Pre_Clonotype",
  "Unique_SequenceID",

  # =========================
  # BCR / VDJ data
  # =========================
  "CDR3",
  "CDR3aa",
  "productive",
  "completeBCR",
  "igVDJ_sequence",
  "igInDels",
  "igSubcloneID_all",
  "igSubcloneID_in_Clonotype_num",

  # =========================
  # Counts / stats
  # =========================
  "n_reads",
  "total_reads_unique_seq",
  "Freq",
  "relFreq",
  "Score_Value",

  # =========================
  # Column names used in data.frames
  # =========================
  "Var1",
  "Var2",
  "Locus",
  "locus",
  "Chain_SbcID",
  "analysis_mode",
  "identity1",
  "outputDir",
  "Contamination_FLAG",

  # =========================
  # Temporary / dplyr NSE symbols
  # =========================
  "x",
  ".",
  "tmp",
  "tmp_col",

  # =========================
  # AIRR / annotation fields
  # =========================
  "sequence_alignment",
  "sequence_alignment_aa",
  "cdr1_start",
  "cdr1_end",
  "cdr2_start",
  "cdr2_end",
  "cdr3_start",
  "cdr3_end",
  "fwr1_start",
  "fwr1_end",
  "fwr2_start",
  "fwr2_end",
  "fwr3_start",
  "fwr3_end",
  "fwr4_start",
  "fwr4_end",

  # =========================
  # Dictionary / metadata
  # =========================
  "colnames_dictionary",

  # =========================
  # External / stats functions used non-namespace
  # =========================
  "aggregate",
  "reshape",
  "ave",
  "as",
  "as.dist",
  "cutree",
  "hclust",
  "setNames",
  "combn",
  "capture.output",
  "write.table",
  "pdf",
  "dev.off",
  "findpeaks",

  # =========================
  # Seurat / SCE slots
  # =========================
  "DataFrame",
  "colData<-"
))
