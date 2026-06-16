# Global variables used in IgScan
# This file is used to satisfy R CMD check NOTES about "no visible binding"

utils::globalVariables(c(
  
  # =========================
  # Core clonotype / VDJ data
  # =========================
  "clonotypeID",
  "Clonotype",
  "ClonotypeID",
  "ClonotypeVariantID",
  "clonotypeVariantID_in_Cltp",
  "ClonotypeVariant_nReads",
  "clonotypeVariant_nreads",
  "Subclone_nReads",
  "Pre_Clonotype",
  
  # =========================
  # Sequence / immunology
  # =========================
  "CDR3",
  "CDR3aa",
  "productive",
  "igVDJ_sequence",
  "igInDels",
  "igSubcloneID_all",
  "igSubcloneID_in_Clonotype_num",
  
  # =========================
  # Read / counts / metrics
  # =========================
  "n_reads",
  "total_reads_unique_seq",
  "Freq",
  "relFreq",
  "Score_Value",
  
  # =========================
  # Generic data.table / dplyr vars
  # =========================
  "Var1",
  "Var2",
  "Locus",
  "locus",
  "Chain_SbcID",
  "analysis_mode",
  "identity1",
  "outputDir",
  
  # =========================
  # Internal temporary variables
  # =========================
  "tmp",
  "tmp_col",
  
  # =========================
  # AIRR conversion fields
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
  # Dictionary / annotation objects
  # =========================
  "colnames_dictionary",
  
  # =========================
  # Misc plotting / stats
  # =========================
  "findpeaks"
))