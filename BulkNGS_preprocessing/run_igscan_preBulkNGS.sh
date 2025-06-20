#!/bin/bash

# ============================
# Preprocessing of raw bulk NGS data for the IgScan workflow
# ============================
# Use:
# BulkNGS_preprocessing/igscan_preprocess_bulkngs.sh sample_R1.fastq sample_R2.fastq output.fasta [min_length] [min_overlap] [max_overlap] [perc_mismatch] [concat_unmerged]

# ============================
# 1. Input arguments
# ============================
PREPRO_DIR="$(cd "$(dirname "$(realpath "$0")")" && pwd)"
R1=$1
R2=$2
OUTPUT=$3
MIN_JOINT_LENGTH=${4:-250}  # Default to 250 if not provided

MIN_OVERLAP=${5:-10}
MAX_OVERLAP=${6:-200}
PERC_MISMATCH=${7:-0.1}

CONCAT_UNMERGED=${8:-"no"}

if [ -z "$R1" ] || [ -z "$R2" ] || [ -z "$OUTPUT" ]; then
  echo "Error: Missing required arguments."
  echo "Use: $0 <R1.fastq> <R2.fastq> <output.fasta> [min_length] [min_overlap] [max_overlap] [perc_mismatch] [concat_unmerged]"
  exit 1
fi

# ============================
# 2. Detect operating system (macOS or Linux)
# ============================
OS=$(uname)
if [[ "$OS" == "Linux" ]]; then
  FLASH="$PREPRO_DIR/FLASH_1.2.11_Linux/flash"
elif [[ "$OS" == "Darwin" ]]; then
  FLASH="$PREPRO_DIR/FLASH_1.2.11_MacOS/flash"
else
  echo "Non compatible operating system (expecting Linux or macOS): $OS"
  exit 1
fi

# ============================
# 3. Execute FLASH to join overlapping paired-end reads
# ============================
mkdir -p flash_output
$FLASH "$R1" "$R2" -m "$MIN_OVERLAP" -M "$MAX_OVERLAP" -x "$PERC_MISMATCH" -o merged -d ./flash_output

if [[ $? -ne 0 ]]; then
  echo "FLASH failed. Exiting."
  exit 1
fi

# ============================
# 4. Convert FASTQ to FASTA
# ============================
awk 'NR%4==1 {printf(">%s\n", substr($0,2))} NR%4==2 {print}' ./flash_output/merged.extendedFrags.fastq > merged.fasta

if [[ "$CONCAT_UNMERGED" == "yes" ]]; then
awk 'NR%4==1 {printf(">%s\n", substr($0,2))} NR%4==2 {print}' ./flash_output/merged.notCombined_1.fastq > ./flash_output/notCombined_1.fasta
awk 'NR%4==1 {printf(">%s\n", substr($0,2))} NR%4==2 {print}' ./flash_output/merged.notCombined_2.fastq > ./flash_output/notCombined_2.fasta

## Make the reverse complement of the R2 sequences
awk '
function revcomp(seq,  i, out) {
    gsub("A", "t", seq)
    gsub("T", "a", seq)
    gsub("G", "c", seq)
    gsub("C", "g", seq)
    out=""
    for (i=length(seq); i>0; i--) {
        out = out substr(seq, i, 1)
    }
    return toupper(out)
}
BEGIN { seq = "" }
{
    if ($0 ~ /^>/) {
        if (seq != "") {
            print revcomp(seq)
            seq = ""
        }
    } else {
        seq = seq $0
    }
}
END {
    if (seq != "") print revcomp(seq)
}
' ./flash_output/notCombined_2.fasta > ./flash_output/notCombined_2_RevComp.fasta

## Paste R1 sequences with revcomp(R2) into the same file and append to merged.fasta
paste - - < ./flash_output/notCombined_1.fasta | awk -v revfile=./flash_output/notCombined_2_RevComp.fasta '
BEGIN {
    while ((getline r < revfile) > 0) {
        revs[++n] = r
    }
    close(revfile)
}
{
    print $1
    print $2 revs[NR]
}
' >> merged.fasta
fi

# ============================
# 5. Collapse and sort sequences by abundance filtering by length
# ============================
awk -v min_len="$MIN_JOINT_LENGTH" '
  /^>/ {next}
  length($0) >= min_len {counts[$0]++}
  END {
    for (seq in counts) {
      printf("%d\t%s\n", counts[seq], seq);
    }
  }' merged.fasta | sort -k1,1nr > sorted.tmp

# ============================
# 6. Output in FASTA with readX_n=COUNT headers
# ============================
i=1
> "$OUTPUT"  # clear output file
while IFS=$'\t' read -r count seq; do
  echo ">read${i}_n=${count}" >> "$OUTPUT"
  echo "$seq" >> "$OUTPUT"
  ((i++))
done < sorted.tmp

# ============================
# 7. Cleanup
# ============================
rm -rf flash_output merged.fasta sorted.tmp

echo "IgScan BulkNGS pre-processing completed!"
echo "Fasta file generated and sorted by abundance: $OUTPUT"
