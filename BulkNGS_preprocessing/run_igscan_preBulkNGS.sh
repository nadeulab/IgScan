#!/bin/bash

# ============================
# Preprocessing of raw bulk NGS data for the IgScan workflow
# ============================
# Use:
# IgScan_directory/tools/igscan_preprocess_bulkngs.sh sample_R1.fastq sample_R2.fastq output.fasta [min_length]

# ============================
# 1. Input arguments
# ============================
PrePro_DIR="$(cd "$(dirname "$(realpath "$0")")" && pwd)"
R1=$1
R2=$2
OUTPUT=$3
MIN_JOINT_LENGTH=${4:-250}  # Default to 250 if not provided

if [ -z "$R1" ] || [ -z "$R2" ] || [ -z "$OUTPUT" ]; then
  echo "Use: $0 <R1.fastq> <R2.fastq> <output.fasta> [min_length]"
  exit 1
fi

# ============================
# 2. Detect operating system (macOS or Linux)
# ============================
OS=$(uname)
if [[ "$OS" == "Linux" ]]; then
  FLASH="$PrePro_DIR/FLASH_1.2.11_Linux/flash"
elif [[ "$OS" == "Darwin" ]]; then
  FLASH="$PrePro_DIR/FLASH_1.2.11_MacOS/flash"
else
  echo "Non compatible operating system (expecting Linux or macOS): $OS"
  exit 1
fi

# ============================
# 3. Execute FLASH to join overlapping paired-end reads
# ============================
mkdir -p flash_output
$FLASH "$R1" "$R2" -m 10 -M 200 -x 0.1 -o merged -d ./flash_output

# ============================
# 4. Convert FASTQ to FASTA
# ============================
awk 'NR%4==1 {printf(">%s\n", substr($0,2))} NR%4==2 {print}' ./flash_output/merged.extendedFrags.fastq > merged.fasta

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
