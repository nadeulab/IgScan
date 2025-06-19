---
editor_options: 
  markdown: 
    wrap: 72
---

# IgScan BulkNGS Preprocessing Script

This script preprocesses raw bulk NGS paired-end sequencing data for the
IgScan workflow by merging reads, converting formats, filtering
sequences by length, collapsing duplicates, and sorting by abundance.

In this workflow, the following tasks are carried out:

1.  Merge of overlapping paired-end reads using the
    [*FLASH*](https://github.com/ebiggers/flash?tab=readme-ov-file) tool
    (version 1.2.11) [Magoc T, et al. *Bioinformatics* (2011)].

2.  Filter merged sequences by a user-defined minimum length (default:
    250 bp).

3.  Collapse identical sequences and counts their abundances.

4.  Output a sorted by abundance FASTA file with read counts in the
    headers.

## Requirements

-   Linux or macOS operating system

-   *FLASH* (included in the script directory under FLASH_1.2.11_Linux
    or FLASH_1.2.11_MacOS)

-   awk and sort utilities

## Usage

`bash ./BulkNGS_preprocessing/igscan_preprocess_bulkngs.sh <R1.fastq> <R2.fastq> <output.fasta> [min_length]`

`<R1.fastq>`: Forward reads FASTQ file

`<R2.fastq>`: Reverse reads FASTQ file

`<output.fasta>`: Output FASTA file name

`[min_length]`: (Optional) Minimum sequence length to keep (default:
250)

### Example

`bash ./BulkNGS_preprocessing/igscan_preprocess_bulkngs.sh sample_R1.fastq sample_R2.fastq  output.fasta 300`

### Output

A FASTA file (output.fasta) containing sequences longer or equal to the
minimum length (300 bp), sorted by abundance.

Sequence headers are formatted as **\>readX_n=COUNT**, where X is the
sequence index and COUNT is the number of reads collapsed.

### Note

Ensure FLASH binaries have execution permissions:

For linux: `chmod +x ./BulkNGS_preprocessing/FLASH_1.2.11_Linux/flash`

For macOS: `chmod +x ./BulkNGS_preprocessing/FLASH_1.2.11_MacOS/flash`

