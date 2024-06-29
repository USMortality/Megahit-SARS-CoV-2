#!/bin/bash

# # Download the reference genome
# wget https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3 -O MN908947.3.fa

# Extract the first 150bp of the reference genome
head -n 2 MN908947.3.fa | awk 'NR==1{print $0; next} {print substr($0, 1, 150)}' >MN908947.3_150bp.fasta

# Extract read IDs from the SAM file
awk '{if ($0 !~ /^@/) print $1}' filtered_reads/start_aligned_reads_trim_5_mismatches_3.sam | sort | uniq >read_ids.txt

# Extract the original reads from the FASTQ files using seqtk
seqtk subseq SRR10971381_1.fastq read_ids.txt >extracted_reads_1.fastq
seqtk subseq SRR10971381_2.fastq read_ids.txt >extracted_reads_2.fastq

# Convert the extracted reads to FASTA format
seqtk seq -a extracted_reads_1.fastq >extracted_reads_1.fasta
seqtk seq -a extracted_reads_2.fastq >extracted_reads_2.fasta

# Combine the reference genome and reads into a single FASTA file
cat MN908947.3_150bp.fasta extracted_reads_1.fasta extracted_reads_2.fasta >combined_150bp.fasta

# Align the sequences with MAFFT
mafft --auto combined_150bp.fasta >aligned_150bp_sequences.fasta

# Reformat the aligned FASTA file to have each sequence on a single line
awk '/^>/ {if (seq) print seq; print $0; seq=""; next} {seq=seq""$0} END {if (seq) print seq}' aligned_150bp_sequences.fasta >aligned_150bp_sequences_singleline.fasta

# Display the reformatted alignment
cat aligned_150bp_sequences_singleline.fasta
