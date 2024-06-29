#!/bin/bash

# Settings:
GENOME="MN908947.3" # Wuhan-Hu-1 / SARS-CoV-2
SRR="SRR10971381"   # Wu et al. 2020
THREADS=$(nproc)

# Download Target Genome
if ! test -f "${GENOME}.fa"; then
  curl -o "${GENOME}.fa" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${GENOME}"
fi

# Download SARS-CoV-2 contigs. Convert SRA format data to fastq format
if ! [ -f "${SRR}_1.fastq" ]; then
  fasterq-dump -e $THREADS --progress ${SRR}
fi

# Get Latest Complete Human Genome
if ! test -f "ncbi_dataset.zip"; then
  wget --content-disposition 'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_009914755.1/download?include_annotation_type=GENOME_FASTA'
  unzip ncbi_dataset.zip
fi

# # Align Wu et al reads to the latest human genome
# minimap2 -a --sam-hit-only ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna SRR10971381_1.fastq SRR10971381_2.fastq >SRR10971381_human.sam

# Define input FASTQ files and final output file
FASTQ1="SRR10971381_1.fastq"
FASTQ2="SRR10971381_2.fastq"
FINAL_OUTPUT="SRR10971381_non_human.fastq"

# Extract the read IDs from SAM file and save to human_read_ids.txt
# samtools view SRR10971381_human.sam | awk '{print $1}' | sort | uniq >human_read_ids.txt

# Filter reads from the first FASTQ file
awk 'BEGIN {while ((getline < "human_read_ids.txt") > 0) l["@" $1] = 1} /^@/ {f = !l[$1]} f' $FASTQ1 >filtered_1.fastq

# Filter reads from the second FASTQ file
awk 'BEGIN {while ((getline < "human_read_ids.txt") > 0) l["@" $1] = 1} /^@/ {f = !l[$1]} f' $FASTQ2 >filtered_2.fastq

cat filtered_1.fastq filtered_2.fastq >$FINAL_OUTPUT

# Clean up intermediate files
rm filtered_1.fastq filtered_2.fastq human_read_ids.txt

echo "Filtering complete. Non-human reads FASTQ file is $FINAL_OUTPUT"
