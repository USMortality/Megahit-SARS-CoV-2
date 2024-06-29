#!/bin/bash

# Dependencies:
# brew install minimap2 samtools pipx sratoolkit

# Settings:
GENOME="MN908947.3" # Wuhan-Hu-1 / SARS-CoV-2
SRR="SRR10971381"   # Wu et al. 2020
THREADS=$(($(nproc) - 1))

# Download Target Genome
if ! test -f "${GENOME}.fa"; then
  curl -o "${GENOME}.fa" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${GENOME}"
fi

# Download SARS-CoV-2 contigs. Convert SRA format data to fastq format
if ! [ -f "${SRR}_1.fastq" ]; then
  fasterq-dump -e $THREADS --progress ${SRR}
fi

# Align Reads
minimap2 -a -Y --sam-hit-only -t ${THREADS} ${GENOME}.fa ${SRR}_1.fastq ${SRR}_2.fastq | samtools sort -@ ${THREADS} - >temp.bam
# Filter Aligned Reads
samtools view temp.bam | awk '$4==1' | awk '$6 ~ /^([0-5]S)?[0-9]+M$/ {print}' | awk 'BEGIN {OFS="\t"} {split($12, a, ":"); if (a[3] <= 3) print $1, $6, $10}' | cut -d. -f2 | sed 's/\t/:/' | seqkit tab2fx | cat $GENOME.fa - | seqkit subseq -r1:150 | mafft --quiet - | alv -
