#!/bin/bash

# Settings:
GENOME="MN908947.3" # Wuhan-Hu-1 / SARS-CoV-2
SRR="SRR10971381"   # Wu et al. 2020
THREADS=$(nproc)

sam_to_sorted_bam() {
  # Convert SAM to BAM
  samtools view -S -b ${1}.sam >${1}.bam
  # Sort BAM file
  samtools sort ${1}.bam -o ${1}.sorted.bam
  # Index the sorted BAM file
  samtools index ${1}.sorted.bam
}

# Download Target Genome
if ! test -f "${GENOME}.fa"; then
  curl -o "${GENOME}.fa" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${GENOME}"
fi

# Download SARS-CoV-2 contigs. Convert SRA format data to fastq format
if ! [ -f "${SRR}_1.fastq" ]; then
  fasterq-dump -e $THREADS --progress ${SRR}
fi

# Align the reads to the genome
minimap2 -a -Y --sam-hit-only ${GENOME}.fa ${SRR}_1.fastq ${SRR}_2.fastq | samtools sort -@7 - >${SRR}.sam
sam_to_sorted_bam ${SRR}

# Filter reads with 0-9 clippings; less than expected read errors
Q=31
samtools view -h ${SRR}.sam |
  gawk '/^@/ || ($6 ~ /^([0-9]S)?[0-9]+M([0-9]S)?$/ && $6 !~ /I/ &&
      $6 !~ /D/) {print}' |
  gawk 'NR==FNR{a[NR]=$0;next}{match($0,/NM:i:([0-9]+)/,m)}
      m[1]<=a[length($10)]' \
    <(Rscript --no-init-file -e \
      "cat(qbinom(.95, 1:1e4, 10^(${Q}/-10)), sep='\\n')") - \
    >${SRR}.hq.sam

sam_to_sorted_bam ${SRR}.hq

# Convert filtered BAM back to FASTQ
bedtools bamtofastq -i ${SRR}.hq.sam -fq ${SRR}.hq.fastq

# Filter reads with no clippings
samtools view -h ${SRR}.hq.sam |
  gawk '/^@/ || $6 ~ /^[0-9]+M$/ {print}' - >${SRR}.hq.0.sam
bedtools bamtofastq -i ${SRR}.hq.0.sam -fq ${SRR}.hq.0.fastq
sam_to_sorted_bam ${SRR}.hq.0

# Filter reads with 3 clippings
samtools view -h ${SRR}.hq.sam |
  gawk '/^@/ || $6 ~ /^3S[0-9]+M$/ || $6 ~ /^[0-9]+M3S$/ {print}' - >${SRR}.hq.3.sam
bedtools bamtofastq -i ${SRR}.hq.3.sam -fq ${SRR}.hq.3.fastq
sam_to_sorted_bam ${SRR}.hq.3

# Filter reads with 4 clippings
samtools view -h ${SRR}.hq.sam |
  gawk '/^@/ || $6 ~ /^4S[0-9]+M$/ || $6 ~ /^[0-9]+M4S$/ {print}' - >${SRR}.hq.4.sam
bedtools bamtofastq -i ${SRR}.hq.4.sam -fq ${SRR}.hq.4.fastq
sam_to_sorted_bam ${SRR}.hq.4

# Filter reads with 5 clippings
samtools view -h ${SRR}.hq.sam |
  gawk '/^@/ || $6 ~ /^5S[0-9]+M$/ || $6 ~ /^[0-9]+M5S$/ {print}' - >${SRR}.hq.5.sam
bedtools bamtofastq -i ${SRR}.hq.5.sam -fq ${SRR}.hq.5.fastq
sam_to_sorted_bam ${SRR}.hq.5
