#!/bin/bash

# Dependencies:
# brew install bowtie2 anaconda
# conda install -c bioconda parallel-fastq-dump

# Options:
N_THREADS=16

# Align SRR against Viral Genome
# $1: Reads, $2: Target Genome
function align() {
  mkdir -p "./out"
  cd "./out"
  if ! test -f "$1"; then
    echo "Downloading reads.."
    curl "https://sra-pub-run-odp.s3.amazonaws.com/sra/$1/$1" > "$1"
    parallel-fastq-dump --gzip --threads $N_THREADS --split-files --sra-id "$1"
  fi
  if ! test -f "${2}.fa"; then
    echo "Downloading target genome.."
    curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${2}" >"${2}.fa"
    bowtie2-build --threads $N_THREADS "${2}.fa" "${2}"
  fi
  echo "Aligning reads to target genome.."
  bowtie2 -p $N_THREADS -x "${2}" -1 "${1}_1.fastq.gz" -2 "${1}_2.fastq.gz" --no-unal | samtools sort -@2 - >temp.bam
  samtools coverage temp.bam
  echo "Finished!"
}

align "$1" "$2"
