#!/bin/bash

# Dependencies:
# brew install bowtie2 anaconda sratoolkit
# conda install -c bioconda parallel-fastq-dump

# Options:
N_THREADS=10

# Align SRR against Viral Genome
# $1: Reads, $2: Target Genome
function align() {
  mkdir -p "./out"
  cd "./out"
  if ! test -f "$1"; then
    echo "Downloading reads.."
    curl "https://sra-pub-run-odp.s3.amazonaws.com/sra/$1/$1" >"$1"
  fi

  echo "Analyzing reads.."
  parallel-fastq-dump --gzip --threads $N_THREADS --split-files --sra-id "$1"

  if ! test -f "${2}.fa"; then
    echo "Downloading target genome.."
    curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${2}" >"${2}.fa"
  fi

  echo "Indexing target genome.."
  bowtie2-build --threads $N_THREADS "${2}.fa" "${2}"

  echo "Adaptor & Quality trimming.."
  trimmomatic PE \
    ${1}_{1,2}.fastq.gz ${1}_{1,2}_{,un}paired.fastq.gz AVGQUAL:20 HEADCROP:3 \
    LEADING:3 TRAILING:3 MINLEN:75 -threads ${N_THREADS}

  # echo "Aligning reads to target genome.."
  target="${1}_${2}.bam"

  if test -f "${1}_2.fastq.gz"; then
    bowtie2 -p $N_THREADS -x "${2}" -1 "${1}_1_paired.fastq.gz" -2 "${1}_2_paired.fastq.gz" --no-unal | samtools sort -@$((N_THREADS - 1)) - >$target
  else
    bowtie2 -p $N_THREADS -x "${2}" -U "${1}_1.fastq.gz" --no-unal | samtools sort -@$((N_THREADS - 1)) - >$target
  fi

  samtools coverage $target
  samtools index $target
  echo "Finished!"
}

align "$1" "$2"
