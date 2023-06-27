#!/bin/bash

# Dependencies:
# brew install bowtie2 anaconda sratoolkit
# conda install -c bioconda parallel-fastq-dump trimmomatic fastp

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

  echo "Quality trimming reads.."
  # trimmomatic PE \
  #   $1_{1,2}.fastq.gz $1_{1,2}.{,un}paired.fq.gz \
  #   ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 \
  #   -threads $N_THREADS
  fastp --in1 ${1}_1.fastq.gz --in2 ${1}_2.fastq.gz \
    --out1 ${1}_1_qt.fastq.gz --out2 ${1}_2_qt.fastq.gz \
    --thread $N_THREADS --detect_adapter_for_pe 1

  if ! test -f "${2}.fa"; then
    echo "Downloading target genome.."
    curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${2}" >"${2}.fa"
    bowtie2-build --threads $N_THREADS "${2}.fa" "${2}"
  fi

  echo "Aligning reads to target genome.."
  target="${1}_${2}.bam"
  if test -f "${1}_2.fastq.gz"; then
    bowtie2 -p $N_THREADS -x "${2}" -1 "${1}_1_qt.fastq.gz" -2 "${1}_2_qt.fastq.gz" --no-unal | samtools sort -@$((N_THREADS - 1)) - >$target
  else
    bowtie2 -p $N_THREADS -x "${2}" -U "${1}_qt_1.fastq.gz" --no-unal | samtools sort -@$((N_THREADS - 1)) - >$target
  fi

  samtools coverage $target
  samtools index $target
  echo "Finished!"
}

align "$1" "$2"
