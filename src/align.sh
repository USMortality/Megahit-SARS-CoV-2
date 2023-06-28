#!/bin/bash

# Options:
N_THREADS=($(nproc --all) -1)

# Align SRR against Viral Genome
# $1: Reads, $2: Target Genome
function align() {
  mkdir -p "./out"
  cd "./out"

  echo "Indexing target genome.."
  bowtie2-build --threads $N_THREADS "${2}.fa" "${2}"

  echo "Aligning reads to target genome.."
  target="${1}_${2}.bam"

  bowtie2 -p $N_THREADS -x "${2}" -U "${1}.fastq" --no-unal | samtools sort -@$((N_THREADS - 1)) - >$target

  samtools coverage $target
  samtools index $target
  echo "Finished!"
}

align "$1" "$2"
