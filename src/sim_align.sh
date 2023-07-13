#!/bin/bash

# Options:
N_THREADS=($(nproc --all) -1)

# Align SRR against Viral Genome
# $1: Reads, $2: Target Genome
function align() {
  cd "./out"

  if ! test -f "${2}.fa"; then
    echo "Downloading target genome.."
    curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${2}" >"${2}.fa"
  fi

  echo "Generating reads..."
  wgsim -N 100000 "${2}.fa" "${1}_1.fq" "${1}_2.fq"

  echo "Indexing target genome.."
  bowtie2-build --threads $N_THREADS "${2}.fa" "${2}"

  echo "Aligning reads to target genome.."
  target="${1}_${2}.bam"

  bowtie2 -p $N_THREADS -x "${2}" -1 "${1}_1.fq" -2 "${1}_2.fq" --no-unal | samtools sort -@$((N_THREADS - 1)) - >$target

  samtools coverage $target
  samtools index $target
  echo "Finished!"
}

align "$1" "$2"
