#!/bin/bash

# Options:
N_THREADS=$(nproc --all)

# Align SRR against Viral Genome
# $1: Reads, $2: Target Genome
function align() {
  mkdir -p "./out"
  cd "./out"
  if ! test -f "$1"; then
    echo "Downloading reads.."
    curl "https://sra-pub-run-odp.s3.amazonaws.com/sra/$1/$1" >"$1"
  fi

  if ! test -f "${1}_1.fastq.gz"; then
    echo "Analyzing reads.."
    parallel-fastq-dump --gzip --threads $N_THREADS --split-files --sra-id "$1"
  fi

  if ! test -f "${2}.fa"; then
    echo "Downloading target genome.."
    curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${2}" >"${2}.fa"
  fi

  if ! test -f "${2}.rev.1.bt2"; then
    echo "Indexing target genome.."
    bowtie2-build --threads $N_THREADS "${2}.fa" "${2}"
  fi

  if ! test -f "${2}.fa.amb"; then
    echo "Indexing target genome.."
    bwa index ${2}.fa
  fi

  echo "Aligning reads to target genome.."
  target="${1}_${2}"

  if test -f "${1}_2.fastq.gz"; then
    # bowtie2 -p $N_THREADS -x "${2}" -1 "${1}_1.fastq.gz" -2 "${1}_2.fastq.gz" --no-unal | samtools sort -@$((N_THREADS - 1)) - >$target.bam
    bwa mem -t $N_THREADS ${2}.fa ${1}_1.fastq.gz ${1}_2.fastq.gz | samtools view -Sb - >$target.bam
  else
    # bowtie2 -p $N_THREADS -x "${2}" -U "${1}_1.fastq.gz" --no-unal | samtools sort -@$((N_THREADS - 1)) - >$target.bam
    bwa mem -t $N_THREADS ${2}.fa ${1}_1.fastq.gz | samtools view -Sb - >$target.bam
  fi

  # required for bwa mem
  samtools sort $target.bam -o $target.sort.bam
  samtools view -b -F 4 $target.sort.bam >$target.sort.mapped.bam
  samtools index $target.sort.mapped.bam

  COV_PCT=$(samtools coverage ${1}_${2}.sort.mapped.bam -H | cut -f6)
  echo "Coverage: ${COV_PCT}%"
  # TODO: Replace this with different endpoint.
  curl -X POST -H "Content-Type: text/plain" --data "${1};${COV_PCT}" https://csvfilter.mortality.watch/error

  echo "Cleaning up.."
  cd ~-
  # rm -rf out
  echo "Finished!"
}

align "$1" "$2"
