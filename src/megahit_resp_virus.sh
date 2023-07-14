#!/bin/bash

rm -rf out/*
mkdir -p out/megahit

function align_assemble {
  cd out
  tar -czf ${1}_1.fastq.gz ${1}_1.fastq
  tar -czf ${1}_2.fastq.gz ${1}_2.fastq
  cd ~-

  #  Align reads
  ./src/align.sh -s ${1} -g ${2} -a bwa &>./out/align.log
  echo "BWA Alignment: $(cat out/align.log | tail -n2 | head -n1)"
  ./src/align.sh -s ${1} -g ${2} -a bowtie2 &>./out/align.log
  echo "Bowtie2 Alignment: $(cat out/align.log | tail -n2 | head -n1)"
  ./src/align.sh -s ${1} -g ${2} -a minimap2 &>./out/align.log
  echo "Minimap2 Alignment: $(cat out/align.log | tail -n2 | head -n1)"

  # De-Novo Assemble
  cd out
  rm -rf megahit/${1}
  megahit -1 ${1}_1.fastq.gz -2 ${1}_2.fastq.gz -o megahit/${1}/ &>/dev/null
  cat megahit/${1}/log | tail -n2 | head -n1

  cd ~-
  echo "Reference Genome Length: $(cat reference_genomes/${2}.fa | samtools view | awk '{print length($10)}')"
}

for f in ./reference_genomes/*.fa; do
  name=$(basename $f .fa)
  srr="SRR_$name"
  # esrun src/genome.ts -s $srr -g $f -n 1000 -e 0.01
  wgsim -N 100000 $f out/${srr}_1.fastq out/${srr}_2.fastq
  align_assemble $srr $name
  echo "----------------------------------------------------------------------"
done
