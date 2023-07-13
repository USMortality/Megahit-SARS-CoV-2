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

echo "----------"
esrun src/genome.ts -s SRR00001803 -g reference_genomes/NC_001803.fa
align_assemble SRR00001803 NC_001803

echo "----------"
esrun src/genome.ts -s SRR00002205 -g reference_genomes/NC_002205.fa
align_assemble SRR00002205 NC_002205

echo "----------"
esrun src/genome.ts -s SRR00002645 -g reference_genomes/NC_002645.fa
align_assemble SRR00002645 NC_002645

echo "----------"
esrun src/genome.ts -s SRR00005831 -g reference_genomes/NC_005831.fa
align_assemble SRR00005831 NC_005831

echo "----------"
esrun src/genome.ts -s SRR00006213 -g reference_genomes/NC_006213.fa
align_assemble SRR00006213 NC_006213

echo "----------"
esrun src/genome.ts -s SRR00006306 -g reference_genomes/NC_006306.fa
align_assemble SRR00006306 NC_006306

echo "----------"
esrun src/genome.ts -s SRR00006577 -g reference_genomes/NC_006577.fa
align_assemble SRR00006577 NC_006577

echo "----------"
esrun src/genome.ts -s SRR00026431 -g reference_genomes/NC_026431.fa
align_assemble SRR00026431 NC_026431

echo "----------"
esrun src/genome.ts -s SRR00036615 -g reference_genomes/NC_036615.fa
align_assemble SRR00036615 NC_036615

echo "----------"
esrun src/genome.ts -s SRR00038311 -g reference_genomes/NC_038311.fa
align_assemble SRR00038311 NC_038311

# SARS
echo "----------"
esrun src/genome.ts -s SRR00004718 -g reference_genomes/NC_004718.fa
align_assemble SRR00004718 NC_004718

# SARS2
echo "----------"
esrun src/genome.ts -s SRR00045512 -g reference_genomes/NC_045512.fa
align_assemble SRR00045512 NC_045512
