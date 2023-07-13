#!/bin/bash

# rm -rf out/*

# Generate Reads based on existing 10 resp. viruses:
esrun genome.ts -s SRR00001803 -g reference_genomes/NC_001803.fa
esrun genome.ts -s SRR00002205 -g reference_genomes/NC_002205.fa
esrun genome.ts -s SRR00002645 -g reference_genomes/NC_002645.fa
esrun genome.ts -s SRR00005831 -g reference_genomes/NC_005831.fa
esrun genome.ts -s SRR00006213 -g reference_genomes/NC_006213.fa
esrun genome.ts -s SRR00006306 -g reference_genomes/NC_006306.fa
esrun genome.ts -s SRR00006577 -g reference_genomes/NC_006577.fa
esrun genome.ts -s SRR00026431 -g reference_genomes/NC_026431.fa
esrun genome.ts -s SRR00036615 -g reference_genomes/NC_036615.fa
esrun genome.ts -s SRR00038311 -g reference_genomes/NC_038311.fa

# esrun random_reads.ts
# esrun genome.ts -g out/MN908947.3.fa

cd out

# cat *_1.fastq > SRR00000001_1.fastq
# cat *_2.fastq > SRR00000001_2.fastq

# wgsim -N 5000 MN908947.3.fa SRR00000001_1.fastq SRR00000001_2.fastq

# seqkit seq --min-len 50 --max-len 150 SRR00000001_1.fastq >SRR00000001_1_f.fastq
# seqkit seq --min-len 50 --max-len 150 SRR00000001_2.fastq >SRR00000001_2_f.fastq
tar -czf SRR00000001_1.fastq.gz SRR00000001_1.fastq
tar -czf SRR00000001_2.fastq.gz SRR00000001_2.fastq
cd ~-

#  Align reads
# bash -x ./align.sh -s SRR00000001 -g MN908947.3

# De-Novo Assemble
rm -rf out/megahit
megahit -1 out/SRR00000001_1.fastq.gz -2 out/SRR00000001_2.fastq.gz -o out/megahit | tee output.txt
