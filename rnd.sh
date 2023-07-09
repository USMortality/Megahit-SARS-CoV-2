#!/bin/bash

esrun genome.ts -g out/MN908947.3.fa

cd out

# wgsim -N 5000 MN908947.3.fa SRR00000000_1.fastq SRR00000000_2.fastq

# seqkit seq --min-len 50 --max-len 150 SRR00000000_1.fastq >SRR00000000_1_f.fastq
# seqkit seq --min-len 50 --max-len 150 SRR00000000_2.fastq >SRR00000000_2_f.fastq
tar -czf SRR00000000_1.fastq.gz SRR00000000_1.fastq
tar -czf SRR00000000_2.fastq.gz SRR00000000_2.fastq
cd ~-

#  Align reads
bash -x ./align.sh -s SRR00000000 -g MN908947.3

# De-Novo Assemble
rm -rf out/megahit
megahit -1 out/SRR00000000_1.fastq.gz -2 out/SRR00000000_2.fastq.gz -o out/megahit | tee output.txt
