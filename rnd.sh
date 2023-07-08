#!/bin/bash

esrun genome.ts -g out/MN908947.3.fa
# esrun genome.ts
seqkit seq --min-len 50 --max-len 150 out/SRR00000000_1.fastq >out/SRR00000000_1_f.fastq
seqkit seq --min-len 50 --max-len 150 out/SRR00000000_2.fastq >out/SRR00000000_2_f.fastq
cd out
tar -czf SRR00000000_1_f.fastq.gz SRR00000000_1_f.fastq
tar -czf SRR00000000_2_f.fastq.gz SRR00000000_2_f.fastq
mv SRR00000000_1_f.fastq.gz SRR00000000_1.fastq.gz
mv SRR00000000_2_f.fastq.gz SRR00000000_2.fastq.gz
cd ~-
./align.sh -s SRR00000000 -g MN908947.3
megahit -1 out/SRR00000000_1.fastq.gz -2 out/SRR00000000_2.fastq.gz -o out/out | tee output.txt
