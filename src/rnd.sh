#!/bin/bash

rm -rf out/*

esrun ./src/genome.ts -g ./reference_genomes/NC_045512.fa -s SRR00000001

cd out
# wgsim -N 100000 -e 0.01 -1 150 -2 150 -d 500 ../reference_genomes/NC_045512.fa SRR00000001_1.fastq SRR00000001_2.fastq

tar -czf SRR00000001_1.fastq.gz SRR00000001_1.fastq
tar -czf SRR00000001_2.fastq.gz SRR00000001_2.fastq
cd ~-

cp reference_genomes/NC_045512.fa out/.
# Align reads
./src/align.sh -s SRR00000001 -g NC_045512 -a bwa
samtools depth -a out/SRR00000001_NC_045512_bwa.bam | head -n1
samtools depth -a out/SRR00000001_NC_045512_bwa.bam | tail -n1

# # De-Novo Assemble
# rm -rf out/megahit
# megahit -1 out/SRR00000001_1.fastq.gz -2 out/SRR00000001_2.fastq.gz -o out/megahit | tee output.txt
