#!/bin/bash

rm -rf out/*

for f in ./reference_genomes/*.fa; do
  name=$(basename $f .fa)
  srr="SRR_$name"
  esrun src/genome.ts -s $srr -g $f -n 1000 -e 0.01
  # wgsim -N 100000 $f out/${srr}_1.fastq out/${srr}_2.fastq
  echo "----------------------------------------------------------------------"
done

cd out
cat *_1.fastq >SRR00000001_1.fastq
cat *_2.fastq >SRR00000001_2.fastq

tar -czf SRR00000001_1.fastq.gz SRR00000001_1.fastq
tar -czf SRR00000001_2.fastq.gz SRR00000001_2.fastq
cd ~-

# De-Novo Assemble
rm -rf out/megahit
megahit -1 out/SRR00000001_1.fastq.gz -2 out/SRR00000001_2.fastq.gz -o out/megahit | tee output.txt

cd out/megahit
makeblastdb -in final.contigs.fa -dbtype nucl

echo "Blasting reference genomes against contigs..."
for f in ../../reference_genomes/*.fa; do
  name=$(basename $f .fa)
  blastn -query $f -db final.contigs.fa -evalue 1 -task megablast >$name.crunch
done

echo "Done"
