#!/bin/bash

rm -rf out/*

# Generate Reads based on existing 10 resp. viruses:
esrun src/genome.ts -s SRR00001803 -g reference_genomes/NC_001803.fa
esrun src/genome.ts -s SRR00002205 -g reference_genomes/NC_002205.fa
esrun src/genome.ts -s SRR00002645 -g reference_genomes/NC_002645.fa
esrun src/genome.ts -s SRR00005831 -g reference_genomes/NC_005831.fa
esrun src/genome.ts -s SRR00006213 -g reference_genomes/NC_006213.fa
esrun src/genome.ts -s SRR00006306 -g reference_genomes/NC_006306.fa
esrun src/genome.ts -s SRR00006577 -g reference_genomes/NC_006577.fa
esrun src/genome.ts -s SRR00026431 -g reference_genomes/NC_026431.fa
esrun src/genome.ts -s SRR00036615 -g reference_genomes/NC_036615.fa
# esrun src/genome.ts -s SRR00038311 -g reference_genomes/NC_038311.fa

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

blastn -query ../../reference_genomes/NC_001803.fa -db final.contigs.fa -evalue 1 -task megablast >NC_001803.crunch
blastn -query ../../reference_genomes/NC_002205.fa -db final.contigs.fa -evalue 1 -task megablast >NC_002205.crunch
blastn -query ../../reference_genomes/NC_002645.fa -db final.contigs.fa -evalue 1 -task megablast >NC_002645.crunch
blastn -query ../../reference_genomes/NC_005831.fa -db final.contigs.fa -evalue 1 -task megablast >NC_005831.crunch
blastn -query ../../reference_genomes/NC_006213.fa -db final.contigs.fa -evalue 1 -task megablast >NC_006213.crunch
blastn -query ../../reference_genomes/NC_006306.fa -db final.contigs.fa -evalue 1 -task megablast >NC_006306.crunch
blastn -query ../../reference_genomes/NC_006577.fa -db final.contigs.fa -evalue 1 -task megablast >NC_006577.crunch
blastn -query ../../reference_genomes/NC_026431.fa -db final.contigs.fa -evalue 1 -task megablast >NC_026431.crunch
blastn -query ../../reference_genomes/NC_036615.fa -db final.contigs.fa -evalue 1 -task megablast >NC_036615.crunch
