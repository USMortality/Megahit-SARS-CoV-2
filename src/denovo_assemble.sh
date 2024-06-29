#!/bin/bash

# Settings:
GENOME="MN908947.3" # Wuhan-Hu-1 / SARS-CoV-2
SRR="SRR10971381"   # Wu et al. 2020
THREADS=$(nproc)

# Trim Reads
java -jar trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
  SRR10971381.hq.fastq SRR10971381.hq.trimmed.fastq.gz \
  AVGQUAL:20 HEADCROP:12 LEADING:3 TRAILING:3 MINLEN:75 -threads ${THREADS}

# Run megahit assembly (trimmed)
rm -rf out
./megahit/bin/megahit -r SRR10971381.hq.trimmed.fastq.gz -o out 2>&1 | tee output.txt

# Analyze produced contigs
cd out
makeblastdb -in final.contigs.fa -dbtype nucl

echo "Blasting reference genomes against contigs..."
blastn -query ../${GENOME}.fa -db final.contigs.fa -evalue 1 -task megablast >${SRR}.crunch
cd ..

# Zip log & output.
zip -r9 out.zip out/final.contigs.fa output.txt out/${SRR}.crunch

echo "Done"
