#!/bin/sh

# brew install seqtk samtools temurin

cd out

# Create inverse of found sequence
seqtk seq -r k141_13590.fasta > k141_13590r.fasta

# Create index
samtools faidx MN908947.3.fasta
samtools faidx k141_13590r.fasta

# Create db
makeblastdb -in MN908947.3.fasta -dbtype nucl

# Create diff
blastn -query k141_13590r.fasta -db MN908947.3.fasta -evalue 1 -task megablast -outfmt 6 > k141_13590r_MN908947.3.crunch

# Use ACT to display diff
# http://sanger-pathogens.github.io/Artemis/ACT/
