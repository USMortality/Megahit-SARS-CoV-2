#!/bin/bash

# Fetch the sequences, extract the first 100 nt, and concatenate in one go
for id in 1 2 3; do
  curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=MN908947.$id" >MN908947.$id.fa
done

cat MN908947.*.fa >MN908947.fa

# # Align the sequences using mafft
# mafft --globalpair --maxiterate 1000 sequences_100nt.fa >temp.aln

# # Clean up
# rm seq_*_100nt.fasta sequences_100nt.fasta

# # Print the alignment
# cat temp.aln
