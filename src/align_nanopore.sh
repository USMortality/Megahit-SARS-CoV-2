#!/bin/bash

# Nanopore reads:
# DL file from here: https://osf.io/dfvbz

mkdir -p out && cd out/

# SARS-CoV-2 ref genome:
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=NC_045512" >"NC_045512.fa"
bwa index -p NC_045512 NC_045512.fa

minimap2 -a --sam-hit-only -x map-ont NC_045512.fa <(seqkit seq --rna2dna ~/Downloads/VeroInf24h.all.fastq.zst) >alignment.sam
samtools sort alignment.sam >alignment.bam
samtools index alignment.bam
