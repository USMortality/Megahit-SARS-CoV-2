#!/bin/bash

# Create a file with the read ID
echo "SRR10971381.383615" >read_id.txt

# Extract the specific read from the FASTQ files using seqtk
seqtk subseq SRR10971381_1.fastq read_id.txt >extracted_read_1.fastq
seqtk subseq SRR10971381_2.fastq read_id.txt >extracted_read_2.fastq

# Reverse complement the second read using seqkit
seqkit seq --reverse --complement -t DNA extracted_read_1.fastq >extracted_read_1_rev_comp.fastq

# Display the extracted reads together for comparison
echo "Read from SRR10971381_1.fastq:"
cat extracted_read_1.fastq
echo
echo "Reverse complement of the read from SRR10971381_2.fastq:"
cat extracted_read_1_rev_comp.fastq
