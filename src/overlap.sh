#!/bin/bash

# Define the identifier and FASTQ files
IDENTIFIER="SRR10971381.206220"
FASTQ1="SRR10971381_1.fastq"
FASTQ2="SRR10971381_2.fastq"

# Extract and display the forward read
FORWARD_READ=$(seqtk subseq $FASTQ1 <(echo "$IDENTIFIER"))

# Extract and display the reverse read
REVERSE_READ=$(seqtk subseq $FASTQ2 <(echo "$IDENTIFIER"))

# Extract the sequence from the forward read
FORWARD_SEQUENCE=$(echo "$FORWARD_READ" | sed -n '2p')
REVERSE_SEQUENCE=$(echo "$REVERSE_READ" | sed -n '2p')

echo "Forward Sequence: $FORWARD_SEQUENCE"
echo "Reverse Sequence: $REVERSE_SEQUENCE"

# Reverse complement the forward sequence
REVERSE_COMPLEMENT=$(echo "$FORWARD_SEQUENCE" | tr 'ATCG' 'TAGC' | rev)
echo "Reverse Complement of Forward Sequence: $REVERSE_COMPLEMENT"

# Check for overlap
OVERLAP=$(echo "$REVERSE_SEQUENCE" | grep -o "$REVERSE_COMPLEMENT")
if [[ -n "$OVERLAP" ]]; then
  echo "Overlap found: $OVERLAP"
else
  echo "No overlap found."
fi
