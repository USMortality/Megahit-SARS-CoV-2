#!/bin/bash

GENOME="MN908947.3" # Wuhan-Hu-1 / SARS-CoV-2
SRR="SRR10971381"   # Wu et al. 2020
MAX_MISMATCHES=3    # Set to 3 for allowing up to 3 mismatches
MAX_TRIM=10

if ! test -f "${GENOME}.fa"; then
  echo "Downloading target genome..."
  curl -o "${GENOME}.fa" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${GENOME}"
fi

# Download SARS-CoV-2 contigs
# Convert SRA format data to fastq format
if ! [ -f "${SRR}_1.fastq" ]; then
  fasterq-dump --progress ${SRR}
fi

# Get the maximum number of threads available on the machine minus one
threads=$(($(nproc) - 1))

# Create output directory
output_dir="filtered_reads"
mkdir -p $output_dir

# Function to process trimming and alignment
process_reads() {
  local trim_length=$1
  local ref_genome=$2

  echo "Trimming reads by $trim_length nucleotides from both ends"
  seqtk trimfq -b $trim_length -e $trim_length ${SRR}_1.fastq >trimmed_${SRR}_1.fastq
  seqtk trimfq -b $trim_length -e $trim_length ${SRR}_2.fastq >trimmed_${SRR}_2.fastq

  # Verify the trimmed files
  if [ ! -s trimmed_${SRR}_1.fastq ] || [ ! -s trimmed_${SRR}_2.fastq ]; then
    echo "Error: Trimmed fastq files were not created successfully."
    exit 1
  fi

  echo "Aligning the trimmed reads to the reference genome"
  bwa mem -t $threads $ref_genome trimmed_${SRR}_1.fastq trimmed_${SRR}_2.fastq >aligned_reads.sam

  # Verify the aligned reads file
  if [ ! -s aligned_reads.sam ]; then
    echo "Error: aligned_reads.sam was not created successfully."
    exit 1
  fi

  echo "Converting and sorting the SAM file to BAM"
  samtools view -@ $threads -bS aligned_reads.sam | samtools sort -@ $threads -o aligned_reads_sorted.bam -

  # Verify the sorted BAM file
  if [ ! -s aligned_reads_sorted.bam ]; then
    echo "Error: aligned_reads_sorted.bam was not created successfully."
    exit 1
  fi

  echo "Filtering reads that align to the start of the reference genome with up to $MAX_MISMATCHES mismatches"
  output_file="$output_dir/start_aligned_reads_trim_${trim_length}_mismatches_${MAX_MISMATCHES}.sam"
  if [ "$MAX_MISMATCHES" -eq 0 ]; then
    samtools view -@ $threads aligned_reads_sorted.bam | awk '$4 == 1 && $6 ~ /^[0-9]+M$/' >$output_file
  else
    samtools view -@ $threads aligned_reads_sorted.bam | awk -v max_mismatches=$MAX_MISMATCHES '$4 == 1 && $6 ~ /^[0-9]+M$/ && $0 ~ "NM:i:" && substr($0, index($0, "NM:i:") + 5, 1) <= max_mismatches' >$output_file
  fi

  # Count the reads that align to the start
  start_count=0
  if [ -s $output_file ]; then
    start_count=$(wc -l <$output_file)
  fi

  # Store the results in the array
  start_results[$trim_length]=$start_count

  echo "Results for trimming $trim_length nucleotides: Start alignments with up to $MAX_MISMATCHES mismatches: $start_count"
}

# Main script execution
ref_genome="${GENOME}.fa"

# Step 1: Index the reference genome
echo "Step 1: Indexing the reference genome"
bwa index $ref_genome

# Initialize array to store the results
declare -a start_results

# Step 2: Loop over trimming lengths
for trim_length in $(seq 0 $MAX_TRIM); do
  process_reads $trim_length $ref_genome
done

# Print the overall results
echo "Overall results:"
for trim_length in $(seq 0 $MAX_TRIM); do
  echo "Trim length $trim_length: Start alignments with up to $MAX_MISMATCHES mismatches: ${start_results[$trim_length]}"
done

echo "Filtered reads are saved in the '$output_dir' directory."
