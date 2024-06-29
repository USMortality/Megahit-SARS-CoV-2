#!/bin/bash

# Visualization Functions
cohl() (
  awk 'NR==1 {split($0, a, ""); next}
   	{split($0, b, ""); for(i=1; i<=length(b); i++)
    printf "%s", (b[i] == a[i] ? b[i] : "\033[31m" b[i] "\033[0m");
    print ""}' "$@"
)

alir() (
  awk '{a[NR]=$0; if(length($0) > l) l=length($0)}
   	END {for(i=1; i<=NR; i++) printf "%" l "s\n", a[i]}' "$@"
)

cohl0() (
  x=$(cat)
  paste -d' ' <(seqkit seq -n <<<"$x" | alir) <(seqkit seq -s <<<"$x" |
    sed 1p | cohl)
)

# Function to process data for a given study
process_study() {
  local bam_file=$1
  local fastq_files=$2
  local Q=$3

  # Align:
  minimap2 -a -Y --sam-hit-only MN908947.3.fa $fastq_files |
    samtools sort -@7 - >$bam_file

  # Filter HQ Reads:
  samtools view -h $bam_file |
    gawk '/^@/ || ($6 ~ /^([0-9]S)?[0-9]+M([0-9]S)?$/ && $6 !~ /I/ &&
      $6 !~ /D/) {print}' |
    gawk 'NR==FNR{a[NR]=$0;next}{match($0,/NM:i:([0-9]+)/,m)}
      m[1]<=a[length($10)]' \
      <(Rscript --no-init-file -e \
        "cat(qbinom(.95, 1:1e4, 10^(${Q}/-10)), sep='\\n')") - \
      >${bam_file%.sam}.hq.sam

  # Visualize Reads at the start
  echo "All Reads:"
  visualize_reads_start $bam_file
  echo "HQ Reads:"
  visualize_reads_start ${bam_file%.sam}.hq.sam "HQ"

  # Visualize Reads at the end
  echo "All Reads:"
  visualize_reads_end $bam_file
  echo "HQ Reads:"
  visualize_reads_end ${bam_file%.sam}.hq.sam "HQ"
}

# Function to visualize reads at the start
visualize_reads_start() {
  local sam_file=$1
  local prefix=$2
  local output_prefix=${prefix:+${prefix}_}

  local start
  start=$(samtools view ${sam_file} | head -n1 | awk '{print $4}')
  echo "First match at position: ${start}"
  # All
  samtools view $sam_file | awk -v start="$start" '$4 == start' |
    cut -f1,6,10 | cut -d. -f2 | sed $'s/\t/:/' | seqkit tab2fx |
    cat <(cut -d, -f1 MN908947.3.fa | sed 's/ .*//') - |
    seqkit subseq -r1:60 | mafft --quiet - | cohl0
}

# Function to visualize reads at the end
visualize_reads_end() {
  local sam_file=$1
  local prefix=$2
  local output_prefix=${prefix:+${prefix}_}

  # All
  bedtools bamtobed -i $sam_file | sort -k3,3nr | head -n30 |
    awk 'NR<=10 {print $4}' >top_ids.txt
  samtools view $sam_file | grep -F -w -f top_ids.txt | cut -f1,6,10 |
    sed 's/\t/:/' | seqkit tab2fx 2>/dev/null |
    cat <(seqkit subseq -r -200:-1 MN908947.3.fa 2>/dev/null |
      sed 's/ .*//') - | mafft --quiet - |
    seqkit subseq -r -180:-1 2>/dev/null | cohl0
}

# echo "### Wu et al., 2020"
# process_study "wu.bam" "SRR10971381_1.fastq SRR10971381_2.fastq" 31

echo "### Kim et al., 2020"
process_study "kim.bam" "VeroInf24h.all.fastq" 21

echo "### Moreno et al., 2022"
process_study "moreno.bam" "SRR11140745.fastq" 21

# echo "### PacBio, 2021"
# process_study "pacbio.bam" "hifi_reads.fastq" 33
