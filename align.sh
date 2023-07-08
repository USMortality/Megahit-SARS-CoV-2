#!/bin/bash

# Dependencies:
# brew install bowtie2 anaconda sratoolkit
# conda install -c bioconda parallel-fastq-dump minimap2

# Options:
N_THREADS=$(nproc --all)

# Align SRR against Viral Genome with bwa mem
function align_bwa() {
  echo "Aligning with ${3}.."
  if ! test -f "${2}.fa.amb"; then
    echo "Indexing target genome.."
    bwa index ${2}.fa
  fi

  echo "Aligning reads to target genome.."
  if test -f "${1}_2.${4}"; then
    bwa mem -t $N_THREADS ${2}.fa ${1}_1.${4} ${1}_2.${4} | samtools view -b -F 4 | samtools sort -@3 - >$TARGET.bam
  else
    bwa mem -t $N_THREADS ${2}.fa ${1}_1.${4} | samtools view -b -F 4 | samtools sort -@3 - >$TARGET.bam
  fi
}

# Align SRR against Viral Genome with bowtie2
function align_bowtie2() {
  echo "Aligning with ${3}.."
  if ! test -f "${2}.fai"; then
    echo "Indexing target genome.."
    bowtie2-build --threads $N_THREADS "${2}.fa" "${2}"
  fi

  echo "Aligning reads to target genome.."
  if test -f "${1}_2.${4}"; then
    bowtie2 -p $N_THREADS -x "${2}" -1 "${1}_1.${4}" -2 "${1}_2.${4}" --no-unal | samtools sort -@$((N_THREADS - 1)) - >$TARGET.bam
  else
    bowtie2 -p $N_THREADS -x "${2}" -U "${1}_1.${4}" --no-unal | samtools sort -@$((N_THREADS - 1)) - >$TARGET.bam
  fi
}

# Align SRR against Viral Genome with minimap2
function align_minimap2() {
  echo "Aligning with ${3}.."

  echo "Aligning reads to target genome.."
  if test -f "${1}_2.${4}"; then
    minimap2 "${2}.fa" "${1}_1.${4}" "${1}_2.${4}" -a --sam-hit-only | samtools sort -@$((N_THREADS - 1)) - >$TARGET.bam
  else
    minimap2 "${2}.fa" "${1}_1.${4}" -a --sam-hit-only | samtools sort -@$((N_THREADS - 1)) - >$TARGET.bam
  fi
}

# Parse cmd line opts.
ALIGNER=bwa
PARAMS=""
while (("$#")); do
  case "$1" in
  -t | --trim)
    TRIM=1
    shift
    ;;
  -s | --sra)
    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
      SRA=$2
      shift 2
    else
      echo "Error: Argument for $SRA is missing" >&2
      exit 1
    fi
    ;;
  -g | --genome)
    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
      GENOME=$2
      shift 2
    else
      echo "Error: Argument for $SRA is missing" >&2
      exit 1
    fi
    ;;
  -a | --aligner)
    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
      ALIGNER=$2
      shift 2
    else
      echo "Error: Argument for $SRA is missing" >&2
      exit 1
    fi
    ;;
  -h | --help)
    echo "Usage: ./align.sh -s/--sra SRR10971381 -g/--genome MN908947.3 (-a/--aligner [bwa/bowtie2]) (-t/--trim)"
    exit 0
    ;;
  -* | --*=) # unsupported flags
    echo "Error: Unsupported flag $SRA" >&2
    exit 1
    ;;
  *) # preserve positional arguments
    PARAMS="$PARAMS $SRA"
    shift
    ;;
  esac
done # set positional arguments in their proper place
eval set -- "$PARAMS"

# Check for required args
if [ -z "$SRA" ]; then
  echo "Error: Param missing --sra/-s SRX******** | SRA is required." && exit 0
fi
if [ -z "$GENOME" ]; then
  echo "Error: Param missing --genome/-g *** | Genome is required." && exit 0
fi

# Download SRA, genome
mkdir -p "./out" && cd "./out"
if ! test -f "$SRA"; then
  echo "Downloading reads.."
  curl "https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRA}/${SRA}" >"$SRA"
fi

if ! test -f "${SRA}_1.fastq.gz"; then
  echo "Analyzing reads.."
  parallel-fastq-dump --gzip --threads $N_THREADS --split-files --sra-id "$SRA"
fi

if ! test -f "${GENOME}.fa"; then
  echo "Downloading target genome.."
  curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id=${GENOME}" >"${GENOME}.fa"
fi

TARGET="${SRA}_${GENOME}_${ALIGNER}"
EXT="fastq.gz"
if [[ $TRIM == 1 ]]; then
  if ! test -f "${SRA}_1.paired.fq.gz"; then
    echo "Adaptor & Quality trimming.."
    trimmomatic PE ${SRA}_{1,2}.fastq.gz ${SRA}_{1,2}.{,un}paired.fq.gz AVGQUAL:20 HEADCROP:3 \
      LEADING:3 TRAILING:3 MINLEN:75 -threads ${N_THREADS}
  fi
  TARGET="${SRA}_${GENOME}_${ALIGNER}_trimmed"
  EXT="paired.fq.gz"
fi

case $ALIGNER in
bwa) align_bwa $SRA $GENOME $ALIGNER $EXT ;;
bowtie2) align_bowtie2 $SRA $GENOME $ALIGNER $EXT ;;
minimap2) align_minimap2 $SRA $GENOME $ALIGNER $EXT ;;
esac

samtools index $TARGET.bam
COV_PCT=$(samtools coverage ${TARGET}.bam -H | cut -f6)
echo "Coverage: ${COV_PCT}%"
echo "Finished!"
