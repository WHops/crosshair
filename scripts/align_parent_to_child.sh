#!/bin/bash

set -euo pipefail

# Function to display usage information
usage() {
    echo "Usage: $0 <child_fasta> <parent_fasta> <chunk_length> <child_locus>"
    echo "Arguments:"
    echo "  child_fasta    Path to the child FASTA file"
    echo "  parent_fasta   Path to the parent FASTA file"
    echo "  chunk_length   Length of the chunks to split the sequences into"
    echo "  child_locus    Genomic locus (chr:start-end) or 'all' to use entire FASTA"
    exit 1
}

# Check argument count
if [ "$#" -ne 4 ]; then
    usage
fi

# Input arguments
CHILDFA=$1
PARENTFA=$2
CHUNKLEN=$3
CHILDLOCUS=$4

# Create temporary directory and file paths
TMPDIR=$(mktemp -d)
SLIDE_FA="${TMPDIR}/parent_slide.fa"
HEADERS_AND_LENGTHS_TSV="${TMPDIR}/headers_and_lengths.tsv"
SLIDE_WITH_LENGTHS_FA="${TMPDIR}/slide_with_lengths.fa"
CHILD_SUBSEQ_FA="${TMPDIR}/child_subseq.fa"
BED_FILE="${TMPDIR}/region.bed"

# Ensure temporary files are cleaned up on exit
trap "rm -rf $TMPDIR" EXIT

# Handle CHILDLOCUS input
if [ "$CHILDLOCUS" = "all" ]; then
    CHILD_INPUT="$CHILDFA"
else
    # Parse CHILDLOCUS (format: chr:start-end)
    CHR=$(echo "$CHILDLOCUS" | cut -d: -f1)
    RANGE=$(echo "$CHILDLOCUS" | cut -d: -f2)
    START=$(echo "$RANGE" | cut -d- -f1)
    END=$(echo "$RANGE" | cut -d- -f2)

    # Create BED file and extract region
    echo -e "${CHR}\t${START}\t${END}" > "$BED_FILE"
    bedtools getfasta -fi "$CHILDFA" -bed "$BED_FILE" -name > "$CHILD_SUBSEQ_FA"

    CHILD_INPUT="$CHILD_SUBSEQ_FA"
fi

# Step 1: Generate sliding windows from parent FASTA
seqkit sliding -s "$CHUNKLEN" -W "$CHUNKLEN" "$PARENTFA" -g -S "" > "$SLIDE_FA"

# Step 2: Record sequence headers and lengths
seqkit fx2tab -n -l "$PARENTFA" > "$HEADERS_AND_LENGTHS_TSV"

# Step 3: Append original sequence length to headers
awk 'NR==FNR {len[$1]=$2; next} /^>/ {split($1, a, /[>:]/); print $1 "/" len[a[2]]; next} {print}' \
    "$HEADERS_AND_LENGTHS_TSV" "$SLIDE_FA" > "$SLIDE_WITH_LENGTHS_FA"

# Step 4: Map the parent chunks to the child input
minimap2 -x map-hifi -t 8 -c "$CHILD_INPUT" "$SLIDE_WITH_LENGTHS_FA"
