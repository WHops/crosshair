#!/bin/bash

# Usage:
#   ./extract_region.sh "utg000179l:1116501-1136500_-" reference.fa
#
# Description:
#   Extracts a sequence from the given FASTA file using bedtools getfasta,
#   reverse-complementing if the strand is '-'.

# Show usage if arguments are missing or help is requested
if [[ $# -ne 2 || "$1" == "-h" || "$1" == "--help" ]]; then
  echo "Usage: $0 \"<region_string>\" <reference_fasta>"
  echo "Example: $0 \"utg000179l:1116501-1136500_-\" genome.fa"
  echo
  echo "Extracts a region from a FASTA file with strand-aware reverse complementing."
  echo "Requires: bedtools"
  exit 1
fi

region="$1"
fasta="$2"

# Parse components
seqid=$(echo "$region" | cut -d':' -f1)
coords=$(echo "$region" | cut -d':' -f2 | cut -d'_' -f1)
strand=$(echo "$region" | grep -oE '[+-]$')
start=$(echo "$coords" | cut -d'-' -f1)
end=$(echo "$coords" | cut -d'-' -f2)

# Convert to BED format (0-based, half-open)
echo -e "$seqid\t$((start - 1))\t$end\t$region\t0\t$strand" > tmp_region.bed

# Extract the sequence using bedtools
bedtools getfasta -fi "$fasta" -bed tmp_region.bed -s -name

# Clean up
rm tmp_region.bed

