#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: $0 <contigs.fasta> <GTDB_mmseqs_database> <output_prefix> <tmp_dir> [threads]" >&2
  exit 1
fi

contigs=$1
gtdb_db=$2
output_prefix=$3
tmp_dir=$4
threads=${5:-64}

command -v mmseqs >/dev/null 2>&1 || {
  echo "Error: mmseqs is not available in PATH." >&2
  exit 1
}

mkdir -p "$(dirname "$output_prefix")" "$tmp_dir"

# Taxonomy for cFMD-derived smORFs is inherited from the source MAG metadata.
# This command is intended for smORFs derived from assembled metagenomic contigs.
mmseqs easy-taxonomy \
  "$contigs" \
  "$gtdb_db" \
  "$output_prefix" \
  "$tmp_dir" \
  --threads "$threads"

