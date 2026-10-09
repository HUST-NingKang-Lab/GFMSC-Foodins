#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 5 ]]; then
  echo "Usage: $0 <cFMD_smORFs.faa> <fermented_grain_smORFs.faa> <output_dir> [threads] [memory_mb]" >&2
  exit 1
fi

cfmd_fasta=$1
grain_fasta=$2
output_dir=$3
threads=${4:-64}
memory_mb=${5:-64000}

command -v cd-hit >/dev/null 2>&1 || {
  echo "Error: cd-hit is not available in PATH." >&2
  exit 1
}

mkdir -p "$output_dir"
combined="$output_dir/GFMSC_all_smORFs.faa"
representatives="$output_dir/GFMSC_representatives.faa"

cat "$cfmd_fasta" "$grain_fasta" > "$combined"

# This preserves the command recorded in FBP.md. Note that -c 0.5 specifies
# 50% sequence identity; -s 0.95 and -aL 0.95 are coverage constraints.
cd-hit \
  -i "$combined" \
  -o "$representatives" \
  -n 3 \
  -p 1 \
  -c 0.5 \
  -d 0 \
  -M "$memory_mb" \
  -s 0.95 \
  -aL 0.95 \
  -g 1 \
  -T "$threads"

