#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: $0 <candidates.faa> <reference.faa_or_db> <output.m8> <tmp_dir> [threads]" >&2
  exit 1
fi

query=$1
reference=$2
output=$3
tmp_dir=$4
threads=${5:-64}

command -v mmseqs >/dev/null 2>&1 || {
  echo "Error: mmseqs is not available in PATH." >&2
  exit 1
}

mkdir -p "$(dirname "$output")" "$tmp_dir"

mmseqs easy-search \
  "$query" \
  "$reference" \
  "$output" \
  "$tmp_dir" \
  -e 1e-5 \
  --min-seq-id 0.8 \
  --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov,evalue,bits \
  --threads "$threads"

