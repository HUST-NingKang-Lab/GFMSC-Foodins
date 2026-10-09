#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <GFMSC_representatives.faa> <output_dir> <tmp_root> <database> [database ...]" >&2
  exit 1
fi

query=$1
output_dir=$2
tmp_root=$3
shift 3

command -v mmseqs >/dev/null 2>&1 || {
  echo "Error: mmseqs is not available in PATH." >&2
  exit 1
}

mkdir -p "$output_dir" "$tmp_root"

for database in "$@"; do
  database_name=$(basename "$database")
  output="$output_dir/${database_name}.m8"
  tmp_dir="$tmp_root/${database_name}"

  mmseqs easy-search \
    "$query" \
    "$database" \
    "$output" \
    "$tmp_dir" \
    -e 0.001 \
    --min-seq-id 0.3 \
    --cov-mode 2 \
    -c 0.8 \
    --max-accept 1 \
    --start-sens 1 \
    --sens-steps 3 \
    -s 7 \
    --threads 64
done

