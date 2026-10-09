#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <assembly_directory> <output_directory> [threads]" >&2
  exit 1
fi

input_dir=$1
output_dir=$2
threads=${3:-1}

command -v pyrodigal >/dev/null 2>&1 || {
  echo "Error: pyrodigal is not available in PATH." >&2
  exit 1
}

mkdir -p "$output_dir"

# Pyrodigal itself processes one assembly per invocation. The optional thread
# count controls how many assemblies are processed concurrently.
find "$input_dir" -type f \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) -print0 |
  xargs -0 -n 1 -P "$threads" bash -c '
    assembly=$1
    output_dir=$2
    sample=$(basename "$assembly")
    sample=${sample%.fasta}
    sample=${sample%.fna}
    sample=${sample%.fa}

    pyrodigal \
      -i "$assembly" \
      -a "$output_dir/${sample}.proteins.faa" \
      -p meta \
      --min-gene 15 \
      --min-edge-gene 15 \
      --max-overlap 15
  ' _ {} "$output_dir"

