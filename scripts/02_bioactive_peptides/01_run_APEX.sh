#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <APEX_predict.py> <GFMSC_representatives.faa> [gpu_id]" >&2
  exit 1
fi

apex_script=$1
input_fasta=$2
gpu_id=${3:-1}

python "$apex_script" -i "$input_fasta" -g "$gpu_id"

