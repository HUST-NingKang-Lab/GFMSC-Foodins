#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: bash $0 input.faa output.csv" >&2
  exit 1
fi
command -v toxinpred3 >/dev/null || { echo 'toxinpred3 is required.' >&2; exit 1; }
mkdir -p "$(dirname "$2")"
toxinpred3 -i "$1" -o "$2"
