#!/usr/bin/env bash
# Merge compatible raw-count Cooler matrices into one single-resolution .cool.
# Usage: bash scripts/merge_microc_replicates.sh OUTPUT.cool 400 INPUT1.mcool::/resolutions/400 INPUT2...
set -euo pipefail
if [[ $# -lt 4 ]]; then
  echo 'Usage: bash scripts/merge_microc_replicates.sh OUTPUT.cool RESOLUTION INPUT1.cool_or_URI INPUT2.cool_or_URI [INPUT3...]' >&2
  exit 2
fi
OUTPUT="$1"
RESOLUTION="$2"
shift 2
INPUTS=("$@")
if ! [[ "$RESOLUTION" =~ ^[1-9][0-9]*$ ]]; then
  echo 'RESOLUTION must be a positive integer' >&2
  exit 2
fi
for command_name in python cooler; do
  command -v "$command_name" >/dev/null || { echo "Missing command: $command_name" >&2; exit 1; }
done
for input in "${INPUTS[@]}"; do
  base_path="${input%%::*}"
  [[ -f "$base_path" ]] || { echo "Input file not found: $base_path" >&2; exit 1; }
done
[[ ! -e "$OUTPUT" ]] || { echo "Output already exists; will not overwrite: $OUTPUT" >&2; exit 1; }
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$SCRIPT_DIR/check_replicates.py" --expected-resolution "$RESOLUTION" "${INPUTS[@]}"
mkdir -p "$(dirname "$OUTPUT")"
echo 'Merging raw contact counts from compatible samples...'
cooler merge "$OUTPUT" "${INPUTS[@]}"
python "$SCRIPT_DIR/inspect_cool.py" "$OUTPUT" --expected-resolution "$RESOLUTION"
echo "Merged output: $OUTPUT"
echo 'Normalization weights were not propagated. Run cooler balance separately if needed.'
