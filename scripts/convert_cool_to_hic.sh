#!/usr/bin/env bash
# Run under Linux / WSL in an environment with Cooler and Java.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/convert_cool_to_hic.sh INPUT.cool OUTPUT.hic dm6.chrom.sizes juicer_tools.jar [RESOLUTION] [--include-trans]

Defaults:
  RESOLUTION = 400
  cis-only contacts (same-chromosome interactions only)

Optional environment variables:
  JAVA_HEAP=12g     Memory limit for Java (adapt to your RAM)
  SORT_MEMORY=1G    Memory allowed for GNU sort
  TMP_BASE=/path    Directory with substantial free disk space
EOF
}
if [[ $# -lt 4 || $# -gt 6 ]]; then usage; exit 2; fi
INPUT="$1"; OUTPUT="$2"; CHROMS="$3"; JAR="$4"
RESOLUTION="${5:-400}"
INCLUDE_TRANS="${6:-}"
if [[ "$RESOLUTION" == "--include-trans" && -z "$INCLUDE_TRANS" ]]; then
  RESOLUTION=400
  INCLUDE_TRANS="--include-trans"
fi
if [[ -n "$INCLUDE_TRANS" && "$INCLUDE_TRANS" != "--include-trans" ]]; then
  echo "Unrecognized option: $INCLUDE_TRANS" >&2; usage; exit 2
fi
if [[ ! "$RESOLUTION" =~ ^[0-9]+$ || "$RESOLUTION" == 0 ]]; then
  echo "RESOLUTION must be a positive integer" >&2; exit 2
fi
for path in "$INPUT" "$CHROMS" "$JAR"; do
  [[ -f "$path" ]] || { echo "Missing required file: $path" >&2; exit 1; }
done
for cmd in python cooler java gzip sort; do
  command -v "$cmd" >/dev/null || { echo "Missing program: $cmd" >&2; exit 1; }
done
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$(dirname "$OUTPUT")"
python "$SCRIPT_DIR/inspect_cool.py" "$INPUT" --expected-resolution "$RESOLUTION"
TMP_BASE="${TMP_BASE:-${TMPDIR:-/tmp}}"
mkdir -p "$TMP_BASE"
TMP_WORK="$(mktemp -d "$TMP_BASE/microc2hic.XXXXXXXX")"
trap 'rm -rf "$TMP_WORK"' EXIT
CIS_OPTIONS=()
if [[ "$INCLUDE_TRANS" != "--include-trans" ]]; then CIS_OPTIONS=(--cis-only); fi

printf 'Step 1/2: Stream unbalanced counts, convert to Juicer format, sort and compress...\n'
cooler dump --join --no-balance "$INPUT" \
  | python "$SCRIPT_DIR/convert_contacts.py" --chrom-sizes "$CHROMS" "${CIS_OPTIONS[@]}" \
  | LC_ALL=C sort -T "$TMP_WORK" -S "${SORT_MEMORY:-1G}" -k2,2 -k6,6 -k3,3n -k7,7n \
  | gzip -c > "$TMP_WORK/contacts.shortWithScore.sorted.txt.gz"
if [[ ! -s "$TMP_WORK/contacts.shortWithScore.sorted.txt.gz" ]]; then
  echo "No contacts written; check the source and chromosome names" >&2; exit 1
fi

printf 'Step 2/2: Build .hic using Juicer Tools...\n'
JAVA_HEAP="${JAVA_HEAP:-12g}"
JUICER_OPTIONS=(-r "$RESOLUTION")
if [[ "$INCLUDE_TRANS" != "--include-trans" ]]; then JUICER_OPTIONS=(-d "${JUICER_OPTIONS[@]}"); fi
java -Xmx"$JAVA_HEAP" -jar "$JAR" pre "${JUICER_OPTIONS[@]}" \
  "$TMP_WORK/contacts.shortWithScore.sorted.txt.gz" "$OUTPUT" "$CHROMS"
echo "Created: $OUTPUT"
