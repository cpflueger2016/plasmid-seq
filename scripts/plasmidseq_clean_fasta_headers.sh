#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  plasmidseq_clean_fasta_headers.sh -i <input.fa> -o <output.fa>

Behavior:
  - For each FASTA header, keeps only text before the first space.
  - Replaces all non-alphanumeric/underscore characters in the remaining header with "__".
EOF
}

in_fa=""
out_fa=""

while getopts ":i:o:h" opt; do
  case "$opt" in
    i) in_fa="$OPTARG" ;;
    o) out_fa="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 2 ;;
  esac
done

if [[ -z "$in_fa" || -z "$out_fa" ]]; then
  usage
  exit 2
fi

if [[ ! -f "$in_fa" ]]; then
  echo "[clean-fasta][ERROR] input FASTA not found: $in_fa" >&2
  exit 1
fi

awk '
  /^>/ {
    hdr = substr($0, 2)
    sub(/[[:space:]].*$/, "", hdr)
    gsub(/[^A-Za-z0-9_]/, "__", hdr)
    if (hdr == "") {
      hdr = "unnamed"
    }
    print ">" hdr
    next
  }
  { print }
' "$in_fa" > "$out_fa"

echo "[clean-fasta] input=$in_fa"
echo "[clean-fasta] output=$out_fa"
