#!/usr/bin/env bash
set -euo pipefail

# Robust matcher/copy helper:
#   - reads PL_to_fasta.tsv (tab-delimited, PL in col1, fasta name in col2)
#   - finds exactly ONE destination dir matching "${PL}*"
#   - finds exactly ONE reference fasta in fastaRefs (prefers exact basename match)
#   - copies fasta into destination dir
#   - otherwise touches marker files + creates "na"
#
# Usage:
#   match_plasmid_fasta_refs_robust.bash [-v] [-n] [-o] [-l logfile] [-r dest_root] <PL_to_fasta.tsv> <fasta_refs_dir>
#
# Examples:
#   ./match_plasmid_fasta_refs_robust.bash -v PL_to_fasta.tsv /path/to/refs  > match.log 2>&1
#   ./match_plasmid_fasta_refs_robust.bash -v -l match.log PL_to_fasta.tsv /path/to/refs
#   ./match_plasmid_fasta_refs_robust.bash -v -n PL_to_fasta.tsv /path/to/refs   # dry run

VERBOSE=0
DRYRUN=0
OVERWRITE=0
LOGFILE=""
DEST_ROOT="."

usage() {
  cat <<'EOF'
Usage:
  match_plasmid_fasta_refs_robust.bash [options] <PL_to_fasta.tsv> <fasta_refs_dir>

Options:
  -v            verbose (prints INFO lines; errors always print)
  -n            dry run (do not copy/touch anything)
  -o            overwrite existing fasta in destination (default: do not overwrite)
  -l <file>     also append logs to <file>
  -r <dir>      destination search root (default: current directory '.')

Notes:
  - Destination folder is searched as: find <dest_root> -type d -iname "${PL}*"
  - Reference fasta search prefers exact basename match first; then falls back to "${name}*"
  - Allowed reference extensions: .fa or .fasta (case-insensitive)
EOF
}

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

log() {
  local level="$1"; shift
  local msg="[$(timestamp)] [$level] $*"
  if [[ -n "$LOGFILE" ]]; then
    # append to logfile as well as stderr/stdout
    if [[ "$level" == "ERROR" || "$level" == "WARN" ]]; then
      echo "$msg" | tee -a "$LOGFILE" >&2
    else
      echo "$msg" | tee -a "$LOGFILE"
    fi
  else
    if [[ "$level" == "ERROR" || "$level" == "WARN" ]]; then
      echo "$msg" >&2
    else
      echo "$msg"
    fi
  fi
}

vlog() {
  if [[ "$VERBOSE" -eq 1 ]]; then
    log "INFO" "$@"
  fi
}

touch_safe() {
  local path="$1"
  if [[ "$DRYRUN" -eq 1 ]]; then
    vlog "DRYRUN touch $path"
    return 0
  fi
  touch "$path"
}

write_file_safe() {
  local path="$1"; shift
  if [[ "$DRYRUN" -eq 1 ]]; then
    vlog "DRYRUN write $path"
    return 0
  fi
  printf "%s\n" "$@" > "$path"
}

copy_safe() {
  local src="$1"
  local dst_dir="$2"
  local dst="$dst_dir/$(basename "$src")"

  if [[ "$DRYRUN" -eq 1 ]]; then
    vlog "DRYRUN cp '$src' -> '$dst_dir/'"
    return 0
  fi

  if [[ -e "$dst" && "$OVERWRITE" -eq 0 ]]; then
    vlog "Destination already has $(basename "$src"); not overwriting: $dst"
    touch_safe "$dst_dir/FASTA_REF_ALREADY_PRESENT"
    return 0
  fi

  cp -f "$src" "$dst_dir/"
}

trim() {
  # trim leading/trailing whitespace
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf "%s" "$s"
}

# ---- parse options
while getopts ":vno:l:r:" opt; do
  case "$opt" in
    v) VERBOSE=1 ;;
    n) DRYRUN=1 ;;
    o) OVERWRITE=1 ;;
    l) LOGFILE="$OPTARG" ;;
    r) DEST_ROOT="$OPTARG" ;;
    *) usage; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

if [[ $# -ne 2 ]]; then
  usage
  exit 2
fi

TSV="$1"
FASTA_REFS="$2"

if [[ ! -f "$TSV" ]]; then
  log "ERROR" "TSV file not found: $TSV"
  exit 1
fi
if [[ ! -d "$FASTA_REFS" ]]; then
  log "ERROR" "Reference folder not found: $FASTA_REFS"
  exit 1
fi
if [[ ! -d "$DEST_ROOT" ]]; then
  log "ERROR" "Destination root not found: $DEST_ROOT"
  exit 1
fi

vlog "Starting. TSV='$TSV' refs='$FASTA_REFS' dest_root='$DEST_ROOT' dryrun=$DRYRUN overwrite=$OVERWRITE"

# ---- main loop
# read: PL (col1) and fasta name (col2), ignore other cols
while IFS=$'\t' read -r pl_raw fasta_raw _rest || [[ -n "${pl_raw:-}" ]]; do
  # skip empty lines
  [[ -z "${pl_raw:-}" && -z "${fasta_raw:-}" ]] && continue

  # remove Windows CR if present + trim
  pl_raw="${pl_raw//$'\r'/}"
  fasta_raw="${fasta_raw//$'\r'/}"
  pl="$(trim "${pl_raw:-}")"
  fasta_name="$(trim "${fasta_raw:-}")"

  # skip comments
  [[ "${pl:0:1}" == "#" ]] && continue

  if [[ -z "$pl" ]]; then
    log "WARN" "Skipping row with empty PL (fasta='$fasta_name')"
    continue
  fi

  # Find destination directories matching PL*
  mapfile -t dests < <(find "$DEST_ROOT" -type d -iname "${pl}*" -print | sort)

  if [[ ${#dests[@]} -eq 0 ]]; then
    log "ERROR" "PL=$pl : No destination directory found under '$DEST_ROOT' matching '${pl}*' (fasta='$fasta_name')"
    continue
  fi

  if [[ ${#dests[@]} -gt 1 ]]; then
    log "ERROR" "PL=$pl : Ambiguous destination directories (${#dests[@]} matches). Will mark all and skip."
    for d in "${dests[@]}"; do
      log "ERROR" "  dest_match: $d"
      touch_safe "$d/DEST_DIR_AMBIGUOUS"
    done
    continue
  fi

  dest="$(readlink -f "${dests[0]}")"

  # If fasta column is empty / na
  if [[ -z "$fasta_name" || "$fasta_name" == "na" || "$fasta_name" == "NA" ]]; then
    log "WARN" "PL=$pl : No fasta provided; creating '$dest/na'"
    touch_safe "$dest/na"
    touch_safe "$dest/FASTA_REF_NOT_SPECIFIED"
    continue
  fi

  # Validate fasta extension (must be .fa or .fasta)
  shopt -s nocasematch
  if [[ ! "$fasta_name" =~ \.(fa|fasta)$ ]]; then
    shopt -u nocasematch
    log "WARN" "PL=$pl : fasta='$fasta_name' does not end in .fa/.fasta; creating '$dest/na'"
    touch_safe "$dest/na"
    touch_safe "$dest/FASTA_REF_INVALID_EXTENSION"
    continue
  fi
  shopt -u nocasematch

  # Find reference fasta candidates
  # Prefer exact basename match first, then allow prefix matches (like original script "${f}*")
  mapfile -t candidates_all < <(find "$FASTA_REFS" -type f \( -iname "$fasta_name" -o -iname "${fasta_name}*" \) -print | sort)

  # Filter to allowed extensions only
  candidates=()
  for c in "${candidates_all[@]}"; do
    shopt -s nocasematch
    if [[ "$c" =~ \.(fa|fasta)$ ]]; then
      candidates+=("$c")
    fi
    shopt -u nocasematch
  done

  exact=()
  for c in "${candidates[@]}"; do
    bn="$(basename "$c")"
    shopt -s nocasematch
    if [[ "$bn" == "$fasta_name" ]]; then
      exact+=("$c")
    fi
    shopt -u nocasematch
  done

  ref=""
  if [[ ${#exact[@]} -eq 1 ]]; then
    ref="${exact[0]}"
    vlog "PL=$pl : exact match found: $ref"
  elif [[ ${#exact[@]} -gt 1 ]]; then
    log "ERROR" "PL=$pl : Multiple EXACT reference matches for '$fasta_name' (${#exact[@]}). Marking ambiguous."
    touch_safe "$dest/na"
    touch_safe "$dest/FASTA_REF_AMBIGUOUS"
    write_file_safe "$dest/FASTA_REF_AMBIGUOUS.matches.txt" "${exact[@]}"
    continue
  elif [[ ${#candidates[@]} -eq 1 ]]; then
    ref="${candidates[0]}"
    vlog "PL=$pl : single prefix match found: $ref"
  elif [[ ${#candidates[@]} -gt 1 ]]; then
    log "ERROR" "PL=$pl : Multiple reference matches for '$fasta_name*' (${#candidates[@]}). Marking ambiguous."
    touch_safe "$dest/na"
    touch_safe "$dest/FASTA_REF_AMBIGUOUS"
    write_file_safe "$dest/FASTA_REF_AMBIGUOUS.matches.txt" "${candidates[@]}"
    continue
  else
    log "WARN" "PL=$pl : Reference fasta not found for '$fasta_name' in '$FASTA_REFS'. Creating '$dest/na'"
    touch_safe "$dest/na"
    touch_safe "$dest/FASTA_REF_NOT_FOUND"
    continue
  fi

  # Copy and record provenance
  log "INFO" "PL=$pl : Copying ref='${ref}' -> dest='${dest}/'"
  copy_safe "$ref" "$dest"
  write_file_safe "$dest/FASTA_REF_SOURCE.txt" \
    "PL: $pl" \
    "Requested fasta name: $fasta_name" \
    "Selected reference path: $ref" \
    "Timestamp: $(timestamp)"

  touch_safe "$dest/FASTA_REF_FOUND"

done < "$TSV"

vlog "Done."
