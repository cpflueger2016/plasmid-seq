#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") -d <fastq_dir> -t <pl_to_fasta.tsv> [options]

Required:
  -d <dir>    FASTQ source directory (submit -d)
  -t <file>   PL_to_fasta TSV

Optional:
  -w <file>   Plate map CSV
  -c <file>   Config file for plasmidseq_submit.sh
  -f <dir>    Reference folder override
  -p <int>    Max concurrent array tasks (default: 100)
  -C <list>   Core list, comma-separated (default: 1,2,3,4,6,8)
  -r <int>    Replicates per core/mode (default: 3)
  -m <list>   Modes: baseline,variants (default: baseline,variants)
  -o <dir>    Output benchmark dir (default: ./benchmark/runs/<timestamp>)
  -n          Dry run (print commands only)

Notes:
  baseline mode: no variant flags
  variants mode: adds -V -E
USAGE
}

FASTQ_DIR=""
TSV=""
PLATE_MAP=""
CONFIG=""
REFS=""
MAX_CONCURRENT="100"
CORE_LIST="1,2,3,4,6,8"
REPS="3"
MODE_LIST="baseline,variants"
OUTDIR=""
DRY_RUN=0

while getopts ":d:t:w:c:f:p:C:r:m:o:nh" opt; do
  case "$opt" in
    d) FASTQ_DIR="$OPTARG" ;;
    t) TSV="$OPTARG" ;;
    w) PLATE_MAP="$OPTARG" ;;
    c) CONFIG="$OPTARG" ;;
    f) REFS="$OPTARG" ;;
    p) MAX_CONCURRENT="$OPTARG" ;;
    C) CORE_LIST="$OPTARG" ;;
    r) REPS="$OPTARG" ;;
    m) MODE_LIST="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    n) DRY_RUN=1 ;;
    h) usage; exit 0 ;;
    *) usage; exit 2 ;;
  esac
done

if [[ -z "$FASTQ_DIR" || -z "$TSV" ]]; then
  usage
  exit 2
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "[bench][ERROR] FASTQ dir not found: $FASTQ_DIR" >&2
  exit 1
fi
if [[ ! -f "$TSV" ]]; then
  echo "[bench][ERROR] TSV not found: $TSV" >&2
  exit 1
fi
if [[ -n "$PLATE_MAP" && ! -f "$PLATE_MAP" ]]; then
  echo "[bench][ERROR] Plate map not found: $PLATE_MAP" >&2
  exit 1
fi
if [[ -n "$CONFIG" && ! -f "$CONFIG" ]]; then
  echo "[bench][ERROR] Config not found: $CONFIG" >&2
  exit 1
fi
if [[ -n "$REFS" && ! -d "$REFS" ]]; then
  echo "[bench][ERROR] Refs dir not found: $REFS" >&2
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
submit_script="$(cd "$script_dir/.." && pwd -P)/plasmidseq_submit.sh"
if [[ ! -x "$submit_script" ]]; then
  echo "[bench][ERROR] submit script not executable: $submit_script" >&2
  exit 1
fi

if [[ -z "$OUTDIR" ]]; then
  ts="$(date +%Y%m%d_%H%M%S)"
  OUTDIR="$script_dir/runs/$ts"
fi
mkdir -p "$OUTDIR/logs"

manifest="$OUTDIR/benchmark_manifest.tsv"
if [[ ! -f "$manifest" ]]; then
  cat > "$manifest" <<MANIFEST
run_label	mode	cores	replicate	submit_start_iso	submit_end_iso	submit_exit	prep_jobid	map_jobids	gather_jobid	scratch	results	submit_log
MANIFEST
fi

IFS=',' read -r -a CORES <<< "$CORE_LIST"
IFS=',' read -r -a MODES <<< "$MODE_LIST"

awk_is_int='^[0-9]+$'
if ! [[ "$REPS" =~ $awk_is_int ]] || [[ "$REPS" -lt 1 ]]; then
  echo "[bench][ERROR] invalid replicates: $REPS" >&2
  exit 1
fi

for mode in "${MODES[@]}"; do
  if [[ "$mode" != "baseline" && "$mode" != "variants" ]]; then
    echo "[bench][ERROR] invalid mode: $mode (allowed: baseline,variants)" >&2
    exit 1
  fi
done

for c in "${CORES[@]}"; do
  if ! [[ "$c" =~ $awk_is_int ]] || [[ "$c" -lt 1 ]]; then
    echo "[bench][ERROR] invalid core count: $c" >&2
    exit 1
  fi
done

for mode in "${MODES[@]}"; do
  for c in "${CORES[@]}"; do
    for rep in $(seq 1 "$REPS"); do
      run_label="${mode}_c${c}_r${rep}"
      log_file="$OUTDIR/logs/${run_label}.submit.log"
      tmp_out="$OUTDIR/logs/${run_label}.stdout.log"

      flags=()
      if [[ "$mode" == "variants" ]]; then
        flags+=("-V" "-E")
      fi

      cmd=("$submit_script" "-d" "$FASTQ_DIR" "-t" "$TSV" "-p" "$MAX_CONCURRENT" "-J" "$c")
      if [[ -n "$PLATE_MAP" ]]; then cmd+=("-w" "$PLATE_MAP"); fi
      if [[ -n "$CONFIG" ]]; then cmd+=("-c" "$CONFIG"); fi
      if [[ -n "$REFS" ]]; then cmd+=("-f" "$REFS"); fi
      cmd+=("-l" "$log_file")
      cmd+=("${flags[@]}")

      start_iso="$(date --iso-8601=seconds)"
      echo "[bench] run=$run_label mode=$mode cores=$c rep=$rep start=$start_iso"
      echo "[bench] cmd: ${cmd[*]}"

      if [[ "$DRY_RUN" -eq 1 ]]; then
        echo -e "${run_label}\t${mode}\t${c}\t${rep}\t${start_iso}\t\t0\t\t\t\t\t\t${log_file}" >> "$manifest"
        continue
      fi

      set +e
      "${cmd[@]}" > "$tmp_out" 2>&1
      submit_rc=$?
      set -e

      end_iso="$(date --iso-8601=seconds)"

      prep_jobid="$(rg -o "\\[submit\\] prep job: [0-9]+" "$tmp_out" | tail -n1 | awk '{print $4}' || true)"
      map_jobids="$(rg -o "map array chunk [0-9]+: job=[0-9]+" "$tmp_out" | awk -F'job=' '{print $2}' | paste -sd: - || true)"
      gather_jobid="$(rg -o "\\[submit\\] gather job: [0-9]+" "$tmp_out" | tail -n1 | awk '{print $4}' || true)"
      scratch="$(rg -o "\\[submit\\] expected SCRATCH=.*" "$tmp_out" | tail -n1 | sed 's/^\[submit\] expected SCRATCH=//' || true)"
      results="$(rg -o "\\[submit\\] expected RESULTS=.*" "$tmp_out" | tail -n1 | sed 's/^\[submit\] expected RESULTS=//' || true)"

      echo -e "${run_label}\t${mode}\t${c}\t${rep}\t${start_iso}\t${end_iso}\t${submit_rc}\t${prep_jobid}\t${map_jobids}\t${gather_jobid}\t${scratch}\t${results}\t${log_file}" >> "$manifest"

      if [[ "$submit_rc" -ne 0 ]]; then
        echo "[bench][WARN] submit failed for $run_label (rc=$submit_rc); see $tmp_out"
      else
        echo "[bench] run=$run_label prep=${prep_jobid:-na} map=${map_jobids:-na} gather=${gather_jobid:-na}"
      fi
    done
  done
done

echo "[bench] manifest: $manifest"
