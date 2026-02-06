#!/usr/bin/env bash
set -euo pipefail

# Defaults
# Defaults (populated by config)
DEFAULT_TSV=""
DEFAULT_REFS=""
MAX_CONCURRENT=""
PLASMIDSEQ_VERSION="1.0"
log_file="" 
plasmidSeqData=""
tsv=""
refs=""
PLASMIDSEQ_CONFIG=""
PLATE_MAP_CSV=""
CLI_ENABLE_VARIANTS=""
CLI_ENABLE_SNPEFF=""


usage() {
  cat <<EOF
Usage:
  $(basename "$0") -d <plasmidSeqData_dir> [options]

Required:
  -d <dir>   Folder containing Project folders / fastqs to stage into scratch

Optional:
  -t <file>  PL_to_fasta.tsv path (default: ${DEFAULT_TSV})
  -f <dir>   Fasta reference folder path (default: ${DEFAULT_REFS})
  -p <int>   Max concurrent array tasks (default: ${MAX_CONCURRENT})
  -c <file>  Config file path (default precedence: local config, then plasmidseq.config)
  -w <file>  Plate map CSV for run summary (columns: PLid,plate,position)
  -l <file>  Submit log file (default: <plasmidSeqData>/plasmidseq_submit_<date>.log)
  -v         Print version and exit
  -V         Enable VarScan variant calling (overrides config)
  -N         Disable VarScan variant calling (overrides config)
  -E         Enable snpEff annotation (overrides config; implies variants on)
  -S         Disable snpEff annotation (overrides config)

Example:
  $(basename "$0") -d /home/.../fastqs
  $(basename "$0") -d /home/.../fastqs -t ./PL_to_fasta.tsv -f ./Fasta_Reference_Files -p 80
  $(basename "$0") -d /home/.../fastqs -c ./plasmidseq.local.config -l ./submit.log
  $(basename "$0") -d /home/.../fastqs -w /group/llshared/PlasmidSeq/PL_to_plate_position.csv
  $(basename "$0") -d /home/.../fastqs -V -E
EOF
}


while getopts ":d:t:f:p:l:c:w:vVNESh" opt; do
  case "$opt" in
    d) plasmidSeqData="$OPTARG" ;;
    t) tsv="$OPTARG" ;;
    f) refs="$OPTARG" ;;
    p) MAX_CONCURRENT="$OPTARG" ;;
    c) PLASMIDSEQ_CONFIG="$OPTARG" ;;
    w) PLATE_MAP_CSV="$OPTARG" ;;
    l) log_file="$OPTARG" ;;
    v) echo "plasmid-seq version ${PLASMIDSEQ_VERSION}"; exit 0 ;;
    V) CLI_ENABLE_VARIANTS="1" ;;
    N) CLI_ENABLE_VARIANTS="0" ;;
    E) CLI_ENABLE_SNPEFF="1" ;;
    S) CLI_ENABLE_SNPEFF="0" ;;
    h) usage; exit 0 ;;
    *) usage; exit 2 ;;
  esac
done

if [[ -z "${plasmidSeqData}" ]]; then
  usage
  exit 2
fi

if [[ ! -d "$plasmidSeqData" ]]; then
  echo "[submit][ERROR] plasmidSeqData dir not found: $plasmidSeqData" >&2
  exit 1
fi


# Config precedence:
#   1) -c <file> (sets PLASMIDSEQ_CONFIG)
#   2) env var PLASMIDSEQ_CONFIG
#   3) scripts/plasmidseq.local.config (if present)
#   4) scripts/plasmidseq.config (default)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
CONFIG_DEFAULT="${script_dir}/plasmidseq.config"
CONFIG_LOCAL="${script_dir}/plasmidseq.local.config"

# Load the config file
load_config() {
  local cfg="${PLASMIDSEQ_CONFIG:-}"

  if [[ -z "$cfg" ]]; then
    if [[ -f "$CONFIG_LOCAL" ]]; then
      cfg="$CONFIG_LOCAL"
    else
      cfg="$CONFIG_DEFAULT"
    fi
  fi

  if [[ ! -f "$cfg" ]]; then
    echo "[submit][ERROR] Config file not found: $cfg" >&2
    exit 1
  fi

  # shellcheck source=/dev/null
  source "$cfg"
  # Keep the resolved config path for downstream sbatch jobs.
  PLASMIDSEQ_CONFIG="$cfg"
  export PLASMIDSEQ_CONFIG
  echo "[submit] config=$PLASMIDSEQ_CONFIG"
}

jobdate="$(date +%Y-%m-%d)"
jobname="plasmidSeq_${jobdate}"

# Parse args first (so -c can point to a different config)

load_config

# Apply config defaults only if user didn't override via CLI
tsv="${tsv:-${DEFAULT_TSV}}"
refs="${refs:-${DEFAULT_REFS}}"
MAX_CONCURRENT="${MAX_CONCURRENT:-${MAX_CONCURRENT_DEFAULT:-50}}"
RESULTS_BASE="${RESULTS_BASE:-/group/llshared/PlasmidSeq/Results}"
MAX_ARRAY_SIZE="${MAX_ARRAY_SIZE:-1001}"

# CLI flags override config defaults.
if [[ -n "${CLI_ENABLE_VARIANTS}" ]]; then
  ENABLE_VARIANTS="${CLI_ENABLE_VARIANTS}"
fi
if [[ -n "${CLI_ENABLE_SNPEFF}" ]]; then
  ENABLE_SNPEFF="${CLI_ENABLE_SNPEFF}"
fi
# snpEff requires variant calling output; auto-enable unless user explicitly forced off.
if [[ "${ENABLE_SNPEFF:-0}" == "1" && "${ENABLE_VARIANTS:-0}" != "1" ]]; then
  if [[ "${CLI_ENABLE_VARIANTS:-}" == "0" ]]; then
    echo "[submit][ERROR] Conflicting options: snpEff enabled (-E) but variants disabled (-N)." >&2
    exit 2
  fi
  ENABLE_VARIANTS="1"
  echo "[submit] enabling variants because snpEff is enabled"
fi
export ENABLE_VARIANTS ENABLE_SNPEFF

if [[ ! -f "$tsv" ]]; then
  echo "[submit][ERROR] TSV not found: $tsv" >&2
  exit 1
fi
if [[ ! -d "$refs" ]]; then
  echo "[submit][ERROR] Refs dir not found: $refs" >&2
  exit 1
fi
if [[ -n "$PLATE_MAP_CSV" && ! -f "$PLATE_MAP_CSV" ]]; then
  echo "[submit][ERROR] Plate map CSV not found: $PLATE_MAP_CSV" >&2
  exit 1
fi

# --- logging: mirror stdout/stderr to a file via tee
if [[ -z "${log_file}" ]]; then
  log_file="${plasmidSeqData%/}/plasmidseq_submit_${jobdate}.log"
fi
mkdir -p "$(dirname "$log_file")"
exec > >(tee -a "$log_file") 2>&1
echo "[submit] logging to: $log_file"


echo "[submit] jobdate=$jobdate"
echo "[submit] version=${PLASMIDSEQ_VERSION}"
echo "[submit] plasmidSeqData=$plasmidSeqData"
echo "[submit] TSV=$tsv"
echo "[submit] REFS=$refs"
echo "[submit] PLATE_MAP_CSV=${PLATE_MAP_CSV:-<none>}"
echo "[submit] max_concurrent=$MAX_CONCURRENT"
echo "[submit] ENABLE_VARIANTS=${ENABLE_VARIANTS:-0}"
echo "[submit] ENABLE_SNPEFF=${ENABLE_SNPEFF:-0}"

# 1) Submit prep job (writes jobs.tsv into scratch)
prep_jobid=$(sbatch --parsable \
  --export=ALL,PLASMIDSEQ_CONFIG="${PLASMIDSEQ_CONFIG:-}" \
  "${script_dir}/plasmidseq_prepare_SLURM.sh" \
  -d "$plasmidSeqData" -t "$tsv" -f "$refs" -j "$jobdate"
)
echo "[submit] prep job: $prep_jobid"

# We can compute the scratch path deterministically (since MYSCRATCH is stable and prep uses jobid)
SCRATCH="${MYSCRATCH}/${jobname}/${prep_jobid}"
RESULTS="${RESULTS_BASE}/${jobname}/${prep_jobid}"

# If snpEff is enabled, pass run-specific DB config/data paths to mapping jobs.
if [[ "${ENABLE_SNPEFF:-0}" == "1" ]]; then
  SNPEFF_DB="${SNPEFF_DB:-plasmidseq_${prep_jobid}}"
  SNPEFF_DATA_DIR="${SNPEFF_DATA_DIR:-${SCRATCH}/snpEff_data}"
  SNPEFF_CONFIG_FILE="${SNPEFF_CONFIG_FILE:-${SCRATCH}/snpEff.config}"
  export SNPEFF_DB SNPEFF_DATA_DIR SNPEFF_CONFIG_FILE
fi

# best effort: keep a copy of the submit log alongside the run in scratch
# (gather can then copy it into RESULTS automatically)
mkdir -p "${SCRATCH}" 2>/dev/null || true
mkdir -p "${SCRATCH}/Logs" 2>/dev/null || true
cp -f "$log_file" "${SCRATCH}/plasmidseq_submit.log" 2>/dev/null || true

echo "[submit] expected SCRATCH=$SCRATCH"
echo "[submit] expected RESULTS=$RESULTS"
if [[ "${ENABLE_SNPEFF:-0}" == "1" ]]; then
  echo "[submit] snpEff enabled: db=${SNPEFF_DB} config=${SNPEFF_CONFIG_FILE} dataDir=${SNPEFF_DATA_DIR}"
fi

# 2) Wait for jobs.tsv to appear (prep must finish staging and writing it)
jobs_file="${SCRATCH}/jobs.tsv"
echo "[submit] waiting for jobs file: $jobs_file"

# up to ~2 hours (720 * 10s); adjust if you want
for _ in $(seq 1 720); do
  # Success condition: jobs.tsv exists
  if [[ -f "$jobs_file" ]]; then
    break
  fi

  # Failure detection: if prep job finished in a bad state, stop waiting
  state=$(sacct -j "$prep_jobid" --noheader --format=State | head -n 1 | awk '{print $1}' || true)
  case "$state" in
    FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY)
      echo "[submit][ERROR] prep job $prep_jobid ended with state=$state; jobs.tsv will not appear." >&2
      echo "[submit][ERROR] Check: sacct -j $prep_jobid --format=JobID,State,ExitCode,Elapsed" >&2
      echo "[submit][ERROR] And check slurm output (e.g. slurm-${prep_jobid}.out) in the submission directory." >&2
      exit 1
      ;;
  esac

  sleep 10
done


if [[ ! -f "$jobs_file" ]]; then
  echo "[submit][ERROR] jobs.tsv did not appear at $jobs_file. Check prep logs for job $prep_jobid." >&2
  exit 1
fi

n_jobs=$(wc -l < "$jobs_file" | tr -d ' ')
if [[ "$n_jobs" -lt 1 ]]; then
  echo "[submit][ERROR] jobs.tsv exists but has 0 lines: $jobs_file" >&2
  exit 1
fi

echo "[submit] jobs: $n_jobs"

# 3) Submit mapping array (4 cores per plasmid task) — depends on prep
MAX_ARRAY=1001
chunk=0
offset=0
array_jobids_dep=""

while [[ $offset -lt $n_jobs ]]; do
  remaining=$((n_jobs - offset))
  this_chunk=$(( remaining < MAX_ARRAY ? remaining : MAX_ARRAY ))
  last_index=$((this_chunk - 1))

  array_jobid=$(sbatch --parsable \
    --dependency=afterok:"$prep_jobid" \
    --array=0-${last_index}%${MAX_CONCURRENT} \
    --output="$SCRATCH/Logs/slurm-%A_%a.out" \
    --error="$SCRATCH/Logs/slurm-%A_%a.out" \
    "${script_dir}/plasmidseq_map_array_SLURM.sh" \
    "$SCRATCH" "$offset"
  )
  echo "[submit] map array chunk $chunk: job=$array_jobid offset=$offset size=$this_chunk"
  if [[ -z "$array_jobids_dep" ]]; then
    array_jobids_dep="$array_jobid"
  else
    array_jobids_dep="${array_jobids_dep}:$array_jobid"
  fi

  offset=$((offset + this_chunk))
  chunk=$((chunk + 1))
done

# 4) Submit gather job — depends on array completing successfully
gather_jobid=$(sbatch --parsable \
  --dependency=afterok:"$array_jobids_dep" \
  --output="$SCRATCH/Logs/slurm-%A.out" \
  --error="$SCRATCH/Logs/slurm-%A.out" \
  "${script_dir}/plasmidseq_gather_SLURM.sh" \
  "$SCRATCH" "$RESULTS" "$plasmidSeqData" "${PLATE_MAP_CSV:-}"
)
echo "[submit] gather job: $gather_jobid"

echo "[submit] done."
echo "[submit] Track with: squeue -j ${prep_jobid},${array_jobids_dep//:/,},${gather_jobid}"
