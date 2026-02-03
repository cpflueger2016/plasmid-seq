#!/usr/bin/env bash
set -euo pipefail

# Defaults
# Defaults (populated by config)
DEFAULT_TSV=""
DEFAULT_REFS=""
MAX_CONCURRENT=""
log_file="" 


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
  -l <file>  Submit log file (default: <plasmidSeqData>/plasmidseq_submit_<date>.log)

Example:
  $(basename "$0") -d /home/.../fastqs
  $(basename "$0") -d /home/.../fastqs -t ./PL_to_fasta.tsv -f ./Fasta_Reference_Files -p 80
EOF
}


while getopts ":d:t:f:p:l:c:h" opt; do
  case "$opt" in
    d) plasmidSeqData="$OPTARG" ;;
    t) tsv="$OPTARG" ;;
    f) refs="$OPTARG" ;;
    p) MAX_CONCURRENT="$OPTARG" ;;
    c) PLASMIDSEQ_CONFIG="$OPTARG" ;;
    l) log_file="$OPTARG" ;;
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
if [[ ! -f "$tsv" ]]; then
  echo "[submit][ERROR] TSV not found: $tsv" >&2
  exit 1
fi
if [[ ! -d "$refs" ]]; then
  echo "[submit][ERROR] Refs dir not found: $refs" >&2
  exit 1
fi


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
  echo "[submit] config=$cfg"
}

jobdate="$(date +%Y-%m-%d)"
jobname="plasmidSeq_${jobdate}"

# Parse args first (so -c can point to a different config)

load_config

# Define where the scripts location are
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cfg="${PLASMIDSEQ_CONFIG:-${script_dir}/plasmidseq.config}"
if [[ -f "$cfg" ]]; then
  # shellcheck source=/dev/null
  source "$cfg"
fi

# Config precedence:
#   1) -c <file> (sets PLASMIDSEQ_CONFIG)
#   2) env var PLASMIDSEQ_CONFIG
#   3) scripts/plasmidseq.local.config (if present)
#   4) scripts/plasmidseq.config (default)
CONFIG_DEFAULT="${script_dir}/plasmidseq.config"
CONFIG_LOCAL="${script_dir}/plasmidseq.local.config"

# Apply config defaults only if user didn't override via CLI
tsv="${tsv:-${DEFAULT_TSV}}"
refs="${refs:-${DEFAULT_REFS}}"
MAX_CONCURRENT="${MAX_CONCURRENT:-${MAX_CONCURRENT_DEFAULT:-50}}"
RESULTS_BASE="${RESULTS_BASE:-/group/llshared/PlasmidSeq/Results}"
MAX_ARRAY_SIZE="${MAX_ARRAY_SIZE:-1001}"

# --- logging: mirror stdout/stderr to a file via tee
if [[ -z "${log_file}" ]]; then
  log_file="${plasmidSeqData%/}/plasmidseq_submit_${jobdate}.log"
fi
mkdir -p "$(dirname "$log_file")"
exec > >(tee -a "$log_file") 2>&1
echo "[submit] logging to: $log_file"


echo "[submit] jobdate=$jobdate"
echo "[submit] plasmidSeqData=$plasmidSeqData"
echo "[submit] TSV=$tsv"
echo "[submit] REFS=$refs"
echo "[submit] max_concurrent=$MAX_CONCURRENT"

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

# best effort: keep a copy of the submit log alongside the run in scratch
# (gather can then copy it into RESULTS automatically)
mkdir -p "${SCRATCH}" 2>/dev/null || true
cp -f "$log_file" "${SCRATCH}/plasmidseq_submit.log" 2>/dev/null || true

echo "[submit] expected SCRATCH=$SCRATCH"
echo "[submit] expected RESULTS=$RESULTS"

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
array_jobids=()

while [[ $offset -lt $n_jobs ]]; do
  remaining=$((n_jobs - offset))
  this_chunk=$(( remaining < MAX_ARRAY ? remaining : MAX_ARRAY ))
  last_index=$((this_chunk - 1))

  array_jobid=$(sbatch --parsable \
    --dependency=afterok:"$prep_jobid" \
    --array=0-${last_index}%${MAX_CONCURRENT} \
    "${script_dir}/plasmidseq_map_array_SLURM.sh" \
    "$SCRATCH" "$offset"
  )

  echo "[submit] map array chunk $chunk: job=$array_jobid offset=$offset size=$this_chunk"
  array_jobids+=("$array_jobid")

  offset=$((offset + this_chunk))
  chunk=$((chunk + 1))
done

# 4) Submit gather job — depends on array completing successfully
gather_jobid=$(sbatch --parsable \
  --dependency=afterok:"$array_jobid" \
  "${script_dir}/plasmidseq_gather_SLURM.sh" \
  "$SCRATCH" "$RESULTS" "$plasmidSeqData"
)
echo "[submit] gather job: $gather_jobid"

echo "[submit] done."
echo "[submit] Track with: squeue -j ${prep_jobid},${array_jobid},${gather_jobid}"
