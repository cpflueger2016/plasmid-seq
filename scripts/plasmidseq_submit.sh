#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  $(basename "$0") -d <plasmidSeqData_dir> [options]

Required:
  -d <dir>    Folder containing Project folders / fastqs to stage into scratch

Optional:
  -t <file>   PL_to_fasta.tsv path (default from config: DEFAULT_TSV)
  -f <dir>    Fasta reference folder path (default from config: DEFAULT_REFS)
  -p <int>    Max concurrent array tasks per array chunk (default from config: MAX_CONCURRENT_DEFAULT)
  -l <file>   Submit log file (default: <plasmidSeqData>/plasmidseq_submit_<date>.log)
  -c <file>   Config file (default: scripts/plasmidseq.local.config if present, else scripts/plasmidseq.config)
  -h          Help

Examples:
  $(basename "$0") -d /path/to/fastqs
  $(basename "$0") -d /path/to/fastqs -p 100
  $(basename "$0") -d /path/to/fastqs -c ./plasmidseq.config
EOF
}

plasmidSeqData=""
tsv=""
refs=""
max_concurrent=""
log_file=""
PLASMIDSEQ_CONFIG=""

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
CONFIG_DEFAULT="${script_dir}/plasmidseq.config"
CONFIG_LOCAL="${script_dir}/plasmidseq.local.config"

while getopts ":d:t:f:p:l:c:h" opt; do
  case "$opt" in
    d) plasmidSeqData="$OPTARG" ;;
    t) tsv="$OPTARG" ;;
    f) refs="$OPTARG" ;;
    p) max_concurrent="$OPTARG" ;;
    l) log_file="$OPTARG" ;;
    c) PLASMIDSEQ_CONFIG="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 2 ;;
  esac
done

if [[ -z "${plasmidSeqData}" ]]; then
  usage
  exit 2
fi

# --- Load config ---
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
  export PLASMIDSEQ_CONFIG="$cfg"
}

load_config

# Apply config defaults if CLI not provided
tsv="${tsv:-${DEFAULT_TSV}}"
refs="${refs:-${DEFAULT_REFS}}"
max_concurrent="${max_concurrent:-${MAX_CONCURRENT_DEFAULT:-50}}"
RESULTS_BASE="${RESULTS_BASE:-/group/llshared/PlasmidSeq/Results}"
MAX_ARRAY_SIZE="${MAX_ARRAY_SIZE:-1001}"
SCRATCH_ROOT="${SCRATCH_BASE:-${MYSCRATCH:-}}"

jobdate="$(date +%Y-%m-%d)"
jobname="plasmidSeq_${jobdate}"

# default submit log file
if [[ -z "${log_file}" ]]; then
  log_file="${plasmidSeqData%/}/plasmidseq_submit_${jobdate}.log"
fi
mkdir -p "$(dirname "$log_file")"
exec > >(tee -a "$log_file") 2>&1

echo "[submit] logging to: $log_file"
echo "[submit] config=${PLASMIDSEQ_CONFIG}"
echo "[submit] jobdate=$jobdate"
echo "[submit] plasmidSeqData=$plasmidSeqData"
echo "[submit] TSV=$tsv"
echo "[submit] REFS=$refs"
echo "[submit] max_concurrent=$max_concurrent"
echo "[submit] max_array_size=$MAX_ARRAY_SIZE"

# --- Validation ---
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
if [[ -z "${SCRATCH_ROOT:-}" ]]; then
  echo "[submit][ERROR] SCRATCH_ROOT is empty (set MYSCRATCH env var or SCRATCH_BASE in config)." >&2
  exit 1
fi

# --- SBATCH common args ---
sbatch_common=()
if [[ -n "${SBATCH_ACCOUNT:-}" ]]; then sbatch_common+=(--account="${SBATCH_ACCOUNT}"); fi
if [[ -n "${SBATCH_QOS:-}" ]]; then sbatch_common+=(--qos="${SBATCH_QOS}"); fi

export_str="ALL,PLASMIDSEQ_CONFIG=${PLASMIDSEQ_CONFIG}"
# Keep existing env (export=ALL) plus our config path
sbatch_export=(--export="${export_str}")

# Resolve SLURM scripts (relative to this submit script by default)
resolve_slurm() {
  local base="$1"
  local override="$2"
  if [[ -n "$override" && -f "$override" ]]; then
    echo "$override"
  else
    echo "${script_dir}/${base}"
  fi
}

PREP_SLURM="$(resolve_slurm "${PREP_SLURM_BASENAME:-plasmidseq_prepare_SLURM.sh}" "${PREP_SLURM_PATH:-}")"
MAP_SLURM="$(resolve_slurm "${MAP_SLURM_BASENAME:-plasmidseq_map_array_SLURM.sh}" "${MAP_SLURM_PATH:-}")"
GATHER_SLURM="$(resolve_slurm "${GATHER_SLURM_BASENAME:-plasmidseq_gather_SLURM.sh}" "${GATHER_SLURM_PATH:-}")"

# --- 1) Prep job ---
prep_args=("${sbatch_common[@]}" "${sbatch_export[@]}")
if [[ -n "${PARTITION_PREP:-}" ]]; then prep_args+=(--partition="${PARTITION_PREP}"); fi

prep_jobid=$(sbatch --parsable "${prep_args[@]}"   "$PREP_SLURM"   -d "$plasmidSeqData" -t "$tsv" -f "$refs" -j "$jobdate"
)
echo "[submit] prep job: $prep_jobid"

SCRATCH="${SCRATCH_ROOT}/${jobname}/${prep_jobid}"
RESULTS="${RESULTS_BASE}/${jobname}/${prep_jobid}"

echo "[submit] expected SCRATCH=$SCRATCH"
echo "[submit] expected RESULTS=$RESULTS"

# --- 2) Wait for jobs.tsv ---
jobs_file="${SCRATCH}/jobs.tsv"
echo "[submit] waiting for jobs file: $jobs_file"

get_state() {
  sacct -n -P -j "$1" --format=State 2>/dev/null | head -n1 | cut -d'|' -f1
}

for i in $(seq 1 720); do
  if [[ -f "$jobs_file" ]]; then
    break
  fi

  # every ~60s, check if prep has failed
  if (( i % 6 == 0 )); then
    st="$(get_state "$prep_jobid" || true)"
    if [[ "$st" == "FAILED" || "$st" == "CANCELLED" || "$st" == "TIMEOUT" ]]; then
      echo "[submit][ERROR] prep job $prep_jobid ended with state=$st; jobs.tsv will not appear." >&2
      echo "[submit][ERROR] Check: sacct -j $prep_jobid --format=JobID,State,ExitCode,Elapsed" >&2
      exit 1
    fi
  fi

  sleep 10
done

if [[ ! -f "$jobs_file" ]]; then
  st="$(get_state "$prep_jobid" || true)"
  echo "[submit][ERROR] jobs.tsv did not appear at $jobs_file (prep state=${st:-unknown})." >&2
  exit 1
fi

n_jobs=$(wc -l < "$jobs_file" | tr -d ' ')
if [[ "$n_jobs" -lt 1 ]]; then
  echo "[submit][ERROR] jobs.tsv exists but has 0 lines: $jobs_file" >&2
  exit 1
fi
echo "[submit] jobs: $n_jobs"

# --- 3) Submit mapping arrays in chunks (MAX_ARRAY_SIZE) ---
array_jobids=()
chunk=0
offset=0

map_base_args=("${sbatch_common[@]}" "${sbatch_export[@]}" --dependency=afterok:"$prep_jobid")
if [[ -n "${PARTITION_MAP:-}" ]]; then map_base_args+=(--partition="${PARTITION_MAP}"); fi
if [[ -n "${MAP_CPUS_PER_TASK:-}" ]]; then map_base_args+=(--cpus-per-task="${MAP_CPUS_PER_TASK}"); fi
if [[ -n "${MAP_MEM_PER_CPU:-}" ]]; then map_base_args+=(--mem-per-cpu="${MAP_MEM_PER_CPU}"); fi

while (( offset < n_jobs )); do
  remaining=$((n_jobs - offset))
  chunk_size=$(( remaining < MAX_ARRAY_SIZE ? remaining : MAX_ARRAY_SIZE ))

  array_spec="0-$((chunk_size-1))%${max_concurrent}"

  jid=$(sbatch --parsable "${map_base_args[@]}"     --array="$array_spec"     "$MAP_SLURM"     "$SCRATCH" "$offset"
  )

  array_jobids+=("$jid")
  echo "[submit] map array chunk ${chunk}: job=${jid} offset=${offset} size=${chunk_size}"
  offset=$((offset + chunk_size))
  chunk=$((chunk + 1))
done

# --- 4) Gather depends on ALL array chunks ---
dep="$(IFS=:; echo "${array_jobids[*]}")"

gather_args=("${sbatch_common[@]}" "${sbatch_export[@]}" --dependency=afterok:"$dep")
if [[ -n "${PARTITION_GATHER:-}" ]]; then gather_args+=(--partition="${PARTITION_GATHER}"); fi

gather_jobid=$(sbatch --parsable "${gather_args[@]}"   "$GATHER_SLURM"   "$SCRATCH" "$RESULTS" "$plasmidSeqData"
)
echo "[submit] gather job: $gather_jobid"

# --- Tracking output ---
job_list="${prep_jobid},$(IFS=,; echo "${array_jobids[*]}"),${gather_jobid}"
echo "[submit] done."
echo "[submit] Track with: squeue -j ${job_list}"
