#!/bin/bash --login
#SBATCH --job-name=plasmidSeq_map
#SBATCH --partition=work
#SBATCH --mem-per-cpu=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --export=ALL
#SBATCH --output=Logs/%x_%A_%a.out
#SBATCH --error=Logs/%x_%A_%a.err
#SBATCH --mail-type=FAIL

set -euo pipefail

# args:
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cfg="${PLASMIDSEQ_CONFIG:-${script_dir}/plasmidseq.config}"
if [[ -f "$cfg" ]]; then
  # shellcheck source=/dev/null
  source "$cfg"
fi


#   1) SCRATCH dir (where jobs.tsv lives and sample folders are)
SCRATCH="$1"
JOBS="$SCRATCH/jobs.tsv"
OFFSET="${2:-0}"   # offset into jobs.tsv (0-based)

# --- Resolve mapper script ---
resolve_script() {
  local explicit="$1" localname="$2" fallback="$3" base_dir="$4"
  if [[ -n "$explicit" && -f "$explicit" ]]; then echo "$explicit"; return 0; fi
  if [[ -n "$base_dir" && -f "$base_dir/$localname" ]]; then echo "$base_dir/$localname"; return 0; fi
  if [[ -f "$script_dir/$localname" ]]; then echo "$script_dir/$localname"; return 0; fi
  if [[ -n "$fallback" && -f "$fallback" ]]; then echo "$fallback"; return 0; fi
  return 1
}
PIPE_SCRIPTS="${PIPELINE_SCRIPTS_DIR:-}"
MAPPER="$(resolve_script "${MAPPER_PATH:-}" "${MAPPER_BASENAME:-plasmidseq_mapper_PE.sh}" "${MAPPER_FALLBACK:-}" "$PIPE_SCRIPTS" || true)"
if [[ -z "${MAPPER:-}" ]]; then
  echo "[map][ERROR] Cannot locate mapper script (set MAPPER_PATH or PIPELINE_SCRIPTS_DIR in config)." >&2
  echo "[map][ERROR] cfg=${cfg}" >&2
  echo "[map][ERROR] script_dir=${script_dir}" >&2
  echo "[map][ERROR] PIPELINE_SCRIPTS_DIR=${PIPE_SCRIPTS:-<empty>}" >&2
  echo "[map][ERROR] MAPPER_BASENAME=${MAPPER_BASENAME:-plasmidseq_mapper_PE.sh}" >&2
  echo "[map][ERROR] Tried: ${MAPPER_PATH:-<empty>}, ${PIPE_SCRIPTS:-<empty>}/${MAPPER_BASENAME:-plasmidseq_mapper_PE.sh}, ${script_dir}/${MAPPER_BASENAME:-plasmidseq_mapper_PE.sh}, ${MAPPER_FALLBACK:-<empty>}" >&2
  exit 1
fi

# Mapper runs in a child bash process, so export config vars it consumes.
export CONDA_ENV CONDA_ENV_PLANNOTATE CONDA_INIT MODULES
export ENABLE_VARIANTS VARSCAN_BIN VARSCAN_JAR VARSCAN_MIN_COVERAGE VARSCAN_MIN_VAR_FREQ VARSCAN_PVALUE
export ENABLE_SNPEFF SNPEFF_BIN SNPEFF_DB
export SNPEFF_CONFIG_FILE SNPEFF_DATA_DIR

# --- Environment ---
if command -v module >/dev/null 2>&1 && [[ -n "${MODULES:-}" ]]; then
  for m in ${MODULES}; do
    module load "$m"
  done
fi

if [[ -n "${CONDA_INIT:-}" && -f "${CONDA_INIT}" ]]; then
  # shellcheck source=/dev/null
  source "${CONDA_INIT}"
fi

if command -v conda >/dev/null 2>&1 && [[ -n "${CONDA_ENV:-}" ]]; then
  conda activate "${CONDA_ENV}"
fi

mkdir -p "$SCRATCH/Logs"

# SLURM_ARRAY_TASK_ID is 0-based in our submission; sed is 1-based:
line_num=$((OFFSET + SLURM_ARRAY_TASK_ID + 1))

# columns: folder   ref   R1   R2   uID
# Use awk field extraction so empty tab fields (e.g. missing ref) are preserved.
job_line=$(awk -F'\t' -v n="$line_num" '
  NR == n { printf "%s\034%s\034%s\034%s\034%s", $1, $2, $3, $4, $5; exit }
' "$JOBS")
IFS=$'\034' read -r folder ref r1 r2 uid <<< "${job_line:-}"

if [[ -z "${folder:-}" || -z "${uid:-}" ]]; then
  echo "[map][ERROR] No job entry found for line ${line_num} (OFFSET=${OFFSET}, task=${SLURM_ARRAY_TASK_ID})." >&2
  exit 1
fi

echo "[map] task=$SLURM_ARRAY_TASK_ID line=$line_num uid=$uid folder=$folder ref=$ref"
cd "$SCRATCH/$folder"

bash "$MAPPER" \
  -1 "$r1" -2 "$r2" -r "$ref" \
  -c -m 300 -u "$uid" -y -s -q 30
