#!/bin/bash --login
#SBATCH --job-name=plasmidSeq_gather
#SBATCH --partition=ll
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cfg="${PLASMIDSEQ_CONFIG:-${script_dir}/plasmidseq.config}"
if [[ -f "$cfg" ]]; then
  # shellcheck source=/dev/null
  source "$cfg"
fi

SCRATCH="$1"
RESULTS="$2"
plasmidSeqData="$3"
PLATE_MAP_CSV="${4:-}"

echo "[gather] started: $(date)"
echo "[gather] SCRATCH=$SCRATCH"
mkdir -p "$SCRATCH/Logs"
exec > >(tee -a "$SCRATCH/Logs/slurm-${SLURM_JOBID}.out") 2>&1
echo "[gather] RESULTS=$RESULTS"

mkdir -p "$RESULTS"

# --- Keep key run-level files
for f in jobs.tsv plasmid_fasta_match.log; do
  if [[ -f "$SCRATCH/$f" ]]; then
    echo "[gather] saving $f -> $RESULTS/"
    cp -f "$SCRATCH/$f" "$RESULTS/"
  else
    echo "[gather][WARN] missing $SCRATCH/$f"
  fi
done



# Move everything except refs/junk into results
cd "$SCRATCH"
find . -mindepth 1 -maxdepth 1 -type d ! \( -name "Fasta_Reference_Files" -o -name "Stats" -o -name "Reports" -o -name "Logs" \) \
  -exec mv -t "$RESULTS" {} +

# Optionally build run summary (CSV + HTML) before copying to Aligned.
SUMMARY_SCRIPT=""
if [[ -n "${PIPELINE_SCRIPTS_DIR:-}" && -f "${PIPELINE_SCRIPTS_DIR}/plasmidseq_run_summary.py" ]]; then
  SUMMARY_SCRIPT="${PIPELINE_SCRIPTS_DIR}/plasmidseq_run_summary.py"
elif [[ -f "${script_dir}/plasmidseq_run_summary.py" ]]; then
  SUMMARY_SCRIPT="${script_dir}/plasmidseq_run_summary.py"
elif [[ -f "$(cd "$(dirname "$cfg")" && pwd -P)/plasmidseq_run_summary.py" ]]; then
  SUMMARY_SCRIPT="$(cd "$(dirname "$cfg")" && pwd -P)/plasmidseq_run_summary.py"
fi
if [[ -n "$PLATE_MAP_CSV" ]]; then
  if [[ ! -f "$PLATE_MAP_CSV" ]]; then
    echo "[gather][WARN] plate map CSV not found; skipping run summary: $PLATE_MAP_CSV"
  elif [[ -z "$SUMMARY_SCRIPT" || ! -f "$SUMMARY_SCRIPT" ]]; then
    echo "[gather][WARN] summary script not found; skipping run summary: $SUMMARY_SCRIPT"
  else
    echo "[gather] building run summary with plate map: $PLATE_MAP_CSV"
    echo "[gather] summary script: $SUMMARY_SCRIPT"
    if command -v python3 >/dev/null 2>&1; then
      PYTHON_BIN="python3"
    elif command -v python >/dev/null 2>&1; then
      PYTHON_BIN="python"
    else
      PYTHON_BIN=""
    fi

    if [[ -z "$PYTHON_BIN" ]]; then
      echo "[gather][WARN] no python interpreter found; skipping run summary."
    elif ! "$PYTHON_BIN" "$SUMMARY_SCRIPT" -r "$RESULTS" -m "$PLATE_MAP_CSV" > "$SCRATCH/Logs/run_summary.log" 2>&1; then
      echo "[gather][WARN] run summary generation failed; continuing gather."
      tail -n 40 "$SCRATCH/Logs/run_summary.log" || true
    else
      echo "[gather] run summary generated: $RESULTS/run_summary.csv and $RESULTS/run_summary.html"
    fi
  fi
else
  echo "[gather] no plate map CSV provided; skipping run summary."
fi

# Preserve all gather/map logs in RESULTS for easier debugging
if [[ -d "$SCRATCH/Logs" ]]; then
  mkdir -p "$RESULTS/Logs"
  cp -af "$SCRATCH/Logs/." "$RESULTS/Logs/" || true
fi

# Copy results also to Run/Aligned (your original logic)
AlignedData=${plasmidSeqData/fastq*/}"/Aligned"
mkdir -p "$AlignedData"
cp -fr "$RESULTS"/* "$AlignedData"


# --- Cleanup scratch (be paranoid)
echo "[gather] cleanup: removing scratch dir $SCRATCH"

# refuse to delete anything that looks suspicious
if [[ -z "${SCRATCH:-}" || "$SCRATCH" == "/" || "$SCRATCH" == "." ]]; then
  echo "[gather][ERROR] SCRATCH path is unsafe ('$SCRATCH'); refusing to delete." >&2
  exit 1
fi

# Strong safety check: must live under MYSCRATCH if available
if [[ -n "${MYSCRATCH:-}" ]]; then
  case "$SCRATCH" in
    "$MYSCRATCH"/*) ;;
    *)
      echo "[gather][ERROR] SCRATCH ('$SCRATCH') is not under MYSCRATCH ('$MYSCRATCH'); refusing to delete." >&2
      exit 1
      ;;
  esac
fi

rm -rf "$SCRATCH"
echo "[gather] scratch removed."


# Permissions
USER=${USER:-$(whoami)}
chown -R "$USER:${RESULTS_GROUP:-llusers}" "$AlignedData" "$RESULTS" || true
chmod -R 777 "$AlignedData" "$RESULTS"

echo "[gather] finished: $(date)"
