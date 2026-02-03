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

SCRATCH="$1"
RESULTS="$2"
plasmidSeqData="$3"

echo "[gather] started: $(date)"
echo "[gather] SCRATCH=$SCRATCH"
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
find . -mindepth 1 -maxdepth 1 -type d ! \( -name "Fasta_Reference_Files" -o -name "Stats" -o -name "Reports" \) \
  -exec mv -t "$RESULTS" {} +

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
chown -R "$USER:llshared" "$AlignedData" "$RESULTS"
chmod -R 777 "$AlignedData" "$RESULTS"

echo "[gather] finished: $(date)"
