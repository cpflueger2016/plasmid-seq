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

# Normalize pLannotate GBK headers in RESULTS using the snpEff build copies.
if [[ -d "$RESULTS/snpEff_build" ]]; then
  while IFS= read -r gbk; do
    sample_dir="$(basename "$(dirname "$gbk")")"
    sample_name="${sample_dir%_plannotate}"
    target_gbk="$RESULTS/3-20kb_plasmids/${sample_name}/${sample_name}_plannotate/$(basename "$gbk")"
    if [[ -f "$target_gbk" ]]; then
      cp -f "$gbk" "$target_gbk"
      echo "[gather] updated pLannotate GBK header: $target_gbk"
    fi
  done < <(find "$RESULTS/snpEff_build" -type f -name "*.gbk")
fi

# Optionally build run summary (CSV + HTML) before copying to Aligned.
SUMMARY_SCRIPT=""
SAMPLE_REPORT_SCRIPT=""
if [[ -n "${PIPELINE_SCRIPTS_DIR:-}" && -f "${PIPELINE_SCRIPTS_DIR}/plasmidseq_run_summary.py" ]]; then
  SUMMARY_SCRIPT="${PIPELINE_SCRIPTS_DIR}/plasmidseq_run_summary.py"
elif [[ -f "${script_dir}/plasmidseq_run_summary.py" ]]; then
  SUMMARY_SCRIPT="${script_dir}/plasmidseq_run_summary.py"
elif [[ -f "$(cd "$(dirname "$cfg")" && pwd -P)/plasmidseq_run_summary.py" ]]; then
  SUMMARY_SCRIPT="$(cd "$(dirname "$cfg")" && pwd -P)/plasmidseq_run_summary.py"
fi
if [[ -n "${PIPELINE_SCRIPTS_DIR:-}" && -f "${PIPELINE_SCRIPTS_DIR}/plasmidseq_sample_report.py" ]]; then
  SAMPLE_REPORT_SCRIPT="${PIPELINE_SCRIPTS_DIR}/plasmidseq_sample_report.py"
elif [[ -f "${script_dir}/plasmidseq_sample_report.py" ]]; then
  SAMPLE_REPORT_SCRIPT="${script_dir}/plasmidseq_sample_report.py"
elif [[ -f "$(cd "$(dirname "$cfg")" && pwd -P)/plasmidseq_sample_report.py" ]]; then
  SAMPLE_REPORT_SCRIPT="$(cd "$(dirname "$cfg")" && pwd -P)/plasmidseq_sample_report.py"
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

# Build per-sample report (JSON + HTML) with coverage comparison.
if [[ -z "$SAMPLE_REPORT_SCRIPT" || ! -f "$SAMPLE_REPORT_SCRIPT" ]]; then
  echo "[gather][WARN] sample report script not found; skipping per-sample reports: $SAMPLE_REPORT_SCRIPT"
elif ! command -v python3 >/dev/null 2>&1 && ! command -v python >/dev/null 2>&1; then
  echo "[gather][WARN] no python interpreter found; skipping per-sample reports."
else
  echo "[gather] building per-sample reports with: $SAMPLE_REPORT_SCRIPT"
  [[ -x "$SAMPLE_REPORT_SCRIPT" ]] || chmod +x "$SAMPLE_REPORT_SCRIPT" || true
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN_SAMPLE="python3"
  else
    PYTHON_BIN_SAMPLE="python"
  fi

  report_count=0
  while IFS= read -r fastp_json; do
    sample_dir="$(dirname "$fastp_json")"
    sample_name="$(basename "$sample_dir")"
    report_prefix="${sample_dir}/${sample_name}_sample_report"
    echo "[gather] sample report: $sample_dir"
    if ! "$PYTHON_BIN_SAMPLE" "$SAMPLE_REPORT_SCRIPT" -s "$sample_dir" -o "$report_prefix" \
      >> "$SCRATCH/Logs/sample_reports.log" 2>&1; then
      echo "[gather][WARN] sample report failed for $sample_dir (continuing)"
      tail -n 20 "$SCRATCH/Logs/sample_reports.log" || true
    else
      report_count=$((report_count + 1))
    fi
  done < <(find "$RESULTS" -type f -name "*_fastp_report.json" | sort)
  echo "[gather] per-sample reports generated: $report_count"
fi

# Organize per-sample output folders.
organize_sample_dir() {
  local sample_dir="$1"
  [[ -d "$sample_dir" ]] || return 0

  local logs_dir="$sample_dir/Logs"
  local summary_dir="$sample_dir/Summary"
  local snp_dir="$sample_dir/SNP_INDEL"
  local fastq_dir="$sample_dir/FASTQ"
  local bam_dir="$sample_dir/BAM"

  mkdir -p "$logs_dir" "$summary_dir" "$snp_dir" "$fastq_dir" "$bam_dir"

  # Logs and small status files
  mv -f "$sample_dir"/*_log "$logs_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_markdup.log "$logs_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_bbmap.log "$logs_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_varscan.log "$logs_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_fasta_header_cleanup.log "$logs_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_unicycler_vs_ref.log "$logs_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/FASTA_REF_* "$logs_dir"/ 2>/dev/null || true

  # Summary reports
  mv -f "$sample_dir"/*_fastp_report.html "$summary_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_fastp_report.json "$summary_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_sample_report.html "$summary_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_sample_report.json "$summary_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_sample_report_coverage_tracks.tsv "$summary_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/snpEff_summary.html "$summary_dir"/ 2>/dev/null || true

  # SNP / INDEL outputs
  mv -f "$sample_dir"/*_varscan_*.vcf "$snp_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_varscan_summary.tsv "$snp_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/snpEff_genes.txt "$snp_dir"/ 2>/dev/null || true

  # FASTQ (raw)
  mv -f "$sample_dir"/*_R1_001.fastq.gz "$fastq_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*_R2_001.fastq.gz "$fastq_dir"/ 2>/dev/null || true

  # Remove trimmed FASTQ files
  rm -f "$sample_dir"/*_R1_trimmed.fastq "$sample_dir"/*_R2_trimmed.fastq 2>/dev/null || true
  rm -f "$sample_dir"/*_R1_trimmed.fastq.gz "$sample_dir"/*_R2_trimmed.fastq.gz 2>/dev/null || true

  # BAM / BAI / mpileup (compress mpileup)
  if compgen -G "$sample_dir/*.mpileup" >/dev/null; then
    gzip -f "$sample_dir"/*.mpileup 2>/dev/null || true
  fi
  mv -f "$sample_dir"/*.bam "$bam_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*.bam.bai "$bam_dir"/ 2>/dev/null || true
  mv -f "$sample_dir"/*.mpileup.gz "$bam_dir"/ 2>/dev/null || true

  # Keep unicycler assembly dir in place (no move)
}

while IFS= read -r sample_dir; do
  organize_sample_dir "$sample_dir"
done < <(find "$RESULTS" -mindepth 2 -maxdepth 2 -type d -path "*/3-20kb_plasmids/*")

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
