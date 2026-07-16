#!/usr/bin/env bash
# Run plasmid-seq locally, without Slurm. Requires the environment in ../environment.yml.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
matcher="${script_dir}/match_plasmid_fasta_refs_v2.bash"
mapper="${script_dir}/plasmidseq_mapper_PE.sh"
sample_report="${script_dir}/plasmidseq_sample_report.py"
run_summary="${script_dir}/plasmidseq_run_summary.py"

usage() {
  cat <<'EOF'
Usage: plasmidseq_run_local.sh -d FASTQ_DIR -t PL_TO_FASTA.tsv -f REFS_DIR [options]

Required:
  -d DIR     input FASTQ directory (nested project folders are supported)
  -t FILE    PL-to-reference TSV
  -f DIR     reference FASTA directory

Options:
  -o DIR     work/output root (default: ./plasmidseq-local-runs)
  -j INT     concurrent samples (default: 1)
  -J INT     threads per sample (default: 4)
  -w FILE    optional plate map CSV for run summary
  -V         enable VarScan (requires varscan in the active environment)
  -P DIR     optional pLannotate Conda environment prefix
  -h         show this help

Results are written to <output-root>/<run-id>/results. Staged inputs and logs are
retained under <output-root>/<run-id>/work for reproducibility and troubleshooting.
EOF
}

fastq_dir=""; tsv=""; refs=""; output_root="$(pwd -P)/plasmidseq-local-runs"
concurrency=1; threads=4; plate_map=""; enable_variants=0; plannotate_env=""
while getopts ':d:t:f:o:j:J:w:VP:h' opt; do
  case "$opt" in
    d) fastq_dir="$OPTARG" ;; t) tsv="$OPTARG" ;; f) refs="$OPTARG" ;;
    o) output_root="$OPTARG" ;; j) concurrency="$OPTARG" ;; J) threads="$OPTARG" ;;
    w) plate_map="$OPTARG" ;; V) enable_variants=1 ;; P) plannotate_env="$OPTARG" ;;
    h) usage; exit 0 ;; *) usage >&2; exit 2 ;;
  esac
done

[[ -d "$fastq_dir" ]] || { echo "[local][ERROR] FASTQ directory not found: $fastq_dir" >&2; exit 2; }
[[ -f "$tsv" ]] || { echo "[local][ERROR] TSV not found: $tsv" >&2; exit 2; }
[[ -d "$refs" ]] || { echo "[local][ERROR] references directory not found: $refs" >&2; exit 2; }
[[ -z "$plate_map" || -f "$plate_map" ]] || { echo "[local][ERROR] plate map not found: $plate_map" >&2; exit 2; }
[[ -z "$plannotate_env" || -d "$plannotate_env" ]] || { echo "[local][ERROR] pLannotate environment not found: $plannotate_env" >&2; exit 2; }
[[ "$concurrency" =~ ^[1-9][0-9]*$ ]] || { echo "[local][ERROR] -j must be a positive integer" >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "[local][ERROR] -J must be a positive integer" >&2; exit 2; }
[[ -x "$matcher" && -x "$mapper" ]] || { echo "[local][ERROR] pipeline scripts are not executable" >&2; exit 2; }

for tool in rsync perl python; do
  command -v "$tool" >/dev/null || { echo "[local][ERROR] missing required command: $tool" >&2; exit 2; }
done

export ENABLE_VARIANTS="$enable_variants" ENABLE_SNPEFF=0 THREADS="$threads"
export CONDA_ENV_PLANNOTATE="$plannotate_env"

run_id="plasmidSeq_local_$(date +%Y%m%d_%H%M%S)"
run_root="$(mkdir -p "$output_root" && cd "$output_root" && pwd -P)/$run_id"
scratch="$run_root/work"
results="$run_root/results"
logs="$run_root/logs"
mkdir -p "$scratch" "$results" "$logs"

echo "[local] run_id=$run_id"
echo "[local] work=$scratch"
echo "[local] results=$results"
echo "[local] concurrency=$concurrency threads_per_sample=$threads"

cp -f "$tsv" "$scratch/PL_to_fasta.tsv"
cp -R "$refs" "$scratch/Fasta_Reference_Files"
fastq_dir="$(cd "$fastq_dir" && pwd -P)"
rsync -a --prune-empty-dirs \
  --exclude='fastqc' --exclude='multiqc' --exclude='Reports' --exclude='Stats' --exclude='Logs' \
  --exclude='Undetermined*' --include='*/' --include='*.fastq.gz' --exclude='*' \
  "$fastq_dir/" "$scratch/"

n_fastq="$(find "$scratch" -type f -name '*.fastq.gz' | wc -l | tr -d ' ')"
[[ "$n_fastq" -gt 0 ]] || { echo "[local][ERROR] no FASTQs staged" >&2; exit 1; }

# Match the historical staging convention: project/PL... directories hold a paired FASTQ.
shopt -s nullglob
for r1 in "$scratch"/*/PL*_R1_*.fastq.gz; do
  r2="${r1/_R1_/_R2_}"
  [[ -f "$r2" ]] || { echo "[local][WARN] missing R2 for $r1" >&2; continue; }
  sample_dir="${r1%_S*}"
  mkdir -p "$sample_dir"
  mv "$r1" "$r2" "$sample_dir/"
done
shopt -u nullglob

( cd "$scratch" && "$matcher" -r "$scratch" -l "$scratch/plasmid_fasta_match.log" "$scratch/PL_to_fasta.tsv" "$scratch/Fasta_Reference_Files" )

jobs="$scratch/jobs.tsv"
: > "$jobs"
while IFS= read -r sample_dir; do
  r1=( "$sample_dir"/*_R1_*.fastq.gz )
  r2=( "$sample_dir"/*_R2_*.fastq.gz )
  [[ ${#r1[@]} -gt 0 && ${#r2[@]} -gt 0 ]] || continue
  ref="$(find "$sample_dir" -maxdepth 1 -type f \( -name '*_clean.fa' -o -name '*_clean.fasta' -o -name '*.fa' -o -name '*.fasta' \) | head -n 1 || true)"
  uid="$(basename "$sample_dir" | perl -ne '/(PL\d{4,})/ && print $1')"
  [[ -n "$ref" && -n "$uid" ]] || { echo "[local][WARN] skipping unmatched sample: $sample_dir" >&2; continue; }
  printf '%s\t%s\t%s\t%s\t%s\n' "${sample_dir#"$scratch/"}" "${ref##*/}" "${r1[0]##*/}" "${r2[0]##*/}" "$uid" >> "$jobs"
done < <(find "$scratch" -mindepth 2 -maxdepth 2 -type d -name 'PL*' | sort)

n_jobs="$(wc -l < "$jobs" | tr -d ' ')"
[[ "$n_jobs" -gt 0 ]] || { echo "[local][ERROR] no matched paired samples found" >&2; exit 1; }
echo "[local] samples=$n_jobs"

run_one() {
  local folder="$1" ref="$2" r1="$3" r2="$4" uid="$5"
  ( cd "$scratch/$folder" && THREADS="$threads" bash "$mapper" -1 "$r1" -2 "$r2" -r "$ref" -c -m 300 -u "$uid" -y -s -q 30 ) \
    >"$logs/${uid}.log" 2>&1
}

declare -a pids=()
declare -a pid_uids=()
failures=0
reap_one() {
  local pid="${pids[0]}" uid="${pid_uids[0]}"
  if ! wait "$pid"; then
    echo "[local][ERROR] sample failed: $uid (log: $logs/${uid}.log)" >&2
    failures=$((failures + 1))
  fi
  pids=( "${pids[@]:1}" ); pid_uids=( "${pid_uids[@]:1}" )
}
while IFS=$'\t' read -r folder ref r1 r2 uid; do
  run_one "$folder" "$ref" "$r1" "$r2" "$uid" &
  pids+=( "$!" ); pid_uids+=( "$uid" )
  if [[ ${#pids[@]} -ge $concurrency ]]; then reap_one; fi
done < "$jobs"
while [[ ${#pids[@]} -gt 0 ]]; do reap_one; done

cp -p "$scratch/jobs.tsv" "$scratch/plasmid_fasta_match.log" "$results/" 2>/dev/null || true
if [[ -n "$plate_map" ]]; then
  cp -p "$plate_map" "$results/PL_to_plate_position.csv"
fi
find "$scratch" -mindepth 2 -maxdepth 2 -type d -name 'PL*' -exec cp -R {} "$results/" \;

while IFS= read -r report_json; do
  sample_dir="$(dirname "$report_json")"
  sample_name="$(basename "$sample_dir")"
  python "$sample_report" -s "$sample_dir" -o "$sample_dir/${sample_name}_sample_report" \
    >>"$logs/sample_reports.log" 2>&1 || echo "[local][WARN] report failed: $sample_name" >&2
done < <(find "$results" -type f -name '*_fastp_report.json' | sort)
if [[ -n "$plate_map" ]]; then
  python "$run_summary" -r "$results" -m "$plate_map" >"$logs/run_summary.log" 2>&1 || echo "[local][WARN] run summary failed" >&2
fi
cp -R "$logs" "$results/Logs"

if [[ "$failures" -gt 0 ]]; then
  echo "[local][ERROR] completed with $failures failed sample(s); staged work was retained at $scratch" >&2
  exit 1
fi
echo "[local] completed successfully: $results"
