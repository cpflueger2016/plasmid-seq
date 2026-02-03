#!/usr/bin/env bash
set -euo pipefail

# Defaults
DEFAULT_TSV="/group/llshared/PlasmidSeq/PL_to_fasta.tsv"
DEFAULT_REFS="/group/llshared/PlasmidSeq/Fasta_Reference_Files"
MAX_CONCURRENT=50

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

Example:
  $(basename "$0") -d /home/.../fastqs
  $(basename "$0") -d /home/.../fastqs -t ./PL_to_fasta.tsv -f ./Fasta_Reference_Files -p 80
EOF
}

plasmidSeqData=""
tsv="$DEFAULT_TSV"
refs="$DEFAULT_REFS"

while getopts ":d:t:f:p:h" opt; do
  case "$opt" in
    d) plasmidSeqData="$OPTARG" ;;
    t) tsv="$OPTARG" ;;
    f) refs="$OPTARG" ;;
    p) MAX_CONCURRENT="$OPTARG" ;;
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

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
jobdate="$(date +%Y-%m-%d)"
jobname="plasmidSeq_${jobdate}"

echo "[submit] jobdate=$jobdate"
echo "[submit] plasmidSeqData=$plasmidSeqData"
echo "[submit] TSV=$tsv"
echo "[submit] REFS=$refs"
echo "[submit] max_concurrent=$MAX_CONCURRENT"

# 1) Submit prep job (writes jobs.tsv into scratch)
prep_jobid=$(sbatch --parsable \
  "${script_dir}/plasmidseq_prepare_SLURM.sh" \
  -d "$plasmidSeqData" -t "$tsv" -f "$refs" -j "$jobdate"
)
echo "[submit] prep job: $prep_jobid"

# We can compute the scratch path deterministically (since MYSCRATCH is stable and prep uses jobid)
SCRATCH="${MYSCRATCH}/${jobname}/${prep_jobid}"
RESULTS="/group/llshared/PlasmidSeq/Results/${jobname}/${prep_jobid}"

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
