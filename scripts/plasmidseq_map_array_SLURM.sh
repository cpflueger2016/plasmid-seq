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
#   1) SCRATCH dir (where jobs.tsv lives and sample folders are)
SCRATCH="$1"
JOBS="$SCRATCH/jobs.tsv"

module load gcc
conda activate /group/llshared/shared_conda_envs/plasmidseq
PATH=/group/llshared/shared_conda_envs/plasmidseq:$PATH
export PATH

mkdir -p "$SCRATCH/Logs"

# SLURM_ARRAY_TASK_ID is 0-based in our submission; sed is 1-based:
line_num=$((SLURM_ARRAY_TASK_ID + 1))

# columns: folder   ref   R1   R2   uID
IFS=$'\t' read -r folder ref r1 r2 uid < <(sed -n "${line_num}p" "$JOBS")

echo "[map] task=$SLURM_ARRAY_TASK_ID line=$line_num uid=$uid folder=$folder ref=$ref"
cd "$SCRATCH/$folder"

bash /group/llshared/PlasmidSeq/plasmidseq_mapper_PE.sh \
  -1 "$r1" -2 "$r2" -r "$ref" \
  -c -m 300 -u "$uid" -y -s -q 30
