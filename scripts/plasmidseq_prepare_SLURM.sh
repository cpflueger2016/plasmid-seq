#!/bin/bash --login
#SBATCH --job-name=plasmidSeq_prep
#SBATCH --partition=ll
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL

set -euo pipefail

SCRIPT_DIR="/mmfs1/data/group/llshared/PlasmidSeq/plasmid-seq/scripts"
MATCHER="${SCRIPT_DIR}/match_plasmid_fasta_refs_v2.bash"


module load gcc
conda activate /group/llshared/shared_conda_envs/plasmidseq
PATH=/group/llshared/shared_conda_envs/plasmidseq:$PATH
export PATH

usage() {
  echo "Usage: $0 -d <plasmidSeqData_dir> [-t PL_to_fasta.tsv] [-f Fasta_Reference_Files_dir] [-j YYYY-MM-DD]" >&2
  exit 2
}

plasmidSeqData=""
tsv="/group/llshared/PlasmidSeq/PL_to_fasta.tsv"
refs="/group/llshared/PlasmidSeq/Fasta_Reference_Files"
jobdate="$(date +%Y-%m-%d)"

while getopts ":d:t:f:j:" opt; do
  case "$opt" in
    d) plasmidSeqData="$OPTARG" ;;
    t) tsv="$OPTARG" ;;
    f) refs="$OPTARG" ;;
    j) jobdate="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -z "$plasmidSeqData" ]] && usage
[[ ! -d "$plasmidSeqData" ]] && { echo "[prep][ERROR] plasmidSeqData dir not found: $plasmidSeqData" >&2; exit 1; }
[[ ! -f "$tsv" ]] && { echo "[prep][ERROR] TSV not found: $tsv" >&2; exit 1; }
[[ ! -d "$refs" ]] && { echo "[prep][ERROR] Refs dir not found: $refs" >&2; exit 1; }

JOBNAME="plasmidSeq_${jobdate}"
SCRATCH="${MYSCRATCH}/${JOBNAME}/${SLURM_JOBID}"
RESULTS="/group/llshared/PlasmidSeq/Results/${JOBNAME}/${SLURM_JOBID}"

mkdir -p "$SCRATCH" "$RESULTS"
echo "[prep] SCRATCH=$SCRATCH"
echo "[prep] RESULTS=$RESULTS"

# stage inputs into scratch
cp -f "$tsv" "$SCRATCH/PL_to_fasta.tsv"
cp -r "$refs" "$SCRATCH/Fasta_Reference_Files"

cd "$SCRATCH"

# Stage fastqs from plasmidSeqData (NO hard-coded path)
rsync -av \
  --exclude='fastqc' --exclude='multiqc' --exclude='Reports' --exclude='Stats' --exclude='Logs' \
  --exclude='Undetermined*' --prune-empty-dirs \
  --include='*/' \
  --include='*.fastq.gz' \
  --exclude='*' \
  "${plasmidSeqData}/" ./

# Check if fastqs have been copied
n_fastq=$(find . -type f -name '*.fastq.gz' | wc -l | tr -d ' ')
if [[ "$n_fastq" -eq 0 ]]; then
  echo "[prep][ERROR] No *.fastq.gz files were staged from: ${plasmidSeqData}" >&2
  echo "[prep][ERROR] Check the source path and whether files are nested deeper than expected." >&2
  exit 1
fi
echo "[prep] staged fastqs: $n_fastq"


# Move Fastq files to folders
shopt -s nullglob
for i in */PL*R1*; do
  f=${i%_S*}
  mkdir -p "$f"
  n=${i/_R1_/_R2_}
  mv "$i" "$n" "$f"
done
shopt -u nullglob

# Match plasmid fasta to samples
chmod +x "$MATCHER"
"$MATCHER" \
  -r "$SCRATCH" \
  -l "$SCRATCH/plasmid_fasta_match.log" \
  ${VERBOSE_FLAG:-} \
  "$tsv" "$refs"

# Create jobs file
: > jobs.tsv
for i in */PL*/; do
  ref=$(find "$(pwd)/${i}" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" -o -name "na" \) | head -n 1)
  ref=${ref##*/}

  one=$(ls ${i}*_R1_*gz); one=${one##*/}
  two=$(ls ${i}*_R2_*gz); two=${two##*/}

  uID=$(echo "${i}" | perl -ne '/(PL\d{4,})/ && print $1."\n"')
  printf "%s\t%s\t%s\t%s\t%s\n" "${i%/}" "$ref" "$one" "$two" "$uID" >> jobs.tsv
done

echo "[prep] jobs: $(wc -l < jobs.tsv | tr -d ' ')"
