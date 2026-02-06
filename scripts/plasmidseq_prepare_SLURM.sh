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

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cfg="${PLASMIDSEQ_CONFIG:-${script_dir}/plasmidseq.config}"
if [[ -f "$cfg" ]]; then
  # shellcheck source=/dev/null
  source "$cfg"
fi

PIPE_SCRIPTS="${PIPELINE_SCRIPTS_DIR:-$script_dir}"
MATCHER="${PIPE_SCRIPTS}/${MATCHER_BASENAME:-match_plasmid_fasta_refs_v2.bash}"
FASTA_CLEANER="${PIPE_SCRIPTS}/plasmidseq_clean_fasta_headers.sh"

if command -v module >/dev/null 2>&1; then
  module load gcc || true
fi

if command -v conda >/dev/null 2>&1; then
  conda activate "${CONDA_ENV:-/group/llshared/shared_conda_envs/plasmidseq}" || true
fi
PATH="${CONDA_ENV:-/group/llshared/shared_conda_envs/plasmidseq}:$PATH"
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

# Normalize to an absolute path so rsync source stays valid after changing cwd.
plasmidSeqData="$(cd "$plasmidSeqData" && pwd -P)"

JOBNAME="plasmidSeq_${jobdate}"
SCRATCH="${MYSCRATCH}/${JOBNAME}/${SLURM_JOBID}"
RESULTS="/group/llshared/PlasmidSeq/Results/${JOBNAME}/${SLURM_JOBID}"

mkdir -p "$SCRATCH" "$RESULTS"
mkdir -p "$SCRATCH/Logs"
# capture stdout/stderr from this point forward into scratch Logs (and still echo to slurm stdout)
exec > >(tee -a "$SCRATCH/Logs/slurm-${SLURM_JOBID}.out") 2>&1
echo "[prep] SCRATCH=$SCRATCH"
echo "[prep] RESULTS=$RESULTS"
echo "[prep] plasmidSeqData=$plasmidSeqData"

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

build_snpeff_db() {
  local db_name="${SNPEFF_DB:-plasmidseq_${SLURM_JOBID}}"
  local data_dir="${SNPEFF_DATA_DIR:-${SCRATCH}/snpEff_data}"
  local cfg_file="${SNPEFF_CONFIG_FILE:-${SCRATCH}/snpEff.config}"
  local db_dir="${data_dir}/${db_name}"
  local genes_gbk="${db_dir}/genes.gbk"
  local build_dir="${SCRATCH}/snpEff_build"

  mkdir -p "$db_dir" "$build_dir"
  : > "$genes_gbk"

  if [[ ! -x "$FASTA_CLEANER" ]]; then
    echo "[prep][ERROR] FASTA cleaner not found/executable: $FASTA_CLEANER" >&2
    return 1
  fi

  if [[ -z "${CONDA_ENV_PLANNOTATE:-}" || ! -d "${CONDA_ENV_PLANNOTATE}" ]]; then
    echo "[prep][ERROR] CONDA_ENV_PLANNOTATE is not set to a valid env: ${CONDA_ENV_PLANNOTATE:-<empty>}" >&2
    return 1
  fi
  if ! command -v conda >/dev/null 2>&1; then
    echo "[prep][ERROR] conda not found in PATH; cannot run pLannotate for snpEff DB build." >&2
    return 1
  fi
  if ! command -v "${SNPEFF_BIN:-snpEff}" >/dev/null 2>&1; then
    echo "[prep][ERROR] snpEff binary not found: ${SNPEFF_BIN:-snpEff}" >&2
    return 1
  fi

  declare -A seen_ref_md5=()
  local locus_i=0
  local n_ann=0

  for sample_dir in */PL*/; do
    local ref
    ref=$(find "$sample_dir" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) ! -name "*_clean.fa" ! -name "*_clean.fasta" | head -n 1 || true)
    [[ -z "$ref" ]] && continue

    local clean_ref
    clean_ref="${ref%.*}_clean.fa"
    "$FASTA_CLEANER" -i "$ref" -o "$clean_ref" > "${sample_dir%/}_fasta_header_cleanup.log" 2>&1

    local md5
    md5=$(md5sum "$clean_ref" | awk '{print $1}')
    if [[ -n "${seen_ref_md5[$md5]:-}" ]]; then
      continue
    fi
    seen_ref_md5[$md5]=1

    local sample_name
    sample_name="$(basename "${sample_dir%/}")"
    local pldir="${build_dir}/${sample_name}_plannotate"
    mkdir -p "$pldir"
    if ! conda run -p "$CONDA_ENV_PLANNOTATE" plannotate batch \
      -i "$clean_ref" \
      -o "$pldir" \
      -h -f "$sample_name" > "${pldir}/plannotate_build.log" 2>&1; then
      echo "[prep][WARN] pLannotate failed for $clean_ref; skipping this reference for snpEff DB."
      continue
    fi

    local gbk
    gbk=$(find "$pldir" -maxdepth 1 -type f -name "*.gbk" | head -n 1 || true)
    if [[ -z "$gbk" ]]; then
      echo "[prep][WARN] no gbk produced by pLannotate for $clean_ref; skipping."
      continue
    fi

    locus_i=$((locus_i + 1))
    local locus_id
    locus_id=$(grep -m 1 '^>' "$clean_ref" | sed 's/^>//' | awk '{print $1}')
    if [[ -z "$locus_id" ]]; then
      locus_id=$(printf "ref%06d" "$locus_i")
      echo "[prep][WARN] could not read contig name from $clean_ref; using ${locus_id} for snpEff locus" >&2
    fi
    awk -v locus="$locus_id" '
      NR == 1 && $1 == "LOCUS" { $2 = locus; print; next }
      { print }
    ' "$gbk" >> "$genes_gbk"
    echo >> "$genes_gbk"
    n_ann=$((n_ann + 1))
  done

  if [[ "$n_ann" -eq 0 ]]; then
    echo "[prep][WARN] no pLannotate GenBank records built; skipping snpEff DB build."
    return 0
  fi

  cat > "$cfg_file" <<EOF
${db_name}.genome : plasmidSeq dynamic run DB
EOF

  echo "[prep] building snpEff DB: ${db_name}"
  local snpeff_log="${SCRATCH}/Logs/snpeff_build.log"
  local snpeff_env=()
  if [[ -n "${SNPEFF_JAVA_HOME:-}" && -x "${SNPEFF_JAVA_HOME}/bin/java" ]]; then
    echo "[prep] using snpEff JAVA_HOME=${SNPEFF_JAVA_HOME}"
    snpeff_env=(env "JAVA_HOME=${SNPEFF_JAVA_HOME}" "PATH=${SNPEFF_JAVA_HOME}/bin:${PATH}")
  fi
  if ! "${snpeff_env[@]}" "${SNPEFF_BIN:-snpEff}" build -genbank -v -noCheckProtein -c "$cfg_file" -dataDir "$data_dir" "$db_name" \
    > "$snpeff_log" 2>&1; then
    echo "[prep][ERROR] snpEff DB build command failed for ${db_name}. See $snpeff_log" >&2
    return 1
  fi
  if [[ ! -s "${db_dir}/snpEffectPredictor.bin" ]]; then
    echo "[prep][ERROR] snpEff DB build did not produce snpEffectPredictor.bin for ${db_name}. See $snpeff_log" >&2
    return 1
  fi
  echo "[prep] snpEff DB ready: db=${db_name} cfg=${cfg_file} dataDir=${data_dir}"
}

if [[ "${ENABLE_SNPEFF:-0}" == "1" ]]; then
  if ! build_snpeff_db; then
    echo "[prep][ERROR] snpEff pre-map DB build failed." >&2
    exit 1
  fi
fi

# Create jobs file
: > jobs.tsv
shopt -s nullglob
for i in */PL*/; do
  # Exclude helper folders created during snpEff DB build.
  [[ "$i" == snpEff_build/* ]] && continue

  one_matches=( ${i}*_R1_*gz )
  two_matches=( ${i}*_R2_*gz )
  if [[ ${#one_matches[@]} -eq 0 || ${#two_matches[@]} -eq 0 ]]; then
    echo "[prep][WARN] skipping folder without paired FASTQs: ${i%/}"
    continue
  fi

  ref=$(find "$(pwd)/${i}" -maxdepth 1 -type f \( -name "*_clean.fa" -o -name "*_clean.fasta" \) | head -n 1)
  if [[ -z "$ref" ]]; then
    ref=$(find "$(pwd)/${i}" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" -o -name "na" \) | head -n 1)
  fi
  ref=${ref##*/}

  one=${one_matches[0]##*/}
  two=${two_matches[0]##*/}

  uID=$(echo "${i}" | perl -ne '/(PL\d{4,})/ && print $1."\n"')
  printf "%s\t%s\t%s\t%s\t%s\n" "${i%/}" "$ref" "$one" "$two" "$uID" >> jobs.tsv
done
shopt -u nullglob

echo "[prep] jobs: $(wc -l < jobs.tsv | tr -d ' ')"
