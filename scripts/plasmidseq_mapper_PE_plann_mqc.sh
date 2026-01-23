#!/bin/bash

### Maps a plasmid against bowtie2 reference as well as de novo assembles fastq file
### requires: fastp, bowtie2, sambamba, SPAdes, samtools, seqkit, pigz, bcftools, tabix, unicycler, pilon, java, plannotate, multiqc, python
### version 0.8.2 + pLannotate + MultiQC

version="0.8.2"
skipConsenus="TRUE"
assembleUnmappedReads="FALSE"
skipSpadesDenovoAssembly=0
qval=20
cleanSPADES=1
minLengthSPADES=0
uniqueID=""
runUnicycler=0
runPLAnnotate=1
runMultiQC=1
plannotateParser="/group/llshared/PlasmidSeq/parse_plannotate_for_multiqc.py"

### Check dependencies

if ! [ -x "$(command -v fastp)" ]; then
  echo 'Error: fastp is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v bowtie2-align-s)" ]; then
  echo 'Error: bowtie2 is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v sambamba)" ]; then
  echo 'Error: sambamba is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v pigz)" ]; then
  echo 'Error: pigz is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v bcftools)" ]; then
  echo 'Error: bcftools is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v tabix)" ]; then
  echo 'Error: tabix is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v plasmidspades.py)" ]; then
  echo 'Error: SPAdes is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v seqkit)" ]; then
  echo 'Error: seqkit is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v unicycler)" ]; then
  echo 'Error: unicycler is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v pilon)" ]; then
  echo 'Error: pilon is not installed.' >&2
  exit 1
fi
if ! command -v plannotate >/dev/null 2>&1; then
  echo 'Error: plannotate is not installed.' >&2
  exit 1
fi
if ! command -v multiqc >/dev/null 2>&1; then
  echo 'Error: multiqc not installed.' >&2
  exit 1
fi
if ! command -v python >/dev/null 2>&1; then
  echo 'Error: python not installed.' >&2
  exit 1
fi
if [ ! -f "${plannotateParser}" ]; then
  echo "Error: pLannotate parser script missing: ${plannotateParser}" >&2
  exit 1
fi
########################################
# Fix pLannotate bugs dynamically (no root required)
########################################

echo "Applying pLannotate runtime patch ..."

# Where to put the patched copy
PATCH_DIR="/tmp/plannotate_${USER}"
export PATCH_DIR

# Fresh start each run
rm -rf "$PATCH_DIR"
mkdir -p "$PATCH_DIR"

# Locate original plannotate package
ORIG_PKG="$(python - << 'EOF'
import plannotate, os
print(os.path.dirname(plannotate.__file__))
EOF
)"

# Copy original package into the patch dir
cp -r "$ORIG_PKG" "$PATCH_DIR/"

########################################
# PATCH 1 — fix any(1) → any(axis=1) in annotate.py
########################################

PATCH_FILE_ANNOTATE="${PATCH_DIR}/plannotate/annotate.py"

sed -i 's/any(1)/any(axis=1)/g' "$PATCH_FILE_ANNOTATE"
sed -i 's/any( 1)/any(axis=1)/g' "$PATCH_FILE_ANNOTATE"  # safe fallback

########################################
# PATCH 2 — replace get_clean_csv_df to avoid KeyError: ['fragment']
########################################

PATCH_FILE_RES="${PATCH_DIR}/plannotate/resources.py"

python - << 'EOF'
import os, re, textwrap

# Get path to patched resources.py from PATCH_DIR
patch_dir = os.environ.get("PATCH_DIR")
if not patch_dir:
    raise SystemExit("PATCH_DIR not set in environment")

path = os.path.join(patch_dir, "plannotate", "resources.py")

with open(path, "r") as f:
    txt = f.read()

# Match the whole original get_clean_csv_df function
pattern = r'def get_clean_csv_df\(recordDf\):[\\s\\S]*?^\\s*return cleaned'

new_func = """
def get_clean_csv_df(recordDf):
    \"""
    Return a cleaned DataFrame with required columns.

    Patched to be tolerant of missing columns by creating them as empty
    strings if absent. This prevents KeyError: "['fragment'] not in index".
    \"""
    required = [
        "sseqid", "qstart", "qend", "sframe", "pident", "slen", "length",
        "abs percmatch", "fragment", "db", "Feature", "Type", "Description", "qseq"
    ]
    for col in required:
        if col not in recordDf.columns:
            recordDf[col] = ""
    return recordDf[required]
"""

txt2, n = re.subn(
    pattern,
    textwrap.dedent(new_func),
    txt,
    flags=re.MULTILINE
)

if n == 0:
    raise SystemExit("Failed to patch get_clean_csv_df — pattern not found in resources.py")

with open(path, "w") as f:
    f.write(txt2)
EOF

########################################
# Activate patched version
########################################

export PYTHONPATH="${PATCH_DIR}:${PYTHONPATH}"

echo "pLannotate patched and PYTHONPATH updated:"
echo "  Using $(python - << 'EOF'
import plannotate, os
print(os.path.dirname(plannotate.__file__))
EOF
)"
echo


usage () {
cat << EOF

########## HELP ########
plasmid_mapper_PE.sh ${version}

OPTIONS:
 -h  Show help
 -1  Read1 fastq
 -2  Read2 fastq
 -r  Reference FASTA (must end in .fa or .fasta)
 -s  Skip SPAdes
 -q  Minimum Q for trimming (default 20)
 -c  Keep full SPAdes (cleanSPADES=0)
 -m  Minimum contig length
 -u  Unique ID prefix
 -k  Skip consensus calling
 -z  Assemble unmapped reads
 -y  Enable Unicycler assembly

EOF
exit
}

# parse options
while getopts "1:2:r:q:m:u:syckz" o; do
    case "${o}" in
        1) fastQ_f=${OPTARG} ;;
        2) fastQ_r=${OPTARG} ;;
        r) fastaRef=${OPTARG} ;;
        s) skipSpadesDenovoAssembly=1 ;;
        c) cleanSPADES=0 ;;
        q) qval=${OPTARG} ;;
        m) minLengthSPADES=${OPTARG} ;;
        u) uniqueID=${OPTARG} ;;
        y) runUnicycler=1 ;;
        k) skipping="TRUE" ;;
        z) assembleUnmappedReads="TRUE" ;;
        *) usage ;;
    esac
done
shift $((OPTIND-1))

echo "Running plasmid-seq with:"
echo "fastQ_f: $fastQ_f"
echo "fastQ_r: $fastQ_r"
echo "fastaRef: $fastaRef"
echo "skipSpadesDenovoAssembly: $skipSpadesDenovoAssembly"
echo "cleanSPADES: $cleanSPADES"
echo "qval: $qval"
echo "minLengthSPADES: $minLengthSPADES"
echo "uniqueID: $uniqueID"
echo "runUnicycler: $runUnicycler"
echo "assembleUnmappedReads: $assembleUnmappedReads"

### validate required inputs
if [ -z "${fastQ_f}" ] || [ -z "${fastQ_r}" ]; then usage; fi

if [ -z "${fastaRef}" ]; then
    skipping="TRUE"
elif [ "${fastaRef}" == "na" ]; then
    skipping="TRUE"
else
    if [[ "${fastaRef##*.}" != "fa" && "${fastaRef##*.}" != "fasta" ]]; then
        echo "Reference must end in .fa or .fasta"; exit 1
    fi
    skipping="FALSE"
fi

fastaRefFile=${fastaRef}
bowtieIndex=${fastaRefFile%.fa*}

########################################
# 1. Trim reads with fastp
########################################

fastp -w 4 -q ${qval} \
  -i ${fastQ_f} \
  -I ${fastQ_r} \
  -o "${fastQ_f%_S[0-9]*}_R1_trimmed.fastq" \
  -O "${fastQ_r%_S[0-9]*}_R2_trimmed.fastq" \
  -h "${fastQ_f%_S[0-9]*}_fastp_report.html" \
  -j "${fastQ_f%_S[0-9]*}_fastp_report.json" \
  -R "${fastQ_f%_S[0-9]*}_fastp_report" \
  --adapter_sequence CTGTCTCTTATACAC \
  --adapter_sequence_r2 TGTATAAGAGACAG \
  -x -y

########################################
# 2–4. Mapping to reference (optional)
########################################
if [ "${skipping}" == "FALSE" ]; then

    bowtie2-build -f ${fastaRef} ${bowtieIndex} &> "${fastQ_f%_S[0-9]*}_bt2_index.log"

    bowtie2-align-s -x ${bowtieIndex} \
       --very-sensitive \
       -1 "${fastQ_f%_S[0-9]*}_R1_trimmed.fastq" \
       -2 "${fastQ_r%_S[0-9]*}_R2_trimmed.fastq" \
       -S "${fastQ_f%_S[0-9]*}.sam" &> "${fastQ_f%_S[0-9]*}_bowtie2.log"

    sambamba view -q -F "not(secondary_alignment)" -S -f bam "${fastQ_f%_S[0-9]*}.sam" \
      | sambamba sort -q /dev/stdin -o "${fastQ_f%_S[0-9]*}.bam"
fi
########################################
# 5. SPAdes de novo assembly
########################################

if [ ${skipSpadesDenovoAssembly} == 0 ]; then

    spades_dir="${fastQ_f%_S[0-9]*}_SPADES_assembly"
    mkdir ${spades_dir}

    plasmidspades.py \
      -1 "${fastQ_f%_S[0-9]*}_R1_trimmed.fastq" \
      -2 "${fastQ_r%_S[0-9]*}_R2_trimmed.fastq" \
      -o ${spades_dir} &> "${spades_dir}/${spades_dir}.log"

    rm -r "${spades_dir}/tmp" 2>/dev/null || true

    if [ ${cleanSPADES} == 0 ]; then

        mkdir ${spades_dir/_assembly/_essential}

        find . -name "simplified_contigs.fasta" -size +1 | while read line; do
            n=${line/\/simplified/_simplified}
            n=${n##*/}
            mv ${line} ${line%/*}/${n}
        done

        find ${spades_dir} \( -name "*simplified_contigs.fasta" -o -name "*.log" \) \
          -exec cp {} ${spades_dir/_assembly/_essential} \;

        rm -r ${spades_dir}

        if [ ${minLengthSPADES} -gt 0 ]; then
            for f in ${spades_dir/_assembly/_essential}/*simplified_contigs.fasta; do

                outName="${f%_simplified_contigs.*}_SPADES_${minLengthSPADES}bp.fasta"

                seqkit seq -m ${minLengthSPADES} ${f} > ${outName}

                if [[ ! -z "$uniqueID" ]]; then
                    outNameUnique=${outName/\/K/\/${uniqueID}_K}
                    mv ${outName} ${outNameUnique}
                    export uniqueID
                    perl -i -pe 's/\>/\>$ENV{uniqueID}_/g' ${outNameUnique}
                fi

            done
        fi
    fi
fi


########################################
# 6. Unicycler de novo assembly (patched)
########################################

if [ "${runUnicycler}" == 1 ]; then
    unicycler_dir="${fastQ_f%_S[0-9]*}_unicycler_assembly"
    mkdir -p ${unicycler_dir}

    echo "Running Unicycler ..."

    unicycler \
      -1 "${fastQ_f%_S[0-9]*}_R1_trimmed.fastq" \
      -2 "${fastQ_r%_S[0-9]*}_R2_trimmed.fastq" \
      -t 4 \
      --mode bold \
      --keep 0 \
      --verbosity 0 \
      -o ${unicycler_dir} &> "${unicycler_dir}/${unicycler_dir}.log" || \
      echo "⚠️ WARNING: Unicycler exited with an error for ${uniqueID}. Continuing."

    # === PATCH — prevent crashes if Unicycler output is missing ===
    if [ -f "${unicycler_dir}/assembly.fasta" ]; then
        mv "${unicycler_dir}/assembly.fasta" \
           "${unicycler_dir}/${fastQ_f%_S[0-9]*}_unicycler_assembly.fasta"

        [ -f "${unicycler_dir}/assembly.gfa" ] && \
        mv "${unicycler_dir}/assembly.gfa" \
           "${unicycler_dir}/${fastQ_f%_S[0-9]*}_unicycler_assembly.gfa"
    else
        echo "⚠️ No Unicycler contigs produced for ${uniqueID}. Skipping renaming."
    fi
fi


########################################
# 7. Consensus FASTA (optional)
########################################

if [ "${skipConsenus}" == "FALSE" ]; then

    bcftools mpileup -Ou -f ${fastaRef} "${fastQ_f%_S[0-9]*}.bam" | \
    bcftools call -mv -Oz -o "${fastaRef/.fa/.vcf.gz}"

    tabix "${fastaRef/.fa/.vcf.gz}"

    bcftools consensus "${fastaRef/.fa/.vcf.gz}" < ${fastaRef} \
      > "${fastaRef/.fa/_consensus.fa}"

    if [[ ! -z "$uniqueID" ]]; then
        mv "${fastaRef/.fa/_consensus.fa}" \
           "${uniqueID}_${fastaRef/.fa/_consensus.fa}"
    fi
fi


########################################
# 8. Assemble unmapped reads
########################################

if [ "${assembleUnmappedReads}" == "TRUE" ]; then

    unmapped_reads_dir="${fastQ_f%_S[0-9]*}_unmapped_reads"
    mkdir -p ${unmapped_reads_dir}

    samtools view -@ 2 -f4 -bh "${fastQ_f%_S[0-9]*}.bam" > \
      "${fastQ_f%_S[0-9]*}.unmapped.bam"

    samtools fastq -@ 2 "${fastQ_f%_S[0-9]*}.unmapped.bam" \
        -1 "${fastQ_f%_S[0-9]*}.unmapped.R1.fastq.gz" \
        -2 "${fastQ_f%_S[0-9]*}.unmapped.R2.fastq.gz" \
        -0 /dev/null -s /dev/null -n

    unicycler \
      -1 "${fastQ_f%_S[0-9]*}.unmapped.R1.fastq.gz" \
      -2 "${fastQ_f%_S[0-9]*}.unmapped.R2.fastq.gz" \
      -t 4 \
      --mode normal \
      --keep 0 \
      --verbosity 0 \
      -o ${unmapped_reads_dir} &> "${unmapped_reads_dir}/${unmapped_reads_dir}.log"

    if [ -f "${unmapped_reads_dir}/assembly.fasta" ]; then
        mv "${unmapped_reads_dir}/assembly.fasta" \
           "${unmapped_reads_dir}/${fastQ_f%_S[0-9]*}_unicycler_assembly.fasta"
    fi

    mv "${fastQ_f%_S[0-9]*}.unmapped.R1.fastq.gz" \
       "${fastQ_f%_S[0-9]*}.unmapped.R2.fastq.gz" \
       "${fastQ_f%_S[0-9]*}.unmapped.bam" \
       ${unmapped_reads_dir}
fi


########################################
# 10. Run PLAnnotate (safe mode)
########################################

if [ "${runPLAnnotate}" -eq 1 ]; then
    echo "=== Running pLannotate on assemblies ==="

    plannotate_dir="${fastQ_f%_S[0-9]*}_plannotate"
    mkdir -p "${plannotate_dir}"

    assemblies=()

    # consensus
    if [ "${skipConsenus}" == "FALSE" ]; then
        cfa="${fastaRef/.fa/_consensus.fa}"
        [ -n "$uniqueID" ] && cfa="${uniqueID}_${cfa}"
        [ -f "$cfa" ] && assemblies+=("$cfa")
    fi

    # SPAdes
    spades_ess_dir="${fastQ_f%_S[0-9]*}_SPADES_essential"
    if [ -d "$spades_ess_dir" ]; then
        for f in "$spades_ess_dir"/*.fasta; do assemblies+=("$f"); done
    fi

    # Unicycler
    udir="${fastQ_f%_S[0-9]*}_unicycler_assembly"
    if [ -d "$udir" ]; then
        for f in "$udir"/*.fasta; do assemblies+=("$f"); done
    fi

    # unmapped assemblies
    uudir="${fastQ_f%_S[0-9]*}_unmapped_reads"
    if [ -d "$uudir" ]; then
        for f in "$uudir"/*.fasta; do assemblies+=("$f"); done
    fi

    if [ "${#assemblies[@]}" -eq 0 ]; then
        echo "No assemblies found for pLannotate."
    else
        for fa in "${assemblies[@]}"; do
            base=$(basename "${fa}")
            sample="${base%.*}"

            sample_dir="${plannotate_dir}/${sample}"
            mkdir -p "${sample_dir}"

            one_contig_fa="${sample_dir}/${sample}_largest_contig.fasta"

            # extract largest contig safely
            awk '
                /^>/ {if(seq){print header "\n" seq | "cat >&2"} header=$0; seq=""}
                /^[^>]/ {seq=seq $0}
                END {print header "\n" seq}
            ' "$fa" | \
            awk '
                /^>/ {h=$0; next}
                {len=length($0); print len "\t"h"\t"$0}
            ' | sort -nrk1 | head -n1 | \
            awk -F"\t" '{print $2"\n"$3}' > "${one_contig_fa}"

            if [ ! -s "${one_contig_fa}" ]; then
                echo "⚠️ pLannotate skipped: no contigs in ${fa}"
                continue
            fi

            plannotate batch \
              -i "${one_contig_fa}" \
              -o "${sample_dir}" \
              -f "${sample}" \
              -c -x -l || \
              echo "⚠️ pLannotate failed for ${sample}"
        done
    fi
fi


########################################
# 11. Convert pLannotate → MultiQC
########################################

if [ "${runPLAnnotate}" -eq 1 ]; then
    mqc_data_dir="${fastQ_f%_S[0-9]*}_multiqc_data"
    mkdir -p "${mqc_data_dir}"

    python "${plannotateParser}" "${plannotate_dir}" "${mqc_data_dir}"

    echo "pLannotate MultiQC data → ${mqc_data_dir}"
fi


########################################
# 12. MultiQC
########################################
if [ "${runMultiQC}" -eq 1 ]; then
    multiqc . -o "${fastQ_f%_S[0-9]*}_multiqc_report"
fi


########################################
# 13. Cleanup
########################################

pigz -p 12 "${fastQ_f%_S[0-9]*}_R1_trimmed.fastq" \
     "${fastQ_r%_S[0-9]*}_R2_trimmed.fastq"

find . -name "*.sam" -delete
find . -name "*.bt2" -delete
find . -name "0" -delete

echo "Pipeline complete."

