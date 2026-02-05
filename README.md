# plasmid-seq — Illumina plasmid sequencing QC + reference agreement pipeline

This pipeline processes Illumina FASTQ data generated from plasmid sequencing runs and evaluates
agreement between the sequencing-derived assembly/alignment and a provided reference FASTA.

The workflow is optimized for Slurm/HPC:
- stages inputs to scratch
- matches each plasmid (PL####) to its reference FASTA using `PL_to_fasta.tsv`
- launches mapping/assembly jobs per plasmid (Slurm array)
- gathers results into a single results folder
- cleans up scratch to keep storage tidy

---

## Contents

- [Requirements](#requirements)
- [Inputs](#inputs)
  - [FASTQ layout expectations](#fastq-layout-expectations)
  - [PL_to_fasta.tsv format](#pl_to_fastatsv-format)
  - [Reference FASTA folder](#reference-fasta-folder)
- [Quick start](#quick-start)
- [What the pipeline does](#what-the-pipeline-does)
- [Outputs](#outputs)
- [Logging](#logging)
- [Run Summary Report (CSV + HTML)](#run-summary-report-csv--html)
- [Handling missing/ambiguous references](#handling-missingambiguous-references)
- [Performance + Slurm notes](#performance--slurm-notes)
- [Troubleshooting](#troubleshooting)
- [Development / extending the pipeline](#development--extending-the-pipeline)

---

## Requirements

### Software (typical)
- Slurm (`sbatch`, `squeue`, `sacct`)
- GNU coreutils, bash
- `rsync`, `find`
- A conda env providing the mapper/assembly dependencies used by `plasmidseq_mapper_PE.sh`
  (fastp, minimap2/bwa, samtools, unicycler, etc. — depends on your mapper script)

### HPC assumptions
- A stable scratch root (often `$MYSCRATCH`) is available.
- You have permission to write to:
  - scratch staging area
  - results destination (e.g. `/group/llshared/PlasmidSeq/Results/...`)

### Configuration
See `README_config.md` for portable configuration via:
- `scripts/plasmidseq.config` (defaults)
- `scripts/plasmidseq.local.config` (site/user overrides)

---

## Inputs

### FASTQ layout expectations

The pipeline expects paired-end FASTQs named like:

```
PL18955_<something>_S194_L001_R1_001.fastq.gz
PL18955_<something>_S194_L001_R2_001.fastq.gz
```

The submit/prep staging step is designed to handle nested “project folders”.
Example:

```
/path/to/fastqs/
  ProjectA/
    PL18955_..._R1_001.fastq.gz
    PL18955_..._R2_001.fastq.gz
  ProjectB/
    PL18956_..._R1_001.fastq.gz
    PL18956_..._R2_001.fastq.gz
```

During staging, FASTQs are copied into scratch and then moved into per-plasmid folders.

> If FASTQs are not found after staging, check:
> - the input directory is correct
> - FASTQs are not nested deeper than the include rules
> - filenames end in `.fastq.gz` and include `_R1_` / `_R2_`

### PL_to_fasta.tsv format

A tab-separated file mapping each plasmid ID (PL number) to a reference FASTA filename.

Minimal format:

| Column | Meaning |
|---|---|
| 1 | PL ID (e.g. `PL18955`) |
| 2 | FASTA filename (e.g. `pCDS_B3o_IV2.fa`) |

Example:

```
PL18955  pCDS_B3o_IV2.fa
PL18956  pUBQ10-aGCN4-CENH3.fa
```

### Reference FASTA folder

A directory containing reference FASTA files (often organized in subfolders):

```
Fasta_Reference_Files/
  labA/
    pCDS_B3o_IV2.fa
  labB/
    pCDS_B3o_IV2.fa   # possible duplicate name across labs
```

The matcher searches for the filename in this folder tree.

---

## Quick start

### 1) Configure defaults (recommended)

Create a local config (not committed) on each server:

```bash
cp scripts/plasmidseq.local.config.template scripts/plasmidseq.local.config
# edit scripts/plasmidseq.local.config
```

(Or create `scripts/plasmidseq.local.config` following `README_config.md`.)

### 2) Submit a run

Basic (use defaults for TSV + refs):

```bash
scripts/plasmidseq_submit.sh -d /path/to/fastqs
```

Override mapping TSV and reference folder:

```bash
scripts/plasmidseq_submit.sh \
  -d /path/to/fastqs \
  -t /path/to/PL_to_fasta.tsv \
  -f /path/to/Fasta_Reference_Files
```

Use a specific config file and custom submit log path:

```bash
scripts/plasmidseq_submit.sh \
  -d /path/to/fastqs \
  -c scripts/plasmidseq.local.config \
  -l /path/to/plasmidseq_submit.log
```

Enable automatic run summary generation during gather:

```bash
scripts/plasmidseq_submit.sh \
  -d /path/to/fastqs \
  -w /mmfs1/data/group/llshared/PlasmidSeq/PL_to_plate_position.csv
```

Throttle concurrency (max tasks running at once):

```bash
scripts/plasmidseq_submit.sh -d /path/to/fastqs -p 50
```

### 3) Monitor

```bash
squeue -u $USER
```

If you have job IDs, track directly:

```bash
squeue -j <prep_jobid>,<array_jobid>,<gather_jobid>
```

Submit options (from `scripts/plasmidseq_submit.sh -h`):
- `-d <dir>` input FASTQ folder (required)
- `-t <file>` `PL_to_fasta.tsv` override
- `-f <dir>` reference FASTA folder override
- `-p <int>` max concurrent array tasks
- `-c <file>` config file override
- `-w <file>` plate map CSV for auto run summary
- `-l <file>` submit log file path
- `-h` help

Reference help output:

```text
Usage:
  plasmidseq_submit.sh -d <plasmidSeqData_dir> [options]

Required:
  -d <dir>   Folder containing Project folders / fastqs to stage into scratch

Optional:
  -t <file>  PL_to_fasta.tsv path (default: from config DEFAULT_TSV)
  -f <dir>   Fasta reference folder path (default: from config DEFAULT_REFS)
  -p <int>   Max concurrent array tasks (default: from config MAX_CONCURRENT_DEFAULT)
  -c <file>  Config file path (default precedence: local config, then plasmidseq.config)
  -w <file>  Plate map CSV for run summary (columns: PLid,plate,position)
  -l <file>  Submit log file (default: <plasmidSeqData>/plasmidseq_submit_<date>.log)
  -h         Show help
```

---

## What the pipeline does

At a high level:

1. **Prep job**
   - Creates run scratch folder and results folder
   - Stages `PL_to_fasta.tsv` and the reference folder into scratch
   - Copies FASTQs from the input directory into scratch
   - Moves FASTQs into per-plasmid folders
   - Runs FASTA matching (`match_plasmid_fasta_refs_v2.bash`)
   - Writes `jobs.tsv` (one line per plasmid to process)

2. **Map array job(s)**
   - One Slurm array task per plasmid
   - Runs `plasmidseq_mapper_PE.sh` in the plasmid folder using R1/R2 + reference

3. **Gather job**
   - Collects outputs into a single results directory
   - Optionally generates `run_summary.csv` + `run_summary.html` when `-w <plate_map.csv>` is provided
   - Copies key run-level files (e.g. `jobs.tsv`, `plasmid_fasta_match.log`)
   - Deletes scratch run directory (safely)

---

## Outputs

Results are organized by run date + Slurm prep job ID:

```
RESULTS_BASE/
  plasmidSeq_YYYY-MM-DD/
    <prep_jobid>/
      jobs.tsv
      plasmid_fasta_match.log
      Logs/
      <run_subfolder>/   # e.g. "3-20kb_plasmids"
        PL18955_<name>/
          ... mapper outputs ...
        PL18956_<name>/
          ... mapper outputs ...
```

You should also see per-plasmid artifacts, for example:
- trimmed FASTQs
- fastp reports
- assembly/alignment outputs
- reference FASTA (if matched/copied)
- match status files (e.g. ambiguity markers)

> Exact output files depend on `plasmidseq_mapper_PE.sh`.

---

## Logging

There are multiple logging layers:

- **submit wrapper** (`plasmidseq_submit.sh`)
  - prints a concise “submission plan” (jobdate, input dir, TSV/refs, job IDs)
  - automatically mirrors stdout/stderr to a log file (default: `<plasmidSeqData>/plasmidseq_submit_<date>.log`)
  - use `-l` to set a custom submit log path

- **prep job log** (Slurm output)
  - prints scratch and results paths, number of FASTQs staged, jobs.tsv count

- **map array logs**
  - typically configured to write array task logs under `Logs/`
  - each task prints which folder/ref it is using

- **matcher log**
  - `plasmid_fasta_match.log` in scratch (and copied into results)

---

## Run Summary Report (CSV + HTML)

Generate a per-run summary table + report page with:
- total reads per sample (from `*_fastp_report.json`)
- mapping rate (from `*_bbmap.log`, with Bowtie2 logs still supported for legacy runs)
- reference/plannotate status flags
- bar graph of reads per sample
- 96-well plate issue view (using `PL_to_plate_position.csv`)

Run:

```bash
scripts/plasmidseq_run_summary.py \
  -r /group/llshared/PlasmidSeq/Results/plasmidSeq_YYYY-MM-DD/<prep_jobid> \
  -m /mmfs1/data/group/llshared/PlasmidSeq/PL_to_plate_position.csv
```

Outputs (written into the run folder by default):
- `run_summary.csv`
- `run_summary.html`

Optional:
- `--low-map-threshold 80` sets the warning cutoff for mapping %
- `-o /path/to/output_prefix` writes `<prefix>.csv` and `<prefix>.html`

Automatic mode:
- pass `-w /path/to/PL_to_plate_position.csv` to `plasmidseq_submit.sh`
- gather will run `plasmidseq_run_summary.py` before copying to `Aligned`
- summary files are copied to `Aligned` with the rest of run results

Per-sample report mode (automatic in gather):
- gather runs `plasmidseq_sample_report.py` for each sample folder with fastp output
- outputs per sample:
  - `<sample>_sample_report.json`
  - `<sample>_sample_report.html`
  - `<sample>_sample_report_coverage_tracks.tsv`
- report compares:
  - deduplicated BBMap depth over reference
  - Unicycler assembly mapped back to reference depth
  - includes traffic-light status + QC metrics (reads, mapping, duplicates, coverage)

---

## Handling missing/ambiguous references

The FASTA matcher supports multiple states per plasmid:

- **Found reference**
  - reference FASTA copied into the plasmid folder
  - can write a provenance file (e.g. `FASTA_REF_SOURCE.txt`)

- **Missing reference**
  - marker file may be created (e.g. `na`) to indicate no reference was found
  - pipeline can still run (depending on mapper behavior), but reference alignment may be skipped

- **Ambiguous reference (duplicate FASTA filenames)**
  - if the same FASTA filename exists multiple times under the reference folder,
    the matcher records matches in:
    `FASTA_REF_AMBIGUOUS.matches.txt`

Recommended behavior:
- proceed using the **first match** and emit a **warning**
- keep the `FASTA_REF_AMBIGUOUS.matches.txt` file for auditing

---

## Performance + Slurm notes

### Concurrency (`-p`)
`-p` controls how many array tasks can run simultaneously (array throttle):
- for short plasmid jobs (~minutes), values like 50–200 are usually fine
- for heavier assemblies, lower it to avoid I/O pressure

### Array size limits
Many clusters have a maximum array size (e.g. 1001).
If you routinely run 200–300 plasmids you’re fine in a single array, but
for very large runs you may need chunking into multiple arrays.

If your site enforces `MaxArraySize=1001`, ensure the submit wrapper can split arrays
and set gather dependencies on **all** array jobs.

---

## Troubleshooting

### Prep job fails with `unbound variable`
Common when `set -u` is enabled and a variable name mismatches (e.g. `TSV` vs `tsv`).
Fix by using consistent variable names in the Slurm scripts and submit wrapper.

### `jobs.tsv exists but has 0 lines`
Usually means:
- no FASTQs were staged (bad input path or staging rules)
- FASTQs did not match expected naming pattern
- plasmid folder creation failed

Check the prep job Slurm output.

### No FASTA copied / `na` present
If the plasmid folder already contains a marker like `na` from a previous run,
re-running may “look like it skipped”.
For a clean re-run, remove the plasmid folder or re-stage into a new scratch run.

### Array tasks fail immediately
Check `Logs/%x_%A_%a.err` and `.out` for:
- missing conda env
- missing mapper script
- missing reference FASTA
- wrong working directory

---

## Development / extending the pipeline

Recommended principles:

1. Keep site-specific values in config (`scripts/plasmidseq.config` + local override)
2. Ensure each module:
   - is independently runnable
   - logs clearly (inputs, outputs, key counts)
   - fails fast on missing prerequisites
3. Keep run “glue” files at the top-level results folder:
   - `jobs.tsv`
   - `plasmid_fasta_match.log`
   - submit log (job IDs, resolved paths)
4. Prefer explicit “source tracking” files for anything auto-selected
   (e.g. `FASTA_REF_SOURCE.txt` describing which reference was used and why)

If you add new modules (e.g. coverage QC, variant calling, summary tables, MultiQC),
keep them as separate Slurm steps with well-defined inputs/outputs and a gather step.

---

## License

See `LICENSE`.
