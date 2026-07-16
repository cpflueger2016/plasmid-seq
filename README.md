# plasmid-seq v1.1

Illumina plasmid QC + reference agreement pipeline for Slurm.

## Local execution (no Slurm)

The repository also provides a local runner for a workstation or Mac. It runs the
same prepare, mapping, and reporting stages synchronously, with bounded local
parallelism and no scheduler dependency.

Create the main environment (Apple Silicon is supported through Conda packages):

```bash
mamba env create -f environment.yml
conda activate plasmidseq
```

Then run:

```bash
scripts/plasmidseq_run_local.sh \
  -d /path/to/fastqs \
  -t references/PL_to_fasta.tsv \
  -f /path/to/fasta_references \
  -o /path/to/local-runs \
  -j 1 -J 4
```

`-j` is the number of samples processed concurrently and `-J` is the number of
threads for each sample. Start with `-j 1 -J 4` on a laptop because assembly is
memory-intensive. Results, staged inputs, and logs are retained below the
selected output directory. Use `-w` to generate the optional plate-map summary.

The optional annotation paths are deliberately separate:

```bash
mamba env create -f environment-plannotate.yml
mamba env create -f environment-snpeff.yml
```

They are not required for standard local QC. Enable variants with `-V`; to include
pLannotate on assembled samples, pass the pLannotate environment prefix with `-P`.
The main local runner currently leaves snpEff annotation disabled because it needs a
run-specific database assembled from the references; its environment is supplied for
manual or future extension use.

## What it does

1. Stage FASTQs + references to scratch (`prep` job)
2. Run per-sample map/assembly/QC in Slurm array (`map` jobs)
3. Gather results, build per-sample and run-level summaries (`gather` job)

## Quick start

Basic run:

```bash
scripts/plasmidseq_submit.sh \
  -d /path/to/fastqs \
  -t /group/llshared/PlasmidSeq/PL_to_fasta.tsv \
  -w /group/llshared/PlasmidSeq/PL_to_plate_position.csv \
  -p 100
```

Enable SNP/INDEL + snpEff:

```bash
scripts/plasmidseq_submit.sh \
  -d /path/to/fastqs \
  -t /group/llshared/PlasmidSeq/PL_to_fasta.tsv \
  -w /group/llshared/PlasmidSeq/PL_to_plate_position.csv \
  -p 100 \
  -V -E
```

Runtime note:
- `-V -E` typically adds about **~1 hour** to a run (depends on cluster load and sample count/read depth).

## Required inputs

- FASTQ directory (nested project folders supported)
  - expects paired files containing `_R1_` and `_R2_` and ending in `.fastq.gz`
- `PL_to_fasta.tsv`
  - column 1: `PL####`
  - column 2: reference fasta filename
- Reference FASTA directory

Optional:
- Plate map CSV (`PLid,plate,position`) for run summary plate view

## Core submit options

```text
-d <dir>   input fastq directory (required)
-t <file>  PL_to_fasta.tsv
-f <dir>   fasta reference root
-w <file>  plate map CSV (enables run summary)
-p <int>   max concurrent array tasks
-J <int>   CPUs per mapper task
-c <file>  config file override
-l <file>  submit log path
-v         print version
-V         enable VarScan
-N         disable VarScan
-E         enable snpEff (implies variants)
-S         disable snpEff
```

## Configuration

Defaults live in:
- `scripts/plasmidseq.config`

Recommended local override:
- `scripts/plasmidseq.local.config`

New notification option:

```bash
: "${SLURM_NOTIFY_EMAIL:=your.name@institute.org}"
```

When set:
- prep job mails on `BEGIN,FAIL`
- gather job mails on `END,FAIL`

## Output layout

Run root:

```text
/group/llshared/PlasmidSeq/Results/plasmidSeq_YYYY-MM-DD/<prep_jobid>/
```

Per-sample folders are organized into:
- `FASTQ/`
- `BAM/` (including compressed mpileup)
- `SNP_INDEL/`
- `Summary/` (fastp report, sample report, snpEff summary)
- `Logs/`

Notes:
- trimmed FASTQs are removed during gather cleanup
- unicycler assembly directory is kept in place

## Reports

### Per-sample report
Generated automatically in gather:
- `<sample>/Summary/<sample>_sample_report.html`
- `<sample>/Summary/<sample>_sample_report.json`

### Run summary
Generated automatically when `-w` is provided:
- `run_summary.csv`
- `run_summary.html`

Current run summary features:
- issue status aligned to per-sample QC traffic light
- plate map cells hyperlink to per-sample QC HTML
- table hyperlinks to per-sample QC HTML
- table filters: search, issue, sample QC, plate, min reads, max reads

Manual regeneration:

```bash
python3 scripts/plasmidseq_run_summary.py \
  -r /group/llshared/PlasmidSeq/Results/plasmidSeq_YYYY-MM-DD/<prep_jobid> \
  -m /group/llshared/PlasmidSeq/PL_to_plate_position.csv
```

## Monitoring

```bash
squeue -u $USER
sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed
```

## Common failure mode

If log says a PL has no destination directory (for example `PL19096`), usually:
- PL exists in `PL_to_fasta.tsv`
- but no matching FASTQ pair exists in the input dataset

Fix by syncing TSV and FASTQ set for that run.

## Version

```bash
scripts/plasmidseq_submit.sh -v
# plasmid-seq version 1.1
```
