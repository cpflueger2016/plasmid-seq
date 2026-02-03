# Plasmid-seq pipeline configuration

This pipeline is designed to run on HPC via Slurm and to be portable across servers.
To avoid hard-coding site-specific paths (shared storage, conda env locations, partitions, etc.),
the pipeline supports a simple shell config file.

## Config files

### Tracked default config (safe to commit)
- `scripts/plasmidseq.config`

This should contain **sane defaults** for your “main” environment (e.g., Kaya / llshared).

### Local override config (do NOT commit)
- `scripts/plasmidseq.local.config`

This is for site- or user-specific overrides:
- conda env location
- module load commands
- partitions / account / qos
- alternate `RESULTS_BASE` / `DEFAULT_TSV` / `DEFAULT_REFS`
- alternate script paths if the pipeline is installed elsewhere

Add to `.gitignore`:
```bash
scripts/plasmidseq.local.config
```

## Config precedence

From highest priority to lowest:

1. Command line: `-c /path/to/config` (if supported by your submit wrapper)
2. Environment variable: `PLASMIDSEQ_CONFIG=/path/to/config`
3. `scripts/plasmidseq.local.config` (if it exists)
4. `scripts/plasmidseq.config` (default)

> Note: If your current `plasmidseq_submit.sh` does not yet support `-c`, you can still use
> `PLASMIDSEQ_CONFIG` for the Slurm scripts that source config, and/or edit the defaults in `plasmidseq_submit.sh`.
> The recommended path is to add `-c` support to `plasmidseq_submit.sh` so everything is consistent.

## Variables

### Core paths
- `RESULTS_BASE`  
  Where final run outputs are stored, typically:
  `/group/llshared/PlasmidSeq/Results`

- `DEFAULT_TSV`  
  Default `PL_to_fasta.tsv` mapping file path.

- `DEFAULT_REFS`  
  Default reference FASTA folder containing the reference plasmid sequences.

- `SCRATCH_BASE`  
  Scratch root. Many sites expose a stable `$MYSCRATCH`. You can set:
  `SCRATCH_BASE="${MYSCRATCH:-/scratch/${USER}}"`

### Pipeline scripts / tool paths
These are useful when the scripts are installed outside the repo or in a shared location.

- `PIPELINE_SCRIPTS_DIR`  
  Folder containing the pipeline scripts.

- `MATCHER_PATH`  
  Path to `match_plasmid_fasta_refs_v2.bash` (optional override).

- `MAPPER_PATH`  
  Path to `plasmidseq_mapper_PE.sh` (optional override).

### Conda + modules (site portability)
- `MODULES`  
  Space-separated modules to load inside Slurm jobs (e.g. `MODULES="gcc samtools"`)

- `CONDA_INIT`  
  Path to conda init script for non-interactive shells, e.g.:
  - `$HOME/miniconda3/etc/profile.d/conda.sh`
  - `/path/to/miniconda/etc/profile.d/conda.sh`

- `CONDA_ENV`  
  Conda env name **or** absolute path to env:
  - `CONDA_ENV="plasmidseq"`
  - `CONDA_ENV="/group/llshared/shared_conda_envs/plasmidseq"`

> Recommended pattern in Slurm scripts:
> - load modules (if needed)
> - `source "$CONDA_INIT"` (if set)
> - `conda activate "$CONDA_ENV"` (if set)

### Slurm defaults
- `PARTITION_PREP`, `PARTITION_MAP`, `PARTITION_GATHER`
- `SBATCH_ACCOUNT`, `SBATCH_QOS` (optional)

### Performance knobs
- `MAX_CONCURRENT_DEFAULT`  
  Default array throttle (e.g. 50 means `--array=...%50`)

- `MAX_ARRAY_SIZE`  
  Cluster maximum Slurm array size (Kaya typically `1001`)

- `MAP_CPUS_PER_TASK`, `MAP_MEM_PER_CPU`  
  Mapping job resources. If you submit with `sbatch --cpus-per-task ...`,
  CLI overrides `#SBATCH` headers.

## Example: local config template

Create `scripts/plasmidseq.local.config`:

```bash
RESULTS_BASE="/group/llshared/PlasmidSeq/Results"
DEFAULT_TSV="/group/llshared/PlasmidSeq/PL_to_fasta.tsv"
DEFAULT_REFS="/group/llshared/PlasmidSeq/Fasta_Reference_Files"

SCRATCH_BASE="${MYSCRATCH:-/scratch/${USER}}"

MODULES="gcc"
CONDA_INIT=""
CONDA_ENV="/group/llshared/shared_conda_envs/plasmidseq"

PARTITION_PREP="ll"
PARTITION_MAP="work"
PARTITION_GATHER="ll"

MAX_CONCURRENT_DEFAULT=50
MAX_ARRAY_SIZE=1001
MAP_CPUS_PER_TASK=4
MAP_MEM_PER_CPU="6G"
```

## Quick verification

Inside an interactive shell on the target server:

```bash
# Check defaults resolve
source scripts/plasmidseq.config
echo "$DEFAULT_TSV"
echo "$DEFAULT_REFS"
echo "$RESULTS_BASE"
```

And for Slurm scripts that source config, you can test:

```bash
PLASMIDSEQ_CONFIG=$(pwd)/scripts/plasmidseq.config \
sbatch --test-only scripts/plasmidseq_prepare_SLURM.sh -h
```
