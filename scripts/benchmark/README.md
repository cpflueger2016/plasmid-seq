# Benchmark Toolkit

This folder contains scripts to benchmark plasmid-seq runtime and memory across mapper core counts.

## Files

- `benchmark_submit.sh`: submits benchmark runs and writes a manifest.
- `collect_benchmark_metrics.py`: collects SLURM timing/memory and read-depth metrics.
- `plot_benchmark.py`: generates publication-ready plots and a summary table.

## 1) Submit benchmark runs

Example matrix:
- cores: `1,2,3,4,6,8`
- modes: `baseline,variants` (`variants` means `-V -E`)
- replicates: `3`

```bash
scripts/benchmark/benchmark_submit.sh \
  -d /group/llshared/sequencing_data/MiSeqi100/260130_SH00808_SC2160790_004/fastqs \
  -t /group/llshared/PlasmidSeq/PL_to_fasta.tsv \
  -w /group/llshared/PlasmidSeq/PL_to_plate_position.csv \
  -p 100 \
  -C 1,2,3,4,6,8 \
  -r 3
```

Outputs:
- manifest TSV in `scripts/benchmark/runs/<timestamp>/benchmark_manifest.tsv`
- submit logs in `scripts/benchmark/runs/<timestamp>/logs/`

Use `-n` for dry-run command generation.

## 2) Collect metrics

```bash
python scripts/benchmark/collect_benchmark_metrics.py \
  -m scripts/benchmark/runs/<timestamp>/benchmark_manifest.tsv \
  -o scripts/benchmark/runs/<timestamp>/benchmark_metrics.csv
```

Collected metrics include:
- prep/map/gather states and timings
- mapper task memory (MaxRSS)
- run wall-clock estimate (prep start -> gather end)
- total reads and samples with `fastp` JSON
- normalized runtime (`sec_per_million_reads`)

## 3) Generate plots

```bash
python scripts/benchmark/plot_benchmark.py \
  -i scripts/benchmark/runs/<timestamp>/benchmark_metrics.csv \
  -o scripts/benchmark/runs/<timestamp>/figures
```

Figures produced:
- `runtime_vs_cores.(png|pdf)`
- `sec_per_million_reads_vs_cores.(png|pdf)`
- `peak_memory_vs_cores.(png|pdf)`
- `runtime_vs_reads.(png|pdf)`
- `benchmark_summary_table.csv`

## 4) Per-step bottleneck timing (mapper internals)

Mapper logs now include `[timing]` lines per major stage (fastp, BBMap+markdup, SPAdes, Unicycler, variants, cleanup, total sample).

Extract them into CSV:

```bash
python scripts/benchmark/collect_step_timings.py \
  -i "/group/llshared/PlasmidSeq/Results/plasmidSeq_YYYY-MM-DD/<runid>/Logs/slurm-*_*.out" \
  -o scripts/benchmark/runs/<timestamp>/step_timings.csv
```

Outputs:
- `step_timings.csv`: one row per sample/step with `elapsed_s`
- `step_timings.summary.csv`: per-step `n`, `mean`, `median`, `p95`, `max`

## Notes

- Benchmark validity depends on mapper thread scaling. This pipeline now supports variable mapper CPUs via submit option `-J`.
- Queue wait time is not included in the main compute wall-clock metric.
