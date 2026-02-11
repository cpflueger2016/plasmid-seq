#!/usr/bin/env python3
import argparse
import csv
import math
import os
from collections import defaultdict

import matplotlib.pyplot as plt


def to_float(x):
    try:
        if x is None or x == "":
            return None
        return float(x)
    except Exception:
        return None


def to_int(x):
    try:
        return int(x)
    except Exception:
        return None


def load_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def grouped_stats(rows, y_field):
    grouped = defaultdict(list)
    for r in rows:
        mode = r.get("mode", "")
        cores = to_int(r.get("cores", ""))
        y = to_float(r.get(y_field, ""))
        if not mode or cores is None or y is None:
            continue
        grouped[(mode, cores)].append(y)

    out = defaultdict(list)
    for (mode, cores), vals in grouped.items():
        vals = sorted(vals)
        n = len(vals)
        mean = sum(vals) / n
        if n > 1:
            var = sum((v - mean) ** 2 for v in vals) / (n - 1)
            sd = math.sqrt(var)
        else:
            sd = 0.0
        out[mode].append((cores, mean, sd, n))

    for mode in out:
        out[mode].sort(key=lambda t: t[0])
    return out


def plot_errorbar(stats, ylabel, title, out_png, out_pdf, scale=1.0):
    plt.figure(figsize=(7.5, 5.0))
    styles = {
        "baseline": {"marker": "o", "color": "#1f77b4", "label": "Baseline"},
        "variants": {"marker": "s", "color": "#d62728", "label": "Variants (-V -E)"},
    }
    for mode, series in sorted(stats.items()):
        xs = [x for x, _, _, _ in series]
        ys = [m / scale for _, m, _, _ in series]
        es = [s / scale for _, _, s, _ in series]
        sty = styles.get(mode, {"marker": "o", "color": None, "label": mode})
        plt.errorbar(xs, ys, yerr=es, capsize=4, linewidth=2, **sty)

    plt.xlabel("Mapper cores per sample")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(alpha=0.25, linestyle="--")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()


def plot_reads_scatter(rows, out_png, out_pdf):
    plt.figure(figsize=(7.5, 5.0))
    styles = {
        "baseline": {"marker": "o", "color": "#1f77b4", "label": "Baseline"},
        "variants": {"marker": "s", "color": "#d62728", "label": "Variants (-V -E)"},
    }
    for mode in sorted(set(r.get("mode", "") for r in rows)):
        xs, ys = [], []
        for r in rows:
            if r.get("mode") != mode:
                continue
            reads = to_float(r.get("reads_total", ""))
            wall = to_float(r.get("pipeline_wallclock_s", ""))
            if reads is None or wall is None:
                continue
            xs.append(reads / 1_000_000.0)
            ys.append(wall / 3600.0)
        if not xs:
            continue
        sty = styles.get(mode, {"marker": "o", "color": None, "label": mode})
        plt.scatter(xs, ys, alpha=0.7, s=45, **sty)

    plt.xlabel("Total reads per run (millions)")
    plt.ylabel("Pipeline wall clock (hours)")
    plt.title("Runtime vs input read count")
    plt.grid(alpha=0.25, linestyle="--")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()


def write_summary_table(rows, out_csv):
    fields = [
        "mode",
        "cores",
        "n_runs",
        "wallclock_mean_h",
        "wallclock_sd_h",
        "sec_per_million_reads_mean",
        "sec_per_million_reads_sd",
        "map_maxrss_max_mb_mean",
        "map_maxrss_max_mb_sd",
    ]

    grouped = defaultdict(list)
    for r in rows:
        mode = r.get("mode", "")
        cores = to_int(r.get("cores", ""))
        if not mode or cores is None:
            continue
        grouped[(mode, cores)].append(r)

    with open(out_csv, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        for (mode, cores), grp in sorted(grouped.items(), key=lambda x: (x[0][0], x[0][1])):
            w = [to_float(r.get("pipeline_wallclock_s")) for r in grp]
            w = [x for x in w if x is not None]
            spm = [to_float(r.get("sec_per_million_reads")) for r in grp]
            spm = [x for x in spm if x is not None]
            rss = [to_float(r.get("map_maxrss_max_mb")) for r in grp]
            rss = [x for x in rss if x is not None]

            def mean_sd(vals):
                if not vals:
                    return None, None
                m = sum(vals) / len(vals)
                if len(vals) == 1:
                    return m, 0.0
                var = sum((v - m) ** 2 for v in vals) / (len(vals) - 1)
                return m, math.sqrt(var)

            wm, ws = mean_sd(w)
            sm, ss = mean_sd(spm)
            rm, rs = mean_sd(rss)

            wr.writerow(
                {
                    "mode": mode,
                    "cores": cores,
                    "n_runs": len(grp),
                    "wallclock_mean_h": None if wm is None else wm / 3600.0,
                    "wallclock_sd_h": None if ws is None else ws / 3600.0,
                    "sec_per_million_reads_mean": sm,
                    "sec_per_million_reads_sd": ss,
                    "map_maxrss_max_mb_mean": rm,
                    "map_maxrss_max_mb_sd": rs,
                }
            )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="metrics CSV from collect_benchmark_metrics.py")
    ap.add_argument("-o", "--outdir", required=True, help="output directory for figures/tables")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rows = load_rows(args.input)

    runtime_stats = grouped_stats(rows, "pipeline_wallclock_s")
    if runtime_stats:
        plot_errorbar(
            runtime_stats,
            ylabel="Wall clock runtime (hours)",
            title="Plasmid-seq runtime vs mapper cores",
            out_png=os.path.join(args.outdir, "runtime_vs_cores.png"),
            out_pdf=os.path.join(args.outdir, "runtime_vs_cores.pdf"),
            scale=3600.0,
        )

    spm_stats = grouped_stats(rows, "sec_per_million_reads")
    if spm_stats:
        plot_errorbar(
            spm_stats,
            ylabel="Seconds per million reads",
            title="Normalized runtime vs mapper cores",
            out_png=os.path.join(args.outdir, "sec_per_million_reads_vs_cores.png"),
            out_pdf=os.path.join(args.outdir, "sec_per_million_reads_vs_cores.pdf"),
            scale=1.0,
        )

    mem_stats = grouped_stats(rows, "map_maxrss_max_mb")
    if mem_stats:
        plot_errorbar(
            mem_stats,
            ylabel="Map task peak RSS (GB)",
            title="Peak mapper memory vs cores",
            out_png=os.path.join(args.outdir, "peak_memory_vs_cores.png"),
            out_pdf=os.path.join(args.outdir, "peak_memory_vs_cores.pdf"),
            scale=1024.0,
        )

    plot_reads_scatter(
        rows,
        out_png=os.path.join(args.outdir, "runtime_vs_reads.png"),
        out_pdf=os.path.join(args.outdir, "runtime_vs_reads.pdf"),
    )

    write_summary_table(rows, os.path.join(args.outdir, "benchmark_summary_table.csv"))
    print(f"[bench] wrote figures/table to: {args.outdir}")


if __name__ == "__main__":
    main()
