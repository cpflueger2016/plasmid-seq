#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re
import statistics
from typing import Dict, List


STEP_END_RE = re.compile(
    r"\[timing\]\s+step=(?P<step>\S+)\s+status=end\s+elapsed_s=(?P<elapsed>\d+)"
)
SAMPLE_START_RE = re.compile(r"\[timing\]\s+sample=(?P<sample>\S+)\s+status=start")
SAMPLE_END_RE = re.compile(
    r"\[timing\]\s+sample=(?P<sample>\S+)\s+status=end\s+elapsed_s=(?P<elapsed>\d+)"
)


def p95(vals: List[float]) -> float:
    if len(vals) == 1:
        return vals[0]
    vals_sorted = sorted(vals)
    idx = int(round(0.95 * (len(vals_sorted) - 1)))
    return vals_sorted[idx]


def collect_rows(log_files: List[str]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for log_file in log_files:
        sample = ""
        try:
            with open(log_file, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    m_sample_start = SAMPLE_START_RE.search(line)
                    if m_sample_start:
                        sample = m_sample_start.group("sample")

                    m_step = STEP_END_RE.search(line)
                    if m_step:
                        rows.append(
                            {
                                "log_file": log_file,
                                "sample": sample,
                                "step": m_step.group("step"),
                                "elapsed_s": m_step.group("elapsed"),
                            }
                        )

                    m_sample_end = SAMPLE_END_RE.search(line)
                    if m_sample_end:
                        rows.append(
                            {
                                "log_file": log_file,
                                "sample": m_sample_end.group("sample"),
                                "step": "sample_total",
                                "elapsed_s": m_sample_end.group("elapsed"),
                            }
                        )
        except OSError:
            continue
    return rows


def write_rows(rows: List[Dict[str, str]], out_csv: str) -> None:
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=["log_file", "sample", "step", "elapsed_s"])
        wr.writeheader()
        wr.writerows(rows)


def write_summary(rows: List[Dict[str, str]], out_csv: str) -> None:
    by_step: Dict[str, List[float]] = {}
    for r in rows:
        try:
            x = float(r["elapsed_s"])
        except (KeyError, ValueError):
            continue
        by_step.setdefault(r["step"], []).append(x)

    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(
            fh,
            fieldnames=["step", "n", "mean_s", "median_s", "p95_s", "max_s"],
        )
        wr.writeheader()
        for step in sorted(by_step):
            vals = by_step[step]
            wr.writerow(
                {
                    "step": step,
                    "n": len(vals),
                    "mean_s": f"{statistics.mean(vals):.3f}",
                    "median_s": f"{statistics.median(vals):.3f}",
                    "p95_s": f"{p95(vals):.3f}",
                    "max_s": f"{max(vals):.3f}",
                }
            )


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Collect per-step timing from mapper [timing] log lines."
    )
    ap.add_argument(
        "-i",
        "--input-glob",
        required=True,
        help="Glob for SLURM log files (e.g., '/path/Logs/slurm-*_*.out')",
    )
    ap.add_argument(
        "-o",
        "--output-csv",
        required=True,
        help="Output CSV path for raw step timings.",
    )
    args = ap.parse_args()

    files = sorted(glob.glob(args.input_glob))
    if not files:
        raise SystemExit(f"No files matched glob: {args.input_glob}")

    rows = collect_rows(files)
    write_rows(rows, args.output_csv)

    root, ext = os.path.splitext(args.output_csv)
    summary_csv = f"{root}.summary{ext or '.csv'}"
    write_summary(rows, summary_csv)
    print(f"Wrote {len(rows)} timing rows: {args.output_csv}")
    print(f"Wrote summary: {summary_csv}")


if __name__ == "__main__":
    main()
