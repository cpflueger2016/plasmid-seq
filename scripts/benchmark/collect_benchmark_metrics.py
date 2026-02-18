#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
import statistics
import subprocess
from datetime import datetime
from typing import Dict, List, Optional, Tuple


def run_cmd(cmd: List[str]) -> str:
    p = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if p.returncode != 0:
        return ""
    return p.stdout


def parse_manifest(path: str) -> List[Dict[str, str]]:
    rows = []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            rows.append(r)
    return rows


def parse_elapsed_raw(v: str) -> Optional[float]:
    v = (v or "").strip()
    if not v:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def parse_maxrss_mb(v: str) -> Optional[float]:
    v = (v or "").strip()
    if not v:
        return None
    m = re.match(r"^([0-9.]+)([KMGTP]?)$", v)
    if not m:
        return None
    num = float(m.group(1))
    unit = m.group(2)
    scale = {
        "": 1.0 / (1024.0 * 1024.0),
        "K": 1.0 / 1024.0,
        "M": 1.0,
        "G": 1024.0,
        "T": 1024.0 * 1024.0,
        "P": 1024.0 * 1024.0 * 1024.0,
    }[unit]
    return num * scale


def parse_iso(ts: str) -> Optional[datetime]:
    ts = (ts or "").strip()
    if not ts:
        return None
    try:
        return datetime.fromisoformat(ts)
    except Exception:
        return None


def sacct_rows(job_id: str) -> List[Dict[str, str]]:
    out = run_cmd([
        "sacct",
        "-P",
        "-n",
        "-j",
        job_id,
        "--format=JobIDRaw,JobName,State,ExitCode,ElapsedRaw,MaxRSS,Start,End",
    ])
    rows = []
    for line in out.splitlines():
        parts = line.split("|")
        if len(parts) != 8:
            continue
        rows.append(
            {
                "JobIDRaw": parts[0].strip(),
                "JobName": parts[1].strip(),
                "State": parts[2].strip(),
                "ExitCode": parts[3].strip(),
                "ElapsedRaw": parts[4].strip(),
                "MaxRSS": parts[5].strip(),
                "Start": parts[6].strip(),
                "End": parts[7].strip(),
            }
        )
    return rows


def aggregate_map_jobs(map_jobids: str) -> Dict[str, Optional[float]]:
    if not map_jobids:
        return {
            "map_tasks_total": 0,
            "map_tasks_failed": 0,
            "map_elapsed_mean_s": None,
            "map_elapsed_p95_s": None,
            "map_elapsed_max_s": None,
            "map_maxrss_mean_mb": None,
            "map_maxrss_p95_mb": None,
            "map_maxrss_max_mb": None,
            "map_last_end": None,
        }

    task_elapsed = []
    task_rss = []
    task_states = []
    task_ends = []

    for jobid in map_jobids.split(":"):
        jobid = jobid.strip()
        if not jobid:
            continue
        for row in sacct_rows(jobid):
            jid = row["JobIDRaw"]
            if not re.match(r"^[0-9]+_[0-9]+$", jid):
                continue
            task_states.append(row["State"])
            e = parse_elapsed_raw(row["ElapsedRaw"])
            if e is not None:
                task_elapsed.append(e)
            r = parse_maxrss_mb(row["MaxRSS"])
            if r is not None:
                task_rss.append(r)
            end = parse_iso(row["End"])
            if end is not None:
                task_ends.append(end)

    def p95(vals: List[float]) -> Optional[float]:
        if not vals:
            return None
        if len(vals) == 1:
            return vals[0]
        vals = sorted(vals)
        idx = int(round(0.95 * (len(vals) - 1)))
        return vals[idx]

    failed = sum(1 for s in task_states if s and not s.startswith("COMPLETED"))
    last_end = max(task_ends).isoformat() if task_ends else None

    return {
        "map_tasks_total": len(task_states),
        "map_tasks_failed": failed,
        "map_elapsed_mean_s": statistics.mean(task_elapsed) if task_elapsed else None,
        "map_elapsed_p95_s": p95(task_elapsed),
        "map_elapsed_max_s": max(task_elapsed) if task_elapsed else None,
        "map_maxrss_mean_mb": statistics.mean(task_rss) if task_rss else None,
        "map_maxrss_p95_mb": p95(task_rss),
        "map_maxrss_max_mb": max(task_rss) if task_rss else None,
        "map_last_end": last_end,
    }


def prep_or_gather_metrics(jobid: str) -> Dict[str, Optional[str]]:
    if not jobid:
        return {
            "state": "",
            "elapsed_s": None,
            "maxrss_mb": None,
            "start": None,
            "end": None,
        }
    rows = sacct_rows(jobid)
    root = None
    for r in rows:
        if r["JobIDRaw"] == jobid:
            root = r
            break
    if root is None and rows:
        root = rows[0]
    if root is None:
        return {
            "state": "",
            "elapsed_s": None,
            "maxrss_mb": None,
            "start": None,
            "end": None,
        }
    return {
        "state": root["State"],
        "elapsed_s": parse_elapsed_raw(root["ElapsedRaw"]),
        "maxrss_mb": parse_maxrss_mb(root["MaxRSS"]),
        "start": root["Start"] or None,
        "end": root["End"] or None,
    }


def read_fastp_totals(base_dir: str) -> Tuple[int, int]:
    total_reads = 0
    n_files = 0
    for root, _, files in os.walk(base_dir):
        for fn in files:
            if not fn.endswith("_fastp_report.json"):
                continue
            path = os.path.join(root, fn)
            try:
                with open(path) as fh:
                    data = json.load(fh)
                before = data.get("summary", {}).get("before_filtering", {})
                reads = int(before.get("total_reads", 0))
                total_reads += reads
                n_files += 1
            except Exception:
                continue
    return total_reads, n_files


def compute_wallclock_s(prep_start: Optional[str], gather_end: Optional[str], map_last_end: Optional[str]) -> Optional[float]:
    st = parse_iso(prep_start or "")
    en = parse_iso(gather_end or "")
    if st and en:
        return (en - st).total_seconds()
    if st and map_last_end:
        me = parse_iso(map_last_end)
        if me:
            return (me - st).total_seconds()
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("-m", "--manifest", required=True, help="benchmark_manifest.tsv")
    ap.add_argument("-o", "--output", required=True, help="output metrics CSV")
    args = ap.parse_args()

    rows = parse_manifest(args.manifest)
    out_rows = []

    for r in rows:
        prep = prep_or_gather_metrics(r.get("prep_jobid", ""))
        gather = prep_or_gather_metrics(r.get("gather_jobid", ""))
        mapm = aggregate_map_jobs(r.get("map_jobids", ""))

        results = (r.get("results") or "").strip()
        scratch = (r.get("scratch") or "").strip()

        reads_total = 0
        sample_count = 0
        if results and os.path.isdir(results):
            reads_total, sample_count = read_fastp_totals(results)
        if reads_total == 0 and scratch and os.path.isdir(scratch):
            reads_total, sample_count = read_fastp_totals(scratch)

        wallclock = compute_wallclock_s(prep.get("start"), gather.get("end"), mapm.get("map_last_end"))
        sec_per_million = None
        if wallclock and reads_total > 0:
            sec_per_million = wallclock / (reads_total / 1_000_000.0)

        out_rows.append(
            {
                "run_label": r.get("run_label", ""),
                "mode": r.get("mode", ""),
                "cores": r.get("cores", ""),
                "replicate": r.get("replicate", ""),
                "submit_exit": r.get("submit_exit", ""),
                "prep_jobid": r.get("prep_jobid", ""),
                "map_jobids": r.get("map_jobids", ""),
                "gather_jobid": r.get("gather_jobid", ""),
                "prep_state": prep.get("state", ""),
                "prep_elapsed_s": prep.get("elapsed_s"),
                "prep_maxrss_mb": prep.get("maxrss_mb"),
                "map_tasks_total": mapm.get("map_tasks_total"),
                "map_tasks_failed": mapm.get("map_tasks_failed"),
                "map_elapsed_mean_s": mapm.get("map_elapsed_mean_s"),
                "map_elapsed_p95_s": mapm.get("map_elapsed_p95_s"),
                "map_elapsed_max_s": mapm.get("map_elapsed_max_s"),
                "map_maxrss_mean_mb": mapm.get("map_maxrss_mean_mb"),
                "map_maxrss_p95_mb": mapm.get("map_maxrss_p95_mb"),
                "map_maxrss_max_mb": mapm.get("map_maxrss_max_mb"),
                "gather_state": gather.get("state", ""),
                "gather_elapsed_s": gather.get("elapsed_s"),
                "gather_maxrss_mb": gather.get("maxrss_mb"),
                "pipeline_wallclock_s": wallclock,
                "reads_total": reads_total,
                "samples_with_fastp_json": sample_count,
                "sec_per_million_reads": sec_per_million,
                "results_dir": results,
                "scratch_dir": scratch,
            }
        )

    fieldnames = list(out_rows[0].keys()) if out_rows else []
    with open(args.output, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=fieldnames)
        wr.writeheader()
        for row in out_rows:
            wr.writerow(row)

    print(f"[bench] wrote metrics: {args.output}")


if __name__ == "__main__":
    main()
