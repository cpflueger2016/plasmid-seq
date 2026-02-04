#!/usr/bin/env python3
"""
Create per-run CSV + HTML summary for plasmid-seq outputs.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


PL_RE = re.compile(r"(PL\d{4,})")
OVERALL_ALIGN_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)% overall alignment rate")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate plasmid-seq run summary CSV + HTML."
    )
    parser.add_argument(
        "-r",
        "--run-dir",
        required=True,
        help="Run directory (e.g. .../Results/plasmidSeq_YYYY-MM-DD/<prep_jobid>).",
    )
    parser.add_argument(
        "-m",
        "--plate-map",
        default="/mmfs1/data/group/llshared/PlasmidSeq/PL_to_plate_position.csv",
        help="Plate map CSV with columns: PLid,plate,position (default: shared file).",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        default=None,
        help="Output path prefix (default: <run-dir>/run_summary).",
    )
    parser.add_argument(
        "--low-map-threshold",
        type=float,
        default=80.0,
        help="Warn when overall mapping percent is below this value (default: 80).",
    )
    return parser.parse_args()


def to_float_or_none(value: object) -> Optional[float]:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def load_plate_map(path: Path) -> Dict[str, Tuple[str, str]]:
    mapping: Dict[str, Tuple[str, str]] = {}
    if not path.exists():
        return mapping

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return mapping
        key_map = {name.strip().lower(): name for name in reader.fieldnames}
        pl_key = key_map.get("plid") or key_map.get("pl_id") or key_map.get("pl")
        plate_key = key_map.get("plate")
        pos_key = key_map.get("position") or key_map.get("well")
        if not (pl_key and plate_key and pos_key):
            return mapping

        for row in reader:
            plid = (row.get(pl_key) or "").strip()
            plate = (row.get(plate_key) or "").strip()
            position = (row.get(pos_key) or "").strip().upper()
            if not plid:
                continue
            if plate.lower() == "na" or position.lower() == "na" or not position:
                continue
            mapping[plid] = (plate, position)
    return mapping


def parse_jobs_file(run_dir: Path) -> List[Dict[str, str]]:
    jobs_file = run_dir / "jobs.tsv"
    rows: List[Dict[str, str]] = []
    if not jobs_file.exists():
        return rows

    with jobs_file.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            while len(parts) < 5:
                parts.append("")
            folder, ref, r1, r2, plid = parts[:5]
            rows.append(
                {
                    "folder": folder,
                    "ref": ref,
                    "r1": r1,
                    "r2": r2,
                    "plid": plid,
                }
            )
    return rows


def parse_fastp_json(sample_dir: Path) -> Dict[str, Optional[float]]:
    result = {
        "raw_reads_total": None,
        "reads_after_filtering": None,
        "q30_before": None,
        "q30_after": None,
    }
    files = sorted(sample_dir.glob("*_fastp_report.json"))
    if not files:
        return result

    try:
        with files[0].open(encoding="utf-8") as handle:
            data = json.load(handle)
    except Exception:
        return result

    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    result["raw_reads_total"] = to_float_or_none(before.get("total_reads"))
    result["reads_after_filtering"] = to_float_or_none(after.get("total_reads"))
    result["q30_before"] = to_float_or_none(before.get("q30_rate"))
    result["q30_after"] = to_float_or_none(after.get("q30_rate"))
    return result


def parse_bowtie_log(sample_dir: Path) -> Dict[str, Optional[float]]:
    result = {
        "bowtie2_overall_alignment_pct": None,
    }
    files = sorted(sample_dir.glob("*_bowtie2.log"))
    if not files:
        return result

    text = files[0].read_text(encoding="utf-8", errors="ignore")
    match = OVERALL_ALIGN_RE.search(text)
    if match:
        result["bowtie2_overall_alignment_pct"] = float(match.group(1))
    return result


def detect_reference_status(sample_dir: Path, ref_name_from_jobs: str) -> str:
    if (sample_dir / "FASTA_REF_FOUND").exists():
        if (sample_dir / "FASTA_REF_AMBIGUOUS").exists():
            return "ambiguous"
        return "found"
    if (sample_dir / "na").exists() or not ref_name_from_jobs.strip():
        return "missing"
    return "unknown"


def detect_unicycler_status(sample_dir: Path) -> str:
    if list(sample_dir.glob("*_unicycler_assembly/*_unicycler_assembly.fasta")):
        return "present"
    return "missing"


def detect_plannotate_status(sample_dir: Path) -> str:
    if list(sample_dir.glob("*_plannotate")):
        return "present"
    return "missing"


def discover_extra_sample_dirs(run_dir: Path, known_folders: set[str]) -> List[Dict[str, str]]:
    extras: List[Dict[str, str]] = []
    for fastp_json in run_dir.glob("**/*_fastp_report.json"):
        sample_dir = fastp_json.parent
        rel = str(sample_dir.relative_to(run_dir))
        if rel in known_folders:
            continue
        match = PL_RE.search(sample_dir.name)
        plid = match.group(1) if match else ""
        extras.append({"folder": rel, "ref": "", "r1": "", "r2": "", "plid": plid})
    return extras


def build_rows(
    run_dir: Path,
    job_rows: List[Dict[str, str]],
    plate_map: Dict[str, Tuple[str, str]],
    low_map_threshold: float,
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    known = {row["folder"] for row in job_rows}
    all_rows = list(job_rows) + discover_extra_sample_dirs(run_dir, known)

    for job in all_rows:
        folder = job["folder"]
        sample_dir = run_dir / folder
        sample_name = Path(folder).name
        pl_match = PL_RE.search(sample_name)
        plid = job["plid"].strip() or (pl_match.group(1) if pl_match else "")

        fastp = parse_fastp_json(sample_dir) if sample_dir.exists() else {}
        bowtie = parse_bowtie_log(sample_dir) if sample_dir.exists() else {}
        ref_status = detect_reference_status(sample_dir, job["ref"]) if sample_dir.exists() else "missing"
        unicycler_status = detect_unicycler_status(sample_dir) if sample_dir.exists() else "missing"
        plannotate_status = detect_plannotate_status(sample_dir) if sample_dir.exists() else "missing"

        plate = ""
        well = ""
        if plid and plid in plate_map:
            plate, well = plate_map[plid]

        issues: List[str] = []
        severity = "ok"

        if not sample_dir.exists():
            severity = "fail"
            issues.append("sample_folder_missing")

        if fastp.get("raw_reads_total") is None:
            severity = "fail"
            issues.append("fastp_report_missing")

        mapped = bowtie.get("bowtie2_overall_alignment_pct")
        if ref_status in {"found", "ambiguous"} and mapped is None:
            severity = "fail"
            issues.append("bowtie2_log_missing")

        if ref_status == "missing":
            if severity != "fail":
                severity = "warn"
            issues.append("reference_missing")
        elif ref_status == "ambiguous":
            if severity != "fail":
                severity = "warn"
            issues.append("reference_ambiguous")

        if mapped is not None and mapped < low_map_threshold:
            if severity == "ok":
                severity = "warn"
            issues.append(f"low_mapping<{low_map_threshold:g}")

        if plannotate_status != "present":
            if severity == "ok":
                severity = "warn"
            issues.append("plannotate_missing")

        raw_reads = fastp.get("raw_reads_total")
        reads_after = fastp.get("reads_after_filtering")
        q30_before = fastp.get("q30_before")
        q30_after = fastp.get("q30_after")

        row: Dict[str, object] = {
            "sample_folder": folder,
            "sample_name": sample_name,
            "pl_id": plid,
            "plate": plate,
            "well": well,
            "reference_from_jobs": job["ref"],
            "reference_status": ref_status,
            "raw_reads_total": int(raw_reads) if raw_reads is not None else "",
            "read_pairs_total": int(raw_reads / 2) if raw_reads is not None else "",
            "reads_after_filtering": int(reads_after) if reads_after is not None else "",
            "q30_before_pct": round((q30_before or 0) * 100, 2) if q30_before is not None else "",
            "q30_after_pct": round((q30_after or 0) * 100, 2) if q30_after is not None else "",
            "bowtie2_overall_alignment_pct": round(mapped, 2) if mapped is not None else "",
            "unicycler_status": unicycler_status,
            "plannotate_status": plannotate_status,
            "issue_flag": severity,
            "issue_reason": ";".join(issues),
        }
        rows.append(row)

    rows.sort(key=lambda r: (str(r["pl_id"]), str(r["sample_name"])))
    return rows


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    fieldnames = [
        "sample_folder",
        "sample_name",
        "pl_id",
        "plate",
        "well",
        "reference_from_jobs",
        "reference_status",
        "raw_reads_total",
        "read_pairs_total",
        "reads_after_filtering",
        "q30_before_pct",
        "q30_after_pct",
        "bowtie2_overall_alignment_pct",
        "unicycler_status",
        "plannotate_status",
        "issue_flag",
        "issue_reason",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def fmt_number(value: object) -> str:
    if value in ("", None):
        return ""
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return str(value)


def make_bar_chart(rows: List[Dict[str, object]]) -> str:
    with_reads = [r for r in rows if r["raw_reads_total"] != ""]
    if not with_reads:
        return "<p>No read counts found.</p>"

    with_reads.sort(key=lambda r: int(r["raw_reads_total"]), reverse=True)
    max_reads = max(int(r["raw_reads_total"]) for r in with_reads)

    parts = ['<div class="bars">']
    for row in with_reads:
        reads = int(row["raw_reads_total"])
        width = (reads / max_reads) * 100 if max_reads else 0
        label = html.escape(str(row["pl_id"] or row["sample_name"]))
        parts.append(
            (
                '<div class="bar-row">'
                f'<div class="bar-label">{label}</div>'
                f'<div class="bar-track"><div class="bar-fill" style="width:{width:.2f}%"></div></div>'
                f'<div class="bar-value">{reads:,}</div>'
                "</div>"
            )
        )
    parts.append("</div>")
    return "\n".join(parts)


def make_plate_grids(rows: List[Dict[str, object]], plate_map: Dict[str, Tuple[str, str]]) -> str:
    rows_by_pl = {str(r["pl_id"]): r for r in rows if r["pl_id"]}
    plate_wells: Dict[str, Dict[str, str]] = defaultdict(dict)
    for plid, (plate, well) in plate_map.items():
        plate_wells[plate][well] = plid

    if not plate_wells:
        return "<p>Plate map not available.</p>"

    out: List[str] = []
    row_labels = list("ABCDEFGH")
    col_labels = [str(i) for i in range(1, 13)]

    for plate in sorted(plate_wells):
        out.append(f"<h3>{html.escape(plate)}</h3>")
        out.append('<table class="plate">')
        out.append("<tr><th></th>" + "".join(f"<th>{c}</th>" for c in col_labels) + "</tr>")

        for r in row_labels:
            cells = [f"<th>{r}</th>"]
            for c in col_labels:
                well = f"{r}{c}"
                plid = plate_wells[plate].get(well, "")
                if not plid:
                    cells.append('<td class="cell-empty">-</td>')
                    continue
                sample = rows_by_pl.get(plid)
                if sample is None:
                    cls = "cell-na"
                    title = f"{plid}: not in run output"
                    text = plid
                else:
                    state = str(sample["issue_flag"])
                    cls = {"ok": "cell-ok", "warn": "cell-warn", "fail": "cell-fail"}.get(state, "cell-na")
                    reason = str(sample["issue_reason"])
                    mapping = sample["bowtie2_overall_alignment_pct"]
                    reads = sample["raw_reads_total"]
                    title = f"{plid} | reads={reads} | map%={mapping} | {reason}"
                    text = plid
                cells.append(f'<td class="{cls}" title="{html.escape(title)}">{html.escape(text)}</td>')
            out.append("<tr>" + "".join(cells) + "</tr>")
        out.append("</table>")

    return "\n".join(out)


def render_html(
    path: Path,
    run_dir: Path,
    rows: List[Dict[str, object]],
    plate_map: Dict[str, Tuple[str, str]],
    low_map_threshold: float,
) -> None:
    total = len(rows)
    n_ok = sum(1 for r in rows if r["issue_flag"] == "ok")
    n_warn = sum(1 for r in rows if r["issue_flag"] == "warn")
    n_fail = sum(1 for r in rows if r["issue_flag"] == "fail")

    bar_chart = make_bar_chart(rows)
    plate_html = make_plate_grids(rows, plate_map)

    table_rows: List[str] = []
    for r in rows:
        table_rows.append(
            "<tr>"
            f"<td>{html.escape(str(r['pl_id']))}</td>"
            f"<td>{html.escape(str(r['sample_folder']))}</td>"
            f"<td>{html.escape(str(r['plate']))}</td>"
            f"<td>{html.escape(str(r['well']))}</td>"
            f"<td>{fmt_number(r['raw_reads_total'])}</td>"
            f"<td>{fmt_number(r['reads_after_filtering'])}</td>"
            f"<td>{html.escape(str(r['bowtie2_overall_alignment_pct']))}</td>"
            f"<td>{html.escape(str(r['reference_status']))}</td>"
            f"<td>{html.escape(str(r['plannotate_status']))}</td>"
            f"<td class='flag-{html.escape(str(r['issue_flag']))}'>{html.escape(str(r['issue_flag']))}</td>"
            f"<td>{html.escape(str(r['issue_reason']))}</td>"
            "</tr>"
        )

    html_text = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>plasmid-seq run summary</title>
  <style>
    body {{ font-family: Helvetica, Arial, sans-serif; margin: 20px; color: #1f2937; }}
    h1, h2, h3 {{ margin-bottom: 8px; }}
    .meta {{ margin-bottom: 18px; color: #4b5563; }}
    .cards {{ display: flex; gap: 12px; flex-wrap: wrap; margin: 16px 0; }}
    .card {{ border: 1px solid #d1d5db; border-radius: 8px; padding: 10px 12px; min-width: 120px; background: #f9fafb; }}
    .bars {{ width: 100%; }}
    .bar-row {{ display: grid; grid-template-columns: 250px 1fr 110px; gap: 8px; align-items: center; margin: 3px 0; }}
    .bar-label {{ font-size: 12px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }}
    .bar-track {{ height: 12px; border-radius: 6px; background: #e5e7eb; }}
    .bar-fill {{ height: 100%; border-radius: 6px; background: #2563eb; }}
    .bar-value {{ font-family: monospace; font-size: 12px; text-align: right; }}
    .legend {{ display: flex; gap: 12px; margin: 10px 0; font-size: 12px; }}
    .pill {{ padding: 3px 8px; border-radius: 999px; border: 1px solid #d1d5db; }}
    .ok {{ background: #dcfce7; }}
    .warn {{ background: #fef3c7; }}
    .fail {{ background: #fee2e2; }}
    .na {{ background: #e5e7eb; }}
    table {{ border-collapse: collapse; width: 100%; margin: 12px 0; }}
    th, td {{ border: 1px solid #d1d5db; padding: 5px 6px; font-size: 12px; }}
    th {{ background: #f3f4f6; position: sticky; top: 0; }}
    .table-wrap {{ max-height: 520px; overflow: auto; border: 1px solid #d1d5db; }}
    .plate {{ width: auto; margin-bottom: 20px; }}
    .plate th, .plate td {{ text-align: center; min-width: 52px; }}
    .cell-ok {{ background: #86efac; }}
    .cell-warn {{ background: #fde68a; }}
    .cell-fail {{ background: #fca5a5; }}
    .cell-na {{ background: #d1d5db; }}
    .cell-empty {{ background: #f3f4f6; color: #9ca3af; }}
    .flag-ok {{ color: #166534; font-weight: 700; }}
    .flag-warn {{ color: #92400e; font-weight: 700; }}
    .flag-fail {{ color: #991b1b; font-weight: 700; }}
  </style>
</head>
<body>
  <h1>plasmid-seq run summary</h1>
  <div class="meta">
    run_dir: {html.escape(str(run_dir))}<br>
    low mapping warning threshold: {low_map_threshold:g}%
  </div>

  <div class="cards">
    <div class="card"><b>Samples</b><br>{total}</div>
    <div class="card"><b>OK</b><br>{n_ok}</div>
    <div class="card"><b>Warn</b><br>{n_warn}</div>
    <div class="card"><b>Fail</b><br>{n_fail}</div>
  </div>

  <h2>Reads Per Sample (bar graph)</h2>
  {bar_chart}

  <h2>96-Well Plate Issue View</h2>
  <div class="legend">
    <span class="pill ok">OK</span>
    <span class="pill warn">Warn</span>
    <span class="pill fail">Fail</span>
    <span class="pill na">No run data</span>
  </div>
  {plate_html}

  <h2>Per-Sample Metrics</h2>
  <div class="table-wrap">
    <table>
      <thead>
        <tr>
          <th>PL</th><th>Sample Folder</th><th>Plate</th><th>Well</th>
          <th>Raw Reads</th><th>Reads After Filter</th><th>Map %</th>
          <th>Ref Status</th><th>pLannotate</th><th>Issue</th><th>Reason</th>
        </tr>
      </thead>
      <tbody>
        {"".join(table_rows)}
      </tbody>
    </table>
  </div>
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).expanduser().resolve()
    if not run_dir.exists():
        raise SystemExit(f"Run directory not found: {run_dir}")

    out_prefix = (
        Path(args.output_prefix).expanduser().resolve()
        if args.output_prefix
        else run_dir / "run_summary"
    )
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    plate_map = load_plate_map(Path(args.plate_map).expanduser())
    job_rows = parse_jobs_file(run_dir)
    rows = build_rows(run_dir, job_rows, plate_map, args.low_map_threshold)

    csv_path = Path(f"{out_prefix}.csv")
    html_path = Path(f"{out_prefix}.html")
    write_csv(csv_path, rows)
    render_html(html_path, run_dir, rows, plate_map, args.low_map_threshold)

    print(f"[summary] samples={len(rows)}")
    print(f"[summary] csv={csv_path}")
    print(f"[summary] html={html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
