#!/usr/bin/env python3
"""
Build a per-sample plasmid report (JSON + HTML) with:
- read QC and mapping metrics
- duplicate metrics
- reference read depth track (deduplicated BAM)
- unicycler assembly mapped back to reference depth track
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import statistics
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple


BBMAP_MAPPED_PCT_RE = re.compile(r"(?im)^\s*mapped:\s*.*?([0-9]+(?:\.[0-9]+)?)%")
BBMAP_UNAMB_PCT_RE = re.compile(r"(?im)^\s*unambiguous:\s*.*?([0-9]+(?:\.[0-9]+)?)%")
BBMAP_AMB_PCT_RE = re.compile(r"(?im)^\s*ambiguous:\s*.*?([0-9]+(?:\.[0-9]+)?)%")
PL_RE = re.compile(r"(PL\d{4,})")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate per-sample JSON + HTML report.")
    parser.add_argument("-s", "--sample-dir", required=True, help="Sample output folder.")
    parser.add_argument(
        "-o",
        "--output-prefix",
        default=None,
        help="Output prefix (default: <sample-dir>/<sample-name>_sample_report).",
    )
    parser.add_argument("--min-mapped-pct", type=float, default=90.0)
    parser.add_argument("--min-breadth-1x-pct", type=float, default=99.0)
    parser.add_argument("--max-duplicate-pct", type=float, default=25.0)
    parser.add_argument("--min-mean-depth", type=float, default=50.0)
    parser.add_argument("--plot-max-points", type=int, default=1200)
    return parser.parse_args()


def run_cmd(cmd: List[str], cwd: Optional[Path] = None) -> str:
    proc = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({' '.join(cmd)}): {proc.stderr.strip()}")
    return proc.stdout


def to_float_or_none(value: object) -> Optional[float]:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_fastp_json(sample_dir: Path) -> Dict[str, Optional[float]]:
    out = {
        "raw_reads_total": None,
        "reads_after_filtering": None,
        "q30_before_pct": None,
        "q30_after_pct": None,
    }
    files = sorted(sample_dir.glob("*_fastp_report.json"))
    if not files:
        return out
    try:
        data = json.loads(files[0].read_text(encoding="utf-8"))
    except Exception:
        return out
    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    raw = to_float_or_none(before.get("total_reads"))
    out["raw_reads_total"] = raw
    out["reads_after_filtering"] = to_float_or_none(after.get("total_reads"))
    q30b = to_float_or_none(before.get("q30_rate"))
    q30a = to_float_or_none(after.get("q30_rate"))
    out["q30_before_pct"] = q30b * 100 if q30b is not None else None
    out["q30_after_pct"] = q30a * 100 if q30a is not None else None
    return out


def parse_bbmap_log(sample_dir: Path) -> Dict[str, Optional[float]]:
    out = {
        "mapped_pct": None,
        "unambiguous_pct": None,
        "ambiguous_pct": None,
    }
    logs = sorted(sample_dir.glob("*_bbmap.log"))
    if not logs:
        return out
    text = logs[0].read_text(encoding="utf-8", errors="ignore")
    mapped = [float(x) for x in BBMAP_MAPPED_PCT_RE.findall(text)]
    unamb = [float(x) for x in BBMAP_UNAMB_PCT_RE.findall(text)]
    amb = [float(x) for x in BBMAP_AMB_PCT_RE.findall(text)]
    if mapped:
        out["mapped_pct"] = sum(mapped) / len(mapped)
    if unamb:
        out["unambiguous_pct"] = sum(unamb) / len(unamb)
    if amb:
        out["ambiguous_pct"] = sum(amb) / len(amb)
    return out


def parse_fasta_lengths(path: Path) -> List[Tuple[str, int]]:
    lengths: List[Tuple[str, int]] = []
    name: Optional[str] = None
    size = 0
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths.append((name, size))
                name = line[1:].split()[0]
                size = 0
            else:
                size += len(line)
    if name is not None:
        lengths.append((name, size))
    return lengths


def locate_reference_fasta(sample_dir: Path) -> Optional[Path]:
    preferred = sorted(sample_dir.glob("*_clean.fa")) + sorted(sample_dir.glob("*_clean.fasta"))
    if preferred:
        return preferred[0]
    all_refs = sorted(sample_dir.glob("*.fa")) + sorted(sample_dir.glob("*.fasta"))
    filtered = [p for p in all_refs if "unicycler" not in p.name.lower()]
    return filtered[0] if filtered else None


def locate_primary_bam(sample_dir: Path) -> Optional[Path]:
    bams = [p for p in sorted(sample_dir.glob("*.bam")) if "unicycler_vs_ref" not in p.name]
    return bams[0] if bams else None


def locate_unicycler_fasta(sample_dir: Path) -> Optional[Path]:
    files = sorted(sample_dir.glob("*_unicycler_assembly/*_unicycler_assembly.fasta"))
    return files[0] if files else None


def ensure_bam_index(bam: Path) -> None:
    bai = bam.with_suffix(bam.suffix + ".bai")
    if not bai.exists():
        run_cmd(["samtools", "index", str(bam)])


def count_reads(cmd: List[str]) -> int:
    out = run_cmd(cmd).strip()
    return int(out) if out else 0


def duplicate_metrics(bam: Path) -> Dict[str, float]:
    total_primary = count_reads(["samtools", "view", "-c", "-F", "2308", str(bam)])
    dup_primary = count_reads(["samtools", "view", "-c", "-f", "1024", "-F", "2308", str(bam)])
    dup_pct = (dup_primary * 100.0 / total_primary) if total_primary else 0.0
    return {
        "total_primary_reads": float(total_primary),
        "duplicate_reads": float(dup_primary),
        "duplicate_pct": dup_pct,
    }


def map_unicycler_to_reference(sample_dir: Path, ref_fa: Path, uni_fa: Path, sample_name: str) -> Optional[Path]:
    out_bam = sample_dir / f"{sample_name}_unicycler_vs_ref.bam"
    map_log = sample_dir / f"{sample_name}_unicycler_vs_ref.log"

    if shutil.which("minimap2") is not None:
        cmd = (
            f"minimap2 -a -x asm5 {ref_fa} {uni_fa} 2> {map_log} "
            f"| samtools view -bS - "
            f"| samtools sort -o {out_bam} -"
        )
    elif shutil.which("bbmap.sh") is not None:
        cmd = (
            f"bbmap.sh ref={ref_fa} in={uni_fa} out=stdout.sam nodisk=t overwrite=t "
            f"2> {map_log} "
            f"| samtools view -bS - "
            f"| samtools sort -o {out_bam} -"
        )
    else:
        return None

    proc = subprocess.run(cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"assembly mapping failed: {proc.stderr.strip()}")
    run_cmd(["samtools", "index", str(out_bam)])
    return out_bam


def cumulative_offsets(ref_lengths: List[Tuple[str, int]]) -> Tuple[Dict[str, int], int]:
    offsets: Dict[str, int] = {}
    cur = 0
    for name, length in ref_lengths:
        offsets[name] = cur
        cur += length
    return offsets, cur


def depth_track(bam: Path, ref_lengths: List[Tuple[str, int]]) -> List[int]:
    offsets, total = cumulative_offsets(ref_lengths)
    arr = [0] * total
    out = run_cmd(["samtools", "depth", "-aa", str(bam)])
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        ctg, pos_s, dep_s = parts[:3]
        if ctg not in offsets:
            continue
        idx = offsets[ctg] + int(pos_s) - 1
        if 0 <= idx < total:
            arr[idx] = int(dep_s)
    return arr


def coverage_stats(depths: List[int]) -> Dict[str, float]:
    if not depths:
        return {"mean_depth": 0.0, "median_depth": 0.0, "max_depth": 0.0, "at_1x_pct": 0.0, "at_10x_pct": 0.0, "at_30x_pct": 0.0}
    n = len(depths)
    return {
        "mean_depth": float(sum(depths)) / n,
        "median_depth": float(statistics.median(depths)),
        "max_depth": float(max(depths)),
        "at_1x_pct": 100.0 * sum(1 for d in depths if d >= 1) / n,
        "at_10x_pct": 100.0 * sum(1 for d in depths if d >= 10) / n,
        "at_30x_pct": 100.0 * sum(1 for d in depths if d >= 30) / n,
    }


def contiguous_regions(depths: List[int], predicate) -> List[Dict[str, int]]:
    regions: List[Dict[str, int]] = []
    start = None
    for i, d in enumerate(depths, start=1):
        if predicate(d):
            if start is None:
                start = i
        elif start is not None:
            end = i - 1
            regions.append({"start": start, "end": end, "length_bp": end - start + 1})
            start = None
    if start is not None:
        end = len(depths)
        regions.append({"start": start, "end": end, "length_bp": end - start + 1})
    return regions


def downsample(depths: List[int], max_points: int) -> List[float]:
    if len(depths) <= max_points:
        return [float(x) for x in depths]
    bin_size = len(depths) / max_points
    out: List[float] = []
    for i in range(max_points):
        s = int(i * bin_size)
        e = int((i + 1) * bin_size)
        if e <= s:
            e = s + 1
        chunk = depths[s:e]
        out.append(float(sum(chunk)) / len(chunk))
    return out


def count_vcf_records(path: Path) -> int:
    if not path.exists():
        return 0
    n = 0
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            if line.strip():
                n += 1
    return n


def parse_snpeff_impacts(path: Path) -> Dict[str, int]:
    impacts = {"HIGH": 0, "MODERATE": 0, "LOW": 0, "MODIFIER": 0}
    if not path.exists():
        return impacts
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            info = parts[7]
            ann_items = []
            for token in info.split(";"):
                if token.startswith("ANN="):
                    ann_items = token[4:].split(",")
                    break
            for ann in ann_items:
                cols = ann.split("|")
                if len(cols) > 2 and cols[2] in impacts:
                    impacts[cols[2]] += 1
    return impacts


def parse_variant_outputs(sample_dir: Path, sample_name: str) -> Dict[str, object]:
    snps_vcf = sample_dir / f"{sample_name}_varscan_snps.vcf"
    indels_vcf = sample_dir / f"{sample_name}_varscan_indels.vcf"
    merged_vcf = sample_dir / f"{sample_name}_varscan_merged.vcf"
    snpeff_vcf = sample_dir / f"{sample_name}_varscan_snpeff.vcf"

    snp_count = count_vcf_records(snps_vcf)
    indel_count = count_vcf_records(indels_vcf)
    merged_count = count_vcf_records(merged_vcf)
    snpeff_count = count_vcf_records(snpeff_vcf)
    impact_counts = parse_snpeff_impacts(snpeff_vcf)

    return {
        "enabled_outputs_present": any([snps_vcf.exists(), indels_vcf.exists(), merged_vcf.exists()]),
        "snps_vcf": str(snps_vcf.name) if snps_vcf.exists() else "",
        "indels_vcf": str(indels_vcf.name) if indels_vcf.exists() else "",
        "merged_vcf": str(merged_vcf.name) if merged_vcf.exists() else "",
        "snpeff_vcf": str(snpeff_vcf.name) if snpeff_vcf.exists() else "",
        "snpeff_status": "present" if snpeff_vcf.exists() else "missing",
        "snp_count": snp_count,
        "indel_count": indel_count,
        "total_variants": merged_count if merged_count else (snp_count + indel_count),
        "annotated_variant_count": snpeff_count,
        "snpeff_impact_counts": impact_counts,
    }


def build_status(
    mapped_pct: Optional[float],
    breadth_1x: float,
    dup_pct: float,
    mean_depth: float,
    low_or_zero_region_count: int,
    thresholds: Dict[str, float],
) -> Tuple[float, str, List[str], str]:
    score = 100.0
    flags: List[str] = []

    min_map = thresholds["min_mapped_pct"]
    min_b1x = thresholds["min_breadth_1x_pct"]
    max_dup = thresholds["max_duplicate_pct"]
    min_md = thresholds["min_mean_depth"]

    if mapped_pct is None or mapped_pct < min_map:
        score -= 30.0
        flags.append("low_mappability")
        flags.append("high_unmapped_fraction")
    if breadth_1x < min_b1x:
        score -= 30.0
        flags.append("missing_reference_segments")
        flags.append("low_breadth")
    if dup_pct > max_dup:
        score -= 15.0
        flags.append("high_duplicate_fraction")
    if mean_depth < min_md:
        score -= 20.0
        flags.append("low_depth")
    if low_or_zero_region_count > 0:
        score -= 10.0
        flags.append("coverage_gaps_present")

    score = max(0.0, min(100.0, score))
    if score >= 85:
        light = "green"
    elif score >= 65:
        light = "yellow"
    else:
        light = "red"
    # Any zero/low-coverage region should not be green; force warn at minimum.
    if low_or_zero_region_count > 0 and light == "green":
        light = "yellow"
    summary = (
        f"mapped={mapped_pct if mapped_pct is not None else 'NA'}%, "
        f"breadth>=1x={breadth_1x:.2f}%, dup={dup_pct:.2f}%, mean_depth={mean_depth:.1f}, "
        f"low_or_zero_regions={low_or_zero_region_count}"
    )
    # unique while preserving order
    dedup_flags = list(dict.fromkeys(flags))
    return score, light, dedup_flags, summary


def write_tracks_tsv(path: Path, read_depth: List[int], asm_depth: List[int]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("position\tread_depth\tassembly_depth\n")
        for i in range(len(read_depth)):
            ad = asm_depth[i] if i < len(asm_depth) else 0
            handle.write(f"{i+1}\t{read_depth[i]}\t{ad}\n")


def render_html(
    path: Path,
    report: Dict[str, object],
    read_plot: List[float],
    asm_plot: List[float],
    total_bp: int,
) -> None:
    light = report["status"]["traffic_light"]
    light_color = {"green": "#10b981", "yellow": "#f59e0b", "red": "#ef4444"}.get(light, "#9ca3af")
    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>plasmid sample report</title>
  <style>
    body {{ font-family: Helvetica, Arial, sans-serif; margin: 20px; color: #1f2937; }}
    .top {{ display: flex; gap: 14px; align-items: center; margin-bottom: 14px; }}
    .badge {{ background: {light_color}; color: white; font-weight: 700; border-radius: 999px; padding: 6px 12px; }}
    .card {{ border: 1px solid #d1d5db; border-radius: 8px; padding: 10px 12px; margin-bottom: 12px; }}
    .grid {{ display: grid; grid-template-columns: repeat(4, minmax(140px, 1fr)); gap: 8px; }}
    .k {{ font-size: 12px; color: #6b7280; }}
    .v {{ font-size: 16px; font-weight: 700; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 8px; }}
    th, td {{ border: 1px solid #d1d5db; padding: 6px; font-size: 12px; text-align: left; }}
    th {{ background: #f3f4f6; }}
    #cov {{ width: 100%; height: 300px; border: 1px solid #d1d5db; border-radius: 8px; }}
    .legend {{ font-size: 12px; margin-top: 6px; }}
  </style>
</head>
<body>
  <div class="top">
    <div class="badge">{report["status"]["traffic_light"].upper()}</div>
    <div><b>{report["sample_id"]}</b> | score={report["status"]["score"]}</div>
  </div>

  <div class="card">
    <div class="k">Summary</div>
    <div>{report["status"]["summary"]}</div>
    <div class="k">Flags: {", ".join(report["flags"]) if report["flags"] else "none"}</div>
  </div>

  <div class="card">
    <div class="grid">
      <div><div class="k">Raw Reads</div><div class="v">{report["reads"]["raw_reads_total"]}</div></div>
      <div><div class="k">Mapped %</div><div class="v">{report["mapping"]["mapped_pct"]}</div></div>
      <div><div class="k">Duplicate %</div><div class="v">{report["duplicates"]["duplicate_pct"]}</div></div>
      <div><div class="k">Breadth >=1x %</div><div class="v">{report["coverage"]["breadth"]["at_1x_pct"]}</div></div>
    </div>
  </div>

  <div class="card">
    <b>Coverage comparison on reference ({total_bp} bp)</b>
    <canvas id="cov" width="1600" height="300"></canvas>
    <div class="legend">Blue: deduplicated read depth | Orange: unicycler assembly mapped-to-reference depth</div>
  </div>

  <div class="card">
    <table>
      <tr><th>Metric</th><th>Value</th></tr>
      <tr><td>Reference Mean Depth</td><td>{report["coverage"]["reference_read_depth"]["mean_depth"]}</td></tr>
      <tr><td>Reference Median Depth</td><td>{report["coverage"]["reference_read_depth"]["median_depth"]}</td></tr>
      <tr><td>Assembly Mean Depth</td><td>{report["coverage"]["assembly_to_reference_depth"]["mean_depth"]}</td></tr>
      <tr><td>Assembly Median Depth</td><td>{report["coverage"]["assembly_to_reference_depth"]["median_depth"]}</td></tr>
      <tr><td>Zero Coverage Regions</td><td>{len(report["coverage"]["zero_coverage_regions"])}</td></tr>
      <tr><td>Low Coverage Regions</td><td>{len(report["coverage"]["low_coverage_regions"])}</td></tr>
      <tr><td>VarScan SNPs</td><td>{report["variants"]["snp_count"]}</td></tr>
      <tr><td>VarScan Indels</td><td>{report["variants"]["indel_count"]}</td></tr>
      <tr><td>Total Variants</td><td>{report["variants"]["total_variants"]}</td></tr>
      <tr><td>snpEff Status</td><td>{report["variants"]["snpeff_status"]}</td></tr>
      <tr><td>snpEff HIGH/MODERATE</td><td>{report["variants"]["snpeff_impact_counts"]["HIGH"]}/{report["variants"]["snpeff_impact_counts"]["MODERATE"]}</td></tr>
    </table>
  </div>

  <script>
    const read = {json.dumps(read_plot)};
    const asm = {json.dumps(asm_plot)};
    const c = document.getElementById("cov");
    const ctx = c.getContext("2d");
    const W = c.width, H = c.height, pad = 30;
    const maxY = Math.max(1, ...read, ...asm);
    function draw(arr, color) {{
      if (!arr.length) return;
      ctx.beginPath();
      ctx.strokeStyle = color;
      ctx.lineWidth = 2;
      for (let i = 0; i < arr.length; i++) {{
        const x = pad + (W - 2*pad) * (i / Math.max(1, arr.length - 1));
        const y = H - pad - (H - 2*pad) * (arr[i] / maxY);
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }}
      ctx.stroke();
    }}
    ctx.fillStyle = "white"; ctx.fillRect(0, 0, W, H);
    ctx.strokeStyle = "#d1d5db"; ctx.strokeRect(pad, pad, W-2*pad, H-2*pad);
    draw(read, "#2563eb");
    draw(asm, "#f97316");
  </script>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")


def main() -> int:
    args = parse_args()
    sample_dir = Path(args.sample_dir).expanduser().resolve()
    if not sample_dir.exists():
        raise SystemExit(f"sample dir not found: {sample_dir}")

    sample_name = sample_dir.name
    pl_match = PL_RE.search(sample_name)
    plasmid_id = pl_match.group(1) if pl_match else ""

    ref_fa = locate_reference_fasta(sample_dir)
    bam = locate_primary_bam(sample_dir)
    if ref_fa is None or bam is None:
        raise SystemExit(f"missing reference fasta or bam in {sample_dir}")

    if shutil.which("samtools") is None:
        raise SystemExit("samtools not found in PATH")

    out_prefix = (
        Path(args.output_prefix).expanduser().resolve()
        if args.output_prefix
        else sample_dir / f"{sample_name}_sample_report"
    )
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    ensure_bam_index(bam)

    ref_lengths = parse_fasta_lengths(ref_fa)
    if not ref_lengths:
        raise SystemExit(f"reference has no sequences: {ref_fa}")
    total_bp = sum(length for _, length in ref_lengths)

    read_depth = depth_track(bam, ref_lengths)
    read_cov = coverage_stats(read_depth)

    uni_fa = locate_unicycler_fasta(sample_dir)
    asm_bam: Optional[Path] = None
    asm_depth: List[int] = [0] * len(read_depth)
    asm_cov = {"mean_depth": 0.0, "median_depth": 0.0, "max_depth": 0.0, "at_1x_pct": 0.0, "at_10x_pct": 0.0, "at_30x_pct": 0.0}
    if uni_fa is not None and uni_fa.exists() and uni_fa.stat().st_size > 0:
        asm_bam = map_unicycler_to_reference(sample_dir, ref_fa, uni_fa, sample_name)
        if asm_bam is not None:
            asm_depth = depth_track(asm_bam, ref_lengths)
            asm_cov = coverage_stats(asm_depth)

    fastp = parse_fastp_json(sample_dir)
    bbmap = parse_bbmap_log(sample_dir)
    dup = duplicate_metrics(bam)
    variants = parse_variant_outputs(sample_dir, sample_name)

    zero_regions = contiguous_regions(read_depth, lambda d: d == 0)
    low_regions = contiguous_regions(read_depth, lambda d: d < 10)
    low_regions = [dict(r, threshold=10) for r in low_regions]

    thresholds = {
        "min_mapped_pct": args.min_mapped_pct,
        "min_breadth_1x_pct": args.min_breadth_1x_pct,
        "max_duplicate_pct": args.max_duplicate_pct,
        "min_mean_depth": args.min_mean_depth,
    }
    score, light, flags, summary = build_status(
        bbmap.get("mapped_pct"),
        read_cov["at_1x_pct"],
        dup["duplicate_pct"],
        read_cov["mean_depth"],
        len(low_regions),
        thresholds,
    )

    asm_lengths = parse_fasta_lengths(uni_fa) if uni_fa and uni_fa.exists() else []
    asm_sizes = [x[1] for x in asm_lengths]
    asm_n50 = 0
    if asm_sizes:
        sorted_sizes = sorted(asm_sizes, reverse=True)
        half = sum(sorted_sizes) / 2
        csum = 0
        for s in sorted_sizes:
            csum += s
            if csum >= half:
                asm_n50 = s
                break

    read_pairs_total = int((fastp["raw_reads_total"] or 0) / 2) if fastp["raw_reads_total"] is not None else 0
    mapped_pct = bbmap.get("mapped_pct")
    unmapped_pct = (100.0 - mapped_pct) if mapped_pct is not None else None

    cleaned_header_fasta = str(ref_fa.name) if "_clean." in ref_fa.name else ""
    source_fasta = str(ref_fa.name)
    if "_clean." in ref_fa.name:
        source_guess = sample_dir / ref_fa.name.replace("_clean", "")
        if source_guess.exists():
            source_fasta = source_guess.name

    run_id = ""
    parents = list(sample_dir.parents)
    for i, p in enumerate(parents):
        if p.name.startswith("plasmidSeq_") and i > 0:
            run_id = f"{p.name}/{parents[i-1].name}"
            break

    report: Dict[str, object] = {
        "schema_version": "1.0.0",
        "run_id": run_id,
        "sample_id": sample_name,
        "plasmid_id": plasmid_id,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "reference": {
            "name": ref_lengths[0][0],
            "length_bp": total_bp,
            "cleaned_header_fasta": cleaned_header_fasta,
            "source_fasta": source_fasta,
        },
        "status": {
            "traffic_light": light,
            "score": round(score, 2),
            "summary": summary,
        },
        "reads": {
            "raw_reads_total": int(fastp["raw_reads_total"] or 0),
            "read_pairs_total": int(read_pairs_total),
            "reads_after_filtering": int(fastp["reads_after_filtering"] or 0),
            "q30_before_pct": round(float(fastp["q30_before_pct"] or 0.0), 2),
            "q30_after_pct": round(float(fastp["q30_after_pct"] or 0.0), 2),
        },
        "mapping": {
            "tool": "bbmap",
            "mapped_pct": round(float(mapped_pct or 0.0), 2),
            "unambiguous_pct": round(float(bbmap.get("unambiguous_pct") or 0.0), 2),
            "ambiguous_pct": round(float(bbmap.get("ambiguous_pct") or 0.0), 2),
            "unmapped_pct": round(float(unmapped_pct or 0.0), 2),
            "error_rate_pct": 0.0,
        },
        "duplicates": {
            "duplicate_reads": int(dup["duplicate_reads"]),
            "duplicate_pct": round(float(dup["duplicate_pct"]), 2),
            "deduplicated_bam": str(bam.name),
            "markdup_log": f"{sample_name}_markdup.log",
        },
        "coverage": {
            "reference_read_depth": {
                "mean_depth": round(read_cov["mean_depth"], 3),
                "median_depth": round(read_cov["median_depth"], 3),
                "max_depth": round(read_cov["max_depth"], 3),
            },
            "assembly_to_reference_depth": {
                "mean_depth": round(asm_cov["mean_depth"], 3),
                "median_depth": round(asm_cov["median_depth"], 3),
                "max_depth": round(asm_cov["max_depth"], 3),
            },
            "breadth": {
                "at_1x_pct": round(read_cov["at_1x_pct"], 3),
                "at_10x_pct": round(read_cov["at_10x_pct"], 3),
                "at_30x_pct": round(read_cov["at_30x_pct"], 3),
            },
            "zero_coverage_regions": zero_regions,
            "low_coverage_regions": low_regions,
        },
        "assembly": {
            "unicycler_status": "present" if uni_fa is not None else "missing",
            "contig_count": len(asm_sizes),
            "total_length_bp": int(sum(asm_sizes)),
            "longest_contig_bp": int(max(asm_sizes)) if asm_sizes else 0,
            "n50_bp": int(asm_n50),
            "assembly_fasta": str(uni_fa.name) if uni_fa else "",
            "assembly_map_bam": str(asm_bam.name) if asm_bam else "",
        },
        "variants": variants,
        "flags": flags,
        "artifacts": {
            "html_report": str(f"{out_prefix.name}.html"),
            "json_report": str(f"{out_prefix.name}.json"),
            "coverage_tracks_tsv": str(f"{out_prefix.name}_coverage_tracks.tsv"),
        },
        "thresholds": thresholds,
    }

    tracks_path = Path(f"{out_prefix}_coverage_tracks.tsv")
    json_path = Path(f"{out_prefix}.json")
    html_path = Path(f"{out_prefix}.html")

    write_tracks_tsv(tracks_path, read_depth, asm_depth)
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    read_plot = downsample(read_depth, args.plot_max_points)
    asm_plot = downsample(asm_depth, args.plot_max_points)
    render_html(html_path, report, read_plot, asm_plot, total_bp)

    print(f"[sample-report] sample={sample_name}")
    print(f"[sample-report] json={json_path}")
    print(f"[sample-report] html={html_path}")
    print(f"[sample-report] tracks={tracks_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
