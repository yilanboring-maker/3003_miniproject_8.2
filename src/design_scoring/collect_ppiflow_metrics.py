#!/usr/bin/env python3
"""Collect PPIFlow batch metrics into reusable ranking tables.

This is a first-pass structural/design-quality screen. It does not replace
PRODIGY binding-affinity scoring; PRODIGY can be joined later by sample name.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from collections import defaultdict
from pathlib import Path


def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def parse_float(value: str, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def rank_score(row: dict[str, str]) -> float:
    """Higher is better for this first-pass filter."""
    hotspot = parse_float(row.get("hotspot_coverage", "0"), 0.0)
    cdr_ratio = parse_float(row.get("cdr_interface_ratio", "0"), 0.0)
    dsasa = parse_float(row.get("dsasa", "0"), 0.0)
    rmsd = parse_float(row.get("rmsd_framework", "10"), 10.0)
    has_clash = parse_bool(row.get("has_clash", "False"))

    score = 100.0 * hotspot + 50.0 * cdr_ratio + 0.02 * dsasa - 10.0 * rmsd
    if has_clash:
        score -= 100.0
    return score


def collect_metrics(out_root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    metric_files = sorted(out_root.glob("*/*/sample_metrics.csv"))
    total = len(metric_files)
    start = time.time()

    for index, metric_file in enumerate(metric_files, start=1):
        modality = metric_file.parents[1].name
        batch_id = metric_file.parent.name
        elapsed = time.time() - start
        percent = 100.0 * index / max(total, 1)
        print(
            f"[collect] {index}/{total} {percent:.1f}% elapsed={elapsed:.1f}s "
            f"file={metric_file}",
            flush=True,
        )

        with metric_file.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                sample = row.get("sample", "")
                record: dict[str, str] = {
                    "design_id": sample.removesuffix(".pdb"),
                    "sample": sample,
                    "modality": modality,
                    "batch_id": batch_id,
                    "pdb_path": row.get("pdb_path", ""),
                    "hotspot_coverage": row.get("hotspot_coverage", ""),
                    "cdr_interface_ratio": row.get("cdr_interface_ratio", ""),
                    "dsasa": row.get("dsasa", ""),
                    "rmsd_framework": row.get("rmsd_framework", ""),
                    "has_clash": row.get("has_clash", ""),
                    "coverage_hotspot_list": row.get("coverage_hotspot_list", ""),
                    "interface_residues": row.get("interface_residues", ""),
                    "cdr_interface_residues": row.get("cdr_interface_residues", ""),
                    "first_pass_score": f"{rank_score(row):.4f}",
                    "prodigy_delta_g": "",
                    "prodigy_kd": "",
                    "status": "first_pass",
                    "notes": "PPIFlow structural metrics only; PRODIGY not yet joined",
                }
                rows.append(record)

    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "rank",
        "design_id",
        "sample",
        "modality",
        "batch_id",
        "pdb_path",
        "hotspot_coverage",
        "cdr_interface_ratio",
        "dsasa",
        "rmsd_framework",
        "has_clash",
        "first_pass_score",
        "prodigy_delta_g",
        "prodigy_kd",
        "coverage_hotspot_list",
        "interface_residues",
        "cdr_interface_residues",
        "status",
        "notes",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-root", required=True, help="PPIFlow generation output root")
    parser.add_argument("--results-dir", required=True, help="Directory for collected outputs")
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    out_root = Path(args.out_root)
    results_dir = Path(args.results_dir)
    checkpoint_path = results_dir / "06_collect_checkpoint.json"
    all_results_path = results_dir / "all_results.csv"
    partial_top_path = results_dir / "partial_top.csv"
    summary_path = results_dir / "summary.json"

    print(f"[start] out_root={out_root}", flush=True)
    print(f"[start] results_dir={results_dir}", flush=True)
    rows = collect_metrics(out_root)

    by_modality: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_modality[row["modality"]].append(row)

    ranked_all: list[dict[str, str]] = []
    top_rows: list[dict[str, str]] = []
    summary: dict[str, object] = {
        "out_root": str(out_root),
        "results_dir": str(results_dir),
        "total_designs": len(rows),
        "top_n_per_modality": args.top_n,
        "modalities": {},
        "updated_at_epoch": int(time.time()),
    }

    for modality, modality_rows in sorted(by_modality.items()):
        modality_rows.sort(key=lambda item: parse_float(item["first_pass_score"], -9999), reverse=True)
        clash_free = sum(1 for row in modality_rows if not parse_bool(row.get("has_clash", "False")))
        full_hotspot = sum(1 for row in modality_rows if parse_float(row.get("hotspot_coverage", "0"), 0.0) >= 0.999)
        summary["modalities"][modality] = {
            "count": len(modality_rows),
            "clash_free": clash_free,
            "full_hotspot_coverage": full_hotspot,
        }
        for rank, row in enumerate(modality_rows, start=1):
            row = dict(row)
            row["rank"] = str(rank)
            ranked_all.append(row)
            if rank <= args.top_n:
                top_rows.append(row)

    write_csv(all_results_path, ranked_all)
    write_csv(partial_top_path, top_rows)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    checkpoint_path.write_text(
        json.dumps(
            {
                "status": "completed",
                "all_results": str(all_results_path),
                "partial_top": str(partial_top_path),
                "summary": str(summary_path),
                "total_designs": len(rows),
                "updated_at_epoch": int(time.time()),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    print(f"[done] total_designs={len(rows)}", flush=True)
    print(f"[done] all_results={all_results_path}", flush=True)
    print(f"[done] partial_top={partial_top_path}", flush=True)
    print(f"[done] summary={summary_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
