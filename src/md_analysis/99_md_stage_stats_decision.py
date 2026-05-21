#!/usr/bin/env python3
"""Summarize staged GROMACS MD analysis and assign interval-level decisions."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


METRIC_FILES = {
    "rmsd_complex": "rmsd_complex.xvg",
    "rmsd_antibody": "rmsd_antibody.xvg",
    "rmsd_antigen": "rmsd_antigen.xvg",
    "hbonds": "hbonds_antibody_antigen.xvg",
    "interface_contacts": "interface_contacts.xvg",
    "interface_mindist": "interface_mindist.xvg",
    "buried_sasa": "buried_sasa.csv",
    "rg_complex": "rg_complex.xvg",
    "rg_antibody": "rg_antibody.xvg",
    "rg_antigen": "rg_antigen.xvg",
    "sasa_complex": "sasa_complex.xvg",
    "sasa_antibody": "sasa_antibody.xvg",
    "sasa_antigen": "sasa_antigen.xvg",
    "nis_proxy": "nis_proxy.csv",
}


def read_xy(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    if not path.exists():
        return rows
    for line in path.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    return rows


def percentile(values: list[float], q: float) -> float | None:
    vals = sorted(v for v in values if math.isfinite(v))
    if not vals:
        return None
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return vals[lo]
    return vals[lo] + (vals[hi] - vals[lo]) * (pos - lo)


def slope(rows: list[tuple[float, float]]) -> float | None:
    pairs = [(x, y) for x, y in rows if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return None
    xs = [x for x, _ in pairs]
    ys = [y for _, y in pairs]
    xbar = sum(xs) / len(xs)
    ybar = sum(ys) / len(ys)
    den = sum((x - xbar) ** 2 for x in xs)
    if den == 0:
        return None
    return sum((x - xbar) * (y - ybar) for x, y in pairs) / den


def stats(rows: list[tuple[float, float]], last_window_ns: float = 2.0) -> dict[str, float | None]:
    vals = [value for _, value in rows if math.isfinite(value)]
    if not vals:
        return {
            "mean": None,
            "median": None,
            "min": None,
            "max": None,
            "p95": None,
            "last": None,
            "last2ns_mean": None,
            "trend_slope": None,
            "n": 0,
        }
    max_time = max(t for t, _ in rows)
    last_window = [v for t, v in rows if t >= max_time - last_window_ns and math.isfinite(v)]
    return {
        "mean": sum(vals) / len(vals),
        "median": percentile(vals, 0.50),
        "min": min(vals),
        "max": max(vals),
        "p95": percentile(vals, 0.95),
        "last": vals[-1],
        "last2ns_mean": (sum(last_window) / len(last_window)) if last_window else None,
        "trend_slope": slope(rows),
        "n": float(len(vals)),
    }


def put_metric(row: dict[str, object], prefix: str, item: dict[str, float | None]) -> None:
    for key in ("mean", "median", "min", "max", "p95", "last", "last2ns_mean", "trend_slope", "n"):
        row[f"{prefix}_{key}"] = item.get(key)


def leq(value: float | None, limit: float) -> bool:
    return value is not None and value <= limit


def geq(value: float | None, limit: float) -> bool:
    return value is not None and value >= limit


def decide(interval: dict[str, dict[str, float | None]]) -> tuple[str, bool, str, str]:
    complex_ok = leq(interval["rmsd_complex"]["last2ns_mean"], 3.0)
    antibody_ok = leq(interval["rmsd_antibody"]["last2ns_mean"], 0.5)
    antigen_ok = leq(interval["rmsd_antigen"]["last2ns_mean"], 1.5)
    contacts_ok = geq(interval["interface_contacts"]["mean"], 150.0)
    mindist_ok = leq(interval["interface_mindist"]["mean"], 0.20)
    buried_ok = geq(interval["buried_sasa"]["mean"], 15.0)
    hbond_ok = geq(interval["hbonds"]["mean"], 5.0)

    core_failures: list[str] = []
    if not complex_ok:
        core_failures.append("complex RMSD last2ns mean > 3.0 nm")
    if not antibody_ok:
        core_failures.append("antibody RMSD last2ns mean > 0.5 nm")
    if not antigen_ok:
        core_failures.append("antigen RMSD last2ns mean > 1.5 nm")
    if not contacts_ok:
        core_failures.append("interface contacts mean < 150")
    if not mindist_ok:
        core_failures.append("interface min distance mean > 0.20 nm")
    if not buried_ok:
        core_failures.append("buried SASA mean < 15 nm^2")

    hbond_note = "H-bonds support stability" if hbond_ok else "H-bonds low; treated as supporting evidence only"
    if core_failures:
        if not hbond_ok and (not contacts_ok or not buried_ok):
            hbond_note = "H-bonds low together with ICs/SASA weakness"
        return "fail", False, "; ".join(core_failures), hbond_note
    if not hbond_ok:
        return "needs_review", True, "core interface/RMSD metrics stable but H-bonds are low", hbond_note
    return "pass", True, "core interface/RMSD metrics stable", hbond_note


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def append_csv(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and path.stat().st_size > 0
    fieldnames = list(row.keys())
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-id", required=True)
    parser.add_argument("--stage-ns", type=int, required=True)
    parser.add_argument("--interval-start-ns", type=float, required=True)
    parser.add_argument("--interval-end-ns", type=float, required=True)
    parser.add_argument("--interval-analysis-dir", required=True)
    parser.add_argument("--cumulative-analysis-dir", required=True)
    parser.add_argument("--stage-summary", required=True)
    parser.add_argument("--partial-results", required=True)
    parser.add_argument("--all-results", required=True)
    parser.add_argument("--json-summary", required=True)
    args = parser.parse_args()

    rows_by_scope: dict[str, dict[str, dict[str, float | None]]] = {}
    for scope, directory in {
        "interval": Path(args.interval_analysis_dir),
        "cumulative": Path(args.cumulative_analysis_dir),
    }.items():
        rows_by_scope[scope] = {}
        for metric, filename in METRIC_FILES.items():
            rows_by_scope[scope][metric] = stats(read_xy(directory / filename))

    stage_call, continue_next, reason, hbond_note = decide(rows_by_scope["interval"])
    row: dict[str, object] = {
        "candidate_id": args.candidate_id,
        "stage_ns": args.stage_ns,
        "interval_start_ns": args.interval_start_ns,
        "interval_end_ns": args.interval_end_ns,
        "stage_call": stage_call,
        "continue_next_stage": str(continue_next).lower(),
        "decision_reason": reason,
        "hbond_supporting_note": hbond_note,
    }
    for scope, metrics in rows_by_scope.items():
        for metric, item in metrics.items():
            put_metric(row, f"{scope}_{metric}", item)

    write_csv(Path(args.stage_summary), [row])
    append_csv(Path(args.partial_results), row)
    append_csv(Path(args.all_results), row)
    Path(args.json_summary).write_text(
        json.dumps({"row": row, "stats": rows_by_scope}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps({"stage_call": stage_call, "continue_next_stage": continue_next, "reason": reason}, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
