#!/usr/bin/env python3
"""Collect binder generation outputs into CSV/JSON tables.

The script is intentionally independent from PPIFlow internals: if
sample_metrics.csv exists it preserves those metrics, and it also computes
basic chain/contact checks directly from the PDB files.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path
from typing import Iterable


HOTSPOTS_DEFAULT = "A87,A88,A89,A90,A91,A114,A115,A116,A117"


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def parse_float(value: object, default: float = math.nan) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return default


def parse_hotspots(text: str) -> list[tuple[str, int]]:
    hotspots: list[tuple[str, int]] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        hotspots.append((item[0], int(item[1:])))
    return hotspots


def parse_pdb_atoms(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    with path.open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom = line[12:16].strip()
            if atom.startswith("H"):
                continue
            try:
                atoms.append(
                    {
                        "atom": atom,
                        "resn": line[17:20].strip(),
                        "chain": line[21].strip(),
                        "resi": int(line[22:26]),
                        "icode": line[26].strip(),
                        "xyz": (
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                        ),
                    }
                )
            except ValueError:
                continue
    return atoms


def squared_distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def residue_count(atoms: Iterable[dict[str, object]], chain: str) -> int:
    return len({(atom["resi"], atom["icode"]) for atom in atoms if atom["chain"] == chain})


def chain_ids(atoms: Iterable[dict[str, object]]) -> str:
    return "".join(sorted({str(atom["chain"]) for atom in atoms if atom["chain"]}))


def min_distance_between(a_atoms: list[dict[str, object]], b_atoms: list[dict[str, object]]) -> float:
    if not a_atoms or not b_atoms:
        return math.nan
    best = math.inf
    for atom_a in a_atoms:
        xyz_a = atom_a["xyz"]
        for atom_b in b_atoms:
            d2 = squared_distance(xyz_a, atom_b["xyz"])
            if d2 < best:
                best = d2
    return math.sqrt(best)


def hotspot_contacts(
    atoms: list[dict[str, object]],
    hotspots: list[tuple[str, int]],
    binder_chain: str,
    cutoff: float,
) -> tuple[int, list[str], float]:
    binder_atoms = [atom for atom in atoms if atom["chain"] == binder_chain]
    cutoff2 = cutoff * cutoff
    contacted: list[str] = []
    max_min_distance = 0.0

    for chain, resi in hotspots:
        hotspot_atoms = [atom for atom in atoms if atom["chain"] == chain and atom["resi"] == resi]
        if not hotspot_atoms or not binder_atoms:
            continue
        best = math.inf
        for hotspot_atom in hotspot_atoms:
            xyz_h = hotspot_atom["xyz"]
            for binder_atom in binder_atoms:
                d2 = squared_distance(xyz_h, binder_atom["xyz"])
                if d2 < best:
                    best = d2
        min_dist = math.sqrt(best)
        max_min_distance = max(max_min_distance, min_dist)
        if best <= cutoff2:
            contacted.append(f"{chain}{resi}")
    return len(contacted), contacted, max_min_distance


def interchain_contact_count(
    atoms: list[dict[str, object]],
    target_chain: str,
    binder_chain: str,
    cutoff: float,
) -> int:
    target_atoms = [atom for atom in atoms if atom["chain"] == target_chain]
    binder_atoms = [atom for atom in atoms if atom["chain"] == binder_chain]
    cutoff2 = cutoff * cutoff
    count = 0
    for target_atom in target_atoms:
        xyz_t = target_atom["xyz"]
        for binder_atom in binder_atoms:
            if squared_distance(xyz_t, binder_atom["xyz"]) <= cutoff2:
                count += 1
    return count


def find_metric_rows(out_root: Path, include_smoke: bool) -> dict[Path, dict[str, str]]:
    rows_by_pdb: dict[Path, dict[str, str]] = {}
    metric_files = sorted((out_root / "binder").glob("batch_*/sample_metrics.csv"))
    if include_smoke:
        metric_files.extend(sorted((out_root / "smoke_test").glob("sample_metrics.csv")))

    for metric_file in metric_files:
        with metric_file.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                sample = row.get("sample") or row.get("pdb") or row.get("pdb_file") or ""
                pdb_path = Path(row.get("pdb_path") or "")
                if not str(pdb_path):
                    pdb_path = metric_file.parent / sample
                if not pdb_path.is_absolute():
                    pdb_path = (metric_file.parent / pdb_path).resolve()
                rows_by_pdb[pdb_path] = dict(row)
    return rows_by_pdb


def discover_pdbs(out_root: Path, include_smoke: bool) -> list[Path]:
    pdbs = sorted((out_root / "binder").glob("batch_*/*.pdb"))
    if include_smoke:
        pdbs.extend(sorted((out_root / "smoke_test").glob("*.pdb")))
    return pdbs


def first_pass_score(row: dict[str, str]) -> float:
    coverage = parse_float(row.get("hotspot_coverage_6a", row.get("hotspot_coverage", "0")), 0.0)
    contact_count = parse_float(row.get("interchain_contact_count_5a", "0"), 0.0)
    max_hotspot_min = parse_float(row.get("max_hotspot_min_distance", "20"), 20.0)
    score = 100.0 * coverage + 0.02 * contact_count - max_hotspot_min
    if parse_bool(row.get("has_clash", "False")):
        score -= 100.0
    return score


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_checkpoint(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-root", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--expected-count", type=int, default=300)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--target-chain", default="A")
    parser.add_argument("--binder-chain", default="T")
    parser.add_argument("--hotspots", default=HOTSPOTS_DEFAULT)
    parser.add_argument("--include-smoke", action="store_true")
    args = parser.parse_args()

    out_root = Path(args.out_root)
    results_dir = Path(args.results_dir)
    all_results = results_dir / "binder_all_results.csv"
    partial_top = results_dir / "binder_partial_top.csv"
    summary_path = results_dir / "binder_summary.json"
    checkpoint = results_dir / "15_binder_collect_checkpoint.json"

    print(f"[start] out_root={out_root}", flush=True)
    print(f"[start] results_dir={results_dir}", flush=True)
    print(f"[start] expected_count={args.expected_count}", flush=True)

    hotspots = parse_hotspots(args.hotspots)
    metrics_by_pdb = find_metric_rows(out_root, args.include_smoke)
    pdbs = discover_pdbs(out_root, args.include_smoke)
    total = len(pdbs)
    rows: list[dict[str, str]] = []
    start = time.time()

    write_checkpoint(
        checkpoint,
        {
            "script": "15_collect_binder_metrics.py",
            "status": "running",
            "done": 0,
            "total": total,
            "out_root": str(out_root),
            "updated_at_epoch": int(time.time()),
        },
    )

    for index, pdb_path in enumerate(pdbs, start=1):
        elapsed = time.time() - start
        percent = 100.0 * index / max(total, 1)
        print(
            f"[collect-binder] {index}/{total} {percent:.1f}% elapsed={elapsed:.1f}s pdb={pdb_path.name}",
            flush=True,
        )

        atoms = parse_pdb_atoms(pdb_path)
        target_atoms = [atom for atom in atoms if atom["chain"] == args.target_chain]
        binder_atoms = [atom for atom in atoms if atom["chain"] == args.binder_chain]
        contact_5a, contacted_5a, _ = hotspot_contacts(atoms, hotspots, args.binder_chain, 5.0)
        contact_6a, contacted_6a, max_hotspot_min = hotspot_contacts(atoms, hotspots, args.binder_chain, 6.0)
        min_interchain = min_distance_between(target_atoms, binder_atoms)

        metric_row = metrics_by_pdb.get(pdb_path.resolve(), metrics_by_pdb.get(pdb_path, {}))
        batch_id = pdb_path.parent.name
        design_id = pdb_path.stem
        hotspot_total = len(hotspots)
        row: dict[str, str] = {
            "rank": "",
            "design_id": design_id,
            "sample": pdb_path.name,
            "modality": "binder",
            "batch_id": batch_id,
            "pdb_path": str(pdb_path),
            "target_chain": args.target_chain,
            "binder_chain": args.binder_chain,
            "chains": chain_ids(atoms),
            "binder_length": str(residue_count(atoms, args.binder_chain)),
            "target_length": str(residue_count(atoms, args.target_chain)),
            "hotspots": args.hotspots,
            "hotspot_contacted_5a": str(contact_5a),
            "hotspot_contacted_6a": str(contact_6a),
            "hotspot_total": str(hotspot_total),
            "hotspot_coverage_5a": f"{contact_5a / hotspot_total:.3f}",
            "hotspot_coverage_6a": f"{contact_6a / hotspot_total:.3f}",
            "coverage_hotspot_list_5a": ",".join(contacted_5a),
            "coverage_hotspot_list_6a": ",".join(contacted_6a),
            "max_hotspot_min_distance": f"{max_hotspot_min:.3f}",
            "min_interchain_distance": "" if math.isnan(min_interchain) else f"{min_interchain:.3f}",
            "interchain_contact_count_5a": str(interchain_contact_count(atoms, args.target_chain, args.binder_chain, 5.0)),
            "has_clash": str((not math.isnan(min_interchain)) and min_interchain < 2.0),
            "metric_hotspot_coverage": metric_row.get("hotspot_coverage", ""),
            "metric_dsasa": metric_row.get("dsasa", ""),
            "metric_has_clash": metric_row.get("has_clash", ""),
            "raw_metrics_json": json.dumps(metric_row, sort_keys=True),
            "status": "collected",
            "notes": "Binder structural metrics; PRODIGY not yet joined",
        }
        row["first_pass_score"] = f"{first_pass_score(row):.4f}"
        rows.append(row)

        if index % 10 == 0 or index == total:
            ranked = sorted(rows, key=lambda item: parse_float(item["first_pass_score"], -9999), reverse=True)
            for rank, item in enumerate(ranked, start=1):
                item["rank"] = str(rank)
            fields = list(ranked[0].keys()) if ranked else []
            write_csv(all_results, ranked, fields)
            write_csv(partial_top, ranked[: args.top_n], fields)
            write_checkpoint(
                checkpoint,
                {
                    "script": "15_collect_binder_metrics.py",
                    "status": "running",
                    "done": index,
                    "total": total,
                    "all_results": str(all_results),
                    "partial_top": str(partial_top),
                    "updated_at_epoch": int(time.time()),
                },
            )

    rows.sort(key=lambda item: parse_float(item["first_pass_score"], -9999), reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = str(rank)
    fields = list(rows[0].keys()) if rows else []
    write_csv(all_results, rows, fields)
    write_csv(partial_top, rows[: args.top_n], fields)

    clash_free = sum(1 for row in rows if not parse_bool(row.get("has_clash", "False")))
    full_hotspot = sum(1 for row in rows if parse_float(row.get("hotspot_coverage_6a", "0"), 0.0) >= 0.999)
    summary = {
        "script": "15_collect_binder_metrics.py",
        "status": "completed" if len(rows) >= args.expected_count else "partial",
        "out_root": str(out_root),
        "results_dir": str(results_dir),
        "expected_count": args.expected_count,
        "total_designs": len(rows),
        "clash_free": clash_free,
        "full_hotspot_coverage_6a": full_hotspot,
        "top_n": args.top_n,
        "all_results": str(all_results),
        "partial_top": str(partial_top),
        "updated_at_epoch": int(time.time()),
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_checkpoint(
        checkpoint,
        {
            **summary,
            "summary": str(summary_path),
        },
    )

    print(f"[done] total_designs={len(rows)} expected={args.expected_count}", flush=True)
    print(f"[done] all_results={all_results}", flush=True)
    print(f"[done] partial_top={partial_top}", flush=True)
    print(f"[done] summary={summary_path}", flush=True)
    return 0 if len(rows) >= args.expected_count else 1


if __name__ == "__main__":
    raise SystemExit(main())
