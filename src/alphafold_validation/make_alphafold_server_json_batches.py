#!/usr/bin/env python3
"""Create AlphaFold Server JSON upload files from the Task 4 queue."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TASK4_ROOT = ROOT / "results" / "task4_alphafold_server_validation_20260517_112024"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task4-root", type=Path, default=DEFAULT_TASK4_ROOT)
    parser.add_argument("--batch-size", type=int, default=30)
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = ""
    seq_parts: list[str] = []
    with path.open("r", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records.append((name, "".join(seq_parts)))
                name = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if name:
        records.append((name, "".join(seq_parts)))
    return records


def make_job(row: dict[str, str]) -> dict[str, object]:
    records = read_fasta(Path(row["fasta_path"]))
    if len(records) != 3:
        raise ValueError(f"{row['design_id']} expected 3 FASTA records, found {len(records)}")
    sequences = []
    for header, sequence in records:
        bad = sorted(set(sequence) - VALID_AA)
        if bad:
            raise ValueError(f"{row['design_id']} {header} contains unsupported amino acids: {''.join(bad)}")
        sequences.append({"proteinChain": {"sequence": sequence, "count": 1}})
    return {
        "name": f"task4_{int(row['submission_rank']):03d}_{row['design_id']}",
        "modelSeeds": [],
        "sequences": sequences,
        "dialect": "alphafoldserver",
        "version": 1,
    }


def main() -> int:
    args = parse_args()
    manifest = args.task4_root / "af3_server_job_manifest.csv"
    out_dir = args.task4_root / "server_inputs" / "alphafold_server_json"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(manifest.open("r", encoding="utf-8-sig", newline="")))
    jobs = [make_job(row) for row in rows]

    smoke = jobs[:3]
    smoke_path = out_dir / "task4_smoke_top3_alphafold_server_jobs.json"
    smoke_path.write_text(json.dumps(smoke, indent=2), encoding="ascii")

    full_paths = []
    for start in range(0, len(jobs), args.batch_size):
        batch = jobs[start : start + args.batch_size]
        batch_no = start // args.batch_size + 1
        path = out_dir / f"task4_batch_{batch_no:02d}_jobs_{start + 1:03d}_{start + len(batch):03d}.json"
        path.write_text(json.dumps(batch, indent=2), encoding="ascii")
        full_paths.append(path)

    all_path = out_dir / "task4_all_150_alphafold_server_jobs.json"
    all_path.write_text(json.dumps(jobs, indent=2), encoding="ascii")

    summary = {
        "status": "done",
        "total_jobs": len(jobs),
        "smoke_json": str(smoke_path),
        "batch_size": args.batch_size,
        "batch_json_files": [str(path) for path in full_paths],
        "all_json": str(all_path),
        "dialect": "alphafoldserver",
        "modelSeeds": [],
        "note": "AlphaFold Server JSON format is different from local AlphaFold3 JSON; this file uses proteinChain entities and a list of jobs.",
    }
    (out_dir / "alphafold_server_json_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"[DONE] server_json_jobs={len(jobs)}", flush=True)
    print(f"[DONE] smoke_json={smoke_path}", flush=True)
    for path in full_paths:
        print(f"[DONE] batch_json={path}", flush=True)
    print(f"[DONE] all_json={all_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
